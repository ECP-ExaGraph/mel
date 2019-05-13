/* Maximal edge matching using graph communicator and RMA */

#pragma once
#ifndef MAXEMATCHNRM_HPP
#define MAXEMATCHNRM_HPP

#include "graph.hpp"

#include <numeric>
#include <cstring>
#include <cassert>

#define MATE_REQUEST_TAG        (1) // mate[x] = y, if mate[y] = x, then match (x/y in different processes)
#define MATE_REJECT_TAG         (2) // reject vertex mate, requires to update mate
#define MATE_INVALID_TAG        (3) // invalidate edge

class MaxEdgeMatchNRM
{
    public:
        MaxEdgeMatchNRM(Graph* g): 
            g_(g), D_(0), M_(0), 
            sources_(0), targets_(0), sendbuf_(0),
            nghosts_in_target_(0), nghosts_target_indices_(0), 
            pindex_(0), indegree_(-1), outdegree_(-1), 
            prcounts_(0), scounts_(0), rcounts_(0), rdispls_(0), 
            nwin_(MPI_WIN_NULL), winbuf_(nullptr)
        {
            const GraphElem lnv = g_->get_lnv();
            comm_ = g_->get_comm();

            MPI_Comm_size(comm_, &size_);
            MPI_Comm_rank(comm_, &rank_);
             
            // create communication graph, 
            // update in|out degrees
            create_graph_topo();
                 
               // cache index corresponding to a process
            for (int i = 0; i < outdegree_; i++)
                pindex_.insert({targets_[i], (GraphElem)i}); 
            
            nghosts_in_target_.resize(outdegree_);
            nghosts_target_indices_.resize(outdegree_);          
            rdispls_.resize(indegree_);

            mate_.resize(lnv);
            std::fill(mate_.begin(), mate_.end(), -1);
             
            // populate counter that tracks number
            // of ghosts not owned by me
            GraphElem tot_ghosts = 0;

            for (GraphElem i = 0; i < lnv; i++)
            {
                GraphElem e0, e1;
                g_->edge_range(i, e0, e1);

                for (GraphElem e = e0; e < e1; e++)
                {
                    Edge const& edge = g_->get_edge(e);
                    
                    // all edge is active in the beginning,
                    // so no need to check edge.active_
                    const int target = g_->get_owner(edge.tail_);
                    if (target != rank_)
                    {
                        nghosts_in_target_[pindex_[target]] += 1;
                        tot_ghosts += 1;
                    }
                }                
            }
            
            // initialize input buffer
            
            // sends a pair of vertices with tag,
            // can send at most 2 messages along a
            // cross edge
            sendbuf_ = new GraphElem[tot_ghosts*3*2];

            // allocate MPI windows and lock them
                      
            // TODO FIXME make the same changes for
            // the base RMA version
            MPI_Info info = MPI_INFO_NULL;

#if defined(USE_MPI_ACCUMULATE)
            MPI_Info_create(&info);
            MPI_Info_set(info, "accumulate_ordering", "none");
            MPI_Info_set(info, "accumulate_ops", "same_op");
#endif

            // TODO FIXME report bug, program crashes when g_comm_ used
            MPI_Win_allocate((tot_ghosts*3*2)*sizeof(GraphElem), 
                    sizeof(GraphElem), info, comm_, &winbuf_, &nwin_);             
 
            MPI_Win_lock_all(MPI_MODE_NOCHECK, nwin_);

            // exclusive scan to compute remote 
            // displacements for RMA CALLS
            GraphElem disp = 0;
            for (int t = 0; t < outdegree_; t++)
            {
                nghosts_target_indices_[t] = disp;
                disp += nghosts_in_target_[t]*3*2;
            }
               
            // incoming updated prefix sums
            MPI_Neighbor_alltoall(nghosts_target_indices_.data(), 1, 
                    MPI_GRAPH_TYPE, rdispls_.data(), 1, MPI_GRAPH_TYPE, g_comm_);

            // set neighbor alltoall params
            scounts_.resize(outdegree_, 0);
            rcounts_.resize(indegree_, 0);
            prcounts_.resize(indegree_, 0);
        }

        ~MaxEdgeMatchNRM() {}

        void clear()
        {
            M_.clear();
            mate_.clear();
           
            sources_.clear();
            targets_.clear();

            scounts_.clear();
            rcounts_.clear();
            prcounts_.clear();
            rdispls_.clear();

            nghosts_in_target_.clear();
            nghosts_target_indices_.clear();
            pindex_.clear();
            
            delete []sendbuf_;
            
            // clear windows
            MPI_Win_unlock_all(nwin_);
            MPI_Win_free(&nwin_);

            MPI_Comm_free(&g_comm_);
        }

        // sources
        // indegree -- number of processes for which 
        // the calling process is the destination

        // destinations
        // outdegree -- number of processes for which calling 
        // process is source
        void create_graph_topo()
        {
            const GraphElem lnv = g_->get_lnv();
            for (GraphElem v = 0; v < lnv; v++)
            {
                GraphElem e0, e1;
                g_->edge_range(v, e0, e1);

                for (GraphElem e = e0; e < e1; e++)
                {
                    Edge const& edge = g_->get_edge(e);
                    const int owner = g_->get_owner(edge.tail_); 
                    if (owner != rank_)
                    {
                        // graph topology assumes directed
                        // graph, so edges stored twice
                        if (std::find(targets_.begin(), targets_.end(), owner) 
                                == targets_.end()
                                && std::find(sources_.begin(), sources_.end(), owner) 
                                == sources_.end())
                        {
                            targets_.push_back(owner);
                            sources_.push_back(owner);
                        }
                    }
                }
            }

            // TODO FIXME consider edge weights equivalent to number of 
            // ghost vertices shared between two processes
            // MPI Spec v3.1 p297: "Multiplicity of edges can likewise indicate 
            // more intense communication between pairs of processes."
            MPI_Dist_graph_create_adjacent(comm_, sources_.size(), sources_.data(), 
                    MPI_UNWEIGHTED, targets_.size(), targets_.data(), MPI_UNWEIGHTED, 
                    MPI_INFO_NULL, 0 /*reorder ranks?*/, &g_comm_);

            // indegree/outdegree
            // No need another MPI function call (just size of sources/targets 
            // would do), but just checking...
            int weighted;
            MPI_Dist_graph_neighbors_count(g_comm_, &indegree_, &outdegree_, &weighted);
            
            // to get sources/targets, use MPI_Dist_graph_neighbors
            assert(indegree_ == sources_.size());
            assert(outdegree_ == targets_.size());
        }
                
        /* Validation */
        // if mate[mate[v]] == v then
        // we're good
        void check_results()
        {
            // gather M_ and mate_
            const int lnv = g_->get_lnv();
            unsigned int m_size = M_.size(), m_global_size = 0;
            // i,j
            m_size *= 2;
            GraphElem* M_buf = new GraphElem[m_size];

            GraphElem* M_global = nullptr;
            GraphElem* mate_global = nullptr;
            
            // communication params from M_ and mate_
            int* rcounts = nullptr;
            int* rdispls = nullptr;
            int* m_rcounts = nullptr;
            int* m_rdispls = nullptr;
            
            // communication params for M
            if (rank_ == 0)
            {
                rcounts = new int[size_];
                rdispls = new int[size_];
                m_rcounts = new int[size_];
                m_rdispls = new int[size_];
            }

            // put M_ into a contiguous buffer
            for (int i = 0, j = 0; i < m_size; i+=2, j++)
            {
                M_buf[i]    = M_[j].ij_[0];
                M_buf[i+1]  = M_[j].ij_[1];
            }

            MPI_Gather(&m_size, 1, MPI_INT, rcounts, 1, MPI_INT, 0, comm_);
            MPI_Gather(&lnv, 1, MPI_INT, m_rcounts, 1, MPI_INT, 0, comm_);
            MPI_Reduce(&m_size, &m_global_size, 1, MPI_INT, MPI_SUM, 0, comm_);
            
            // communication params (at root)
            if (rank_ == 0)
            {
                const GraphElem nv = g_->get_nv();
                mate_global = new GraphElem[nv];
                M_global = new GraphElem[m_global_size];

                unsigned int index = 0, m_index = 0;
                for (int p = 0; p < size_; p++)
                {
                    rdispls[p] = index;
                    index += rcounts[p];
                    m_rdispls[p] = m_index;
                    m_index += m_rcounts[p];
                }
            }
            
            MPI_Barrier(comm_);

            // M_
            MPI_Gatherv(M_buf, m_size, MPI_GRAPH_TYPE, M_global, rcounts, rdispls, 
                    MPI_GRAPH_TYPE, 0, comm_);
            // mate
            MPI_Gatherv(mate_.data(), lnv, MPI_LONG, mate_global, m_rcounts, m_rdispls, 
                    MPI_GRAPH_TYPE, 0, comm_);
            
            MPI_Barrier(comm_);

            // data gathered, now validate
            if (rank_ == 0)
            {
                bool success = true;
                for (int i = 0; i < m_global_size; i+=2)
                {
                    if ((mate_global[mate_global[M_global[i]]] != M_global[i])
                            || (mate_global[mate_global[M_global[i+1]]] != M_global[i+1]))
                    {
                        std::cout << "Validation FAILED." << std::endl; 
                        std::cout << "mate_[mate_[" << M_global[i] << "]] != " << M_global[i] << " OR " 
                            << "mate_[mate_[" << M_global[i+1] << "]] != " << M_global[i+1] << std::endl;
                        success = false;
                        break;
                    }
                }
                if (success) 
                    std::cout << "Validation SUCCESS." << std::endl;
            }

            // clear buffers
            delete []M_global;
            delete []mate_global;
            delete []M_buf;

            delete []rcounts;
            delete []rdispls;
            delete []m_rcounts;
            delete []m_rdispls;
        }

        // print the contents of M_
        void print_M() const
        {
            // gather M_
            unsigned int m_size = M_.size(), m_global_size = 0;
            // i,j
            m_size *= 2;
            GraphElem* M_buf = new GraphElem[m_size];

            GraphElem* M_global = nullptr;
            int* rcounts = nullptr;
            int* rdispls = nullptr;

            // communication params
            if (rank_ == 0)
            {
                rcounts = new int[size_];
                rdispls = new int[size_];
            }

            // put M_ into a contiguous buffer
            for (int i = 0, j = 0; i < m_size; i+=2, j++)
            {
                M_buf[i]    = M_[j].ij_[0];
                M_buf[i+1]  = M_[j].ij_[1];
            }

            MPI_Gather(&m_size, 1, MPI_INT, rcounts, 1, MPI_INT, 0, comm_);
            MPI_Reduce(&m_size, &m_global_size, 1, MPI_INT, MPI_SUM, 0, comm_);

            // communication params (at root)
            if (rank_ == 0)
            {
                M_global = new GraphElem[m_global_size];

                unsigned int index = 0;
                for (int p = 0; p < size_; p++)
                {
                    rdispls[p] = index;
                    index += rcounts[p];
                }
            }

            MPI_Gatherv(M_buf, m_size, MPI_GRAPH_TYPE, M_global, rcounts, rdispls, 
                    MPI_GRAPH_TYPE, 0, comm_);
            MPI_Barrier(comm_);

            // print mates
            if (rank_ == 0)
            {
                std::cout << "Matched vertices: " << std::endl;
                for (int i = 0; i < m_global_size; i+=2)
                    std::cout << M_global[i] << " ---- " << M_global[i+1] << std::endl;
            }

            // clear buffers
            delete []M_global;
            delete []M_buf;
            delete []rcounts;
            delete []rdispls;
        }
                 
        // TODO FIXME not expecting a, b to
        // be large, if large then following
        // absolute tolerance test will fail:
        // http://realtimecollisiondetection.net/blog/?p=89
        bool is_same(double a, double b) 
        { return std::abs(a - b) <= std::numeric_limits<double>::epsilon(); }
 
        // expecting v to be local index
        // require global_to_local
        // before passing 
        // local computation
        void compute_mate(GraphElem v, Edge& max_edge)
        {
            GraphElem e0, e1;
            g_->edge_range(v, e0, e1);

            for (GraphElem e = e0; e < e1; e++)
            {
                EdgeActive& edge = g_->get_active_edge(e);
                if (edge.active_)
                {
                    if (edge.edge_.weight_ > max_edge.weight_)
                        max_edge = edge.edge_;

                    // break tie using vertex index
                    if (is_same(edge.edge_.weight_, max_edge.weight_))
                        if (edge.edge_.tail_ > max_edge.tail_)
                            max_edge = edge.edge_;
                }
            }
        }
        
        // maximal edge matching
        std::vector<EdgeTuple> const& operator()()
        {
            maxematch_nrm();
            return M_;
        }
         
        // search v in M_ (local list
        // of matched vertices)
        bool is_matched(GraphElem v)
        {
            auto found = std::find_if(M_.begin(), M_.end(), 
                    [&](EdgeTuple const& et) 
                    { return ((et.ij_[0] == v) || (et.ij_[1] == v)); });
            if (found == std::end(M_))
                return false;
            return true;
        }   

               
        // x is owned by me, y may be a ghost
        // deactivate edge x -- y and decrement
        inline void deactivate_edge(GraphElem x, GraphElem y)
        {
            GraphElem e0, e1;
            const GraphElem lx = g_->global_to_local(x);
            const int y_owner = g_->get_owner(y);
            const GraphElem pidx = pindex_[y_owner];

            g_->edge_range(lx, e0, e1);

            for (GraphElem e = e0; e < e1; e++)
            {
                EdgeActive& edge = g_->get_active_edge(e);
                if (edge.edge_.tail_ == y && edge.active_)
                {
                    edge.active_ = false;

                    if (y_owner != rank_)
                        nghosts_in_target_[pidx] -= 1;
                    
                    break;
                }
            }
        }

        // initiate put
        void Put(int tag, int target, GraphElem data[2])
        {
            const int pidx = pindex_[target];
            const GraphElem curr_count = scounts_[pidx];
            const GraphElem index = nghosts_target_indices_[pidx] + curr_count;

            sendbuf_[index] = data[0];
            sendbuf_[index + 1] = data[1];
            sendbuf_[index + 2] = tag;

            // get displacement
            GraphElem tdisp = rdispls_[pidx] + curr_count;
                    
            // TODO FIXME consider type contiguos or struct
            // to pack these items?

#if defined(USE_MPI_ACCUMULATE)
            MPI_Accumulate(&sendbuf_[index], 3, MPI_GRAPH_TYPE, target, 
                    (MPI_Aint)tdisp, 3, MPI_GRAPH_TYPE, MPI_REPLACE, nwin_);
#else
            MPI_Put(&sendbuf_[index], 3, MPI_GRAPH_TYPE, target, 
                    (MPI_Aint)tdisp, 3, MPI_GRAPH_TYPE, nwin_);
#endif
            scounts_[pidx] += 3;
        }

        // x is owned by me
        // compute y = mate[x], if mate[y] = x, match
        // else if y = -1, invalidate all edges adj(x)
        void find_mate(GraphElem x)
        {
            const GraphElem lx = g_->global_to_local(x);
            Edge x_max_edge;

            compute_mate(lx, x_max_edge);
            const GraphElem y = mate_[lx] = x_max_edge.tail_;

            // initiate matching request
            if (y != -1)
            {
                // check if y can be matched
                const int y_owner = g_->get_owner(y);
                if (y_owner == rank_)
                {
                    if (mate_[g_->global_to_local(y)] == x)
                    {
                        D_.push_back(x);
                        D_.push_back(y);
                        M_.emplace_back(x, y, x_max_edge.weight_);

                        // mark x-y inactive, because its matched
                        deactivate_edge(x, y);
                        deactivate_edge(y, x);
                    }
                }
                else // send REQUEST
                {
                    deactivate_edge(x, y);

                    GraphElem data[2] = {y, x};
                    Put(MATE_REQUEST_TAG, y_owner, data);
                }
            }
            else // mate[x] = -1, deactivate all x - adj(x) edges
            {
                GraphElem e0, e1;
                g_->edge_range(lx, e0, e1);

                for (GraphElem e = e0; e < e1; e++)
                {
                    EdgeActive& edge = g_->get_active_edge(e);
                    
                    // deactivate only if edge is active
                    if (edge.active_) 
                    {   
                        edge.active_ = false;
                        const GraphElem z = edge.edge_.tail_;
                        const int z_owner = g_->get_owner(z);
                        
                        if (z_owner == rank_) 
                            deactivate_edge(z, x); // z - x
                        else // send INVALID (z - x) 
                        {
                            nghosts_in_target_[pindex_[z_owner]] -= 1;

                            GraphElem data[2] = {z, x};
                            Put(MATE_INVALID_TAG, z_owner, data);
                        }
                    }
                }
            }
        }

        // process matched vertices
        // in Part #2
        void process_neighbors(GraphElem v)
        {
            GraphElem e0, e1;
            const GraphElem lv = g_->global_to_local(v);
            g_->edge_range(lv, e0, e1);

            // find unmatched vertices
            // in v's neighborhood
            for (GraphElem e = e0; e < e1; e++)
            {
                EdgeActive& edge = g_->get_active_edge(e);

                if (edge.active_)
                {
                    const GraphElem x = edge.edge_.tail_;

                    if (mate_[lv] != x)
                    {
                        // invalidate v - x, because v
                        // is already matched, and not 
                        // with x
                        edge.active_ = false;
                        const int x_owner = g_->get_owner(x);

                        // find another mate for x, as v 
                        // is already matched
                        if (x_owner == rank_)
                        {
                            // invalidate x - v
                            deactivate_edge(x, v);

                            // find new candidate
                            if (mate_[g_->global_to_local(x)] == v)
                                find_mate(x);
                        }
                        else // send REJECT to invalidate x-v and recompute mate[x]
                        {                                                        
                            nghosts_in_target_[pindex_[x_owner]] -= 1; 

                            GraphElem data[2] = {x, v};
                            Put(MATE_REJECT_TAG, x_owner, data);
                        }
                    }
                }
            }
        }

        // process remote operations, need to flush pending 
        // RMA ops before progressing
        // ------------------------------------------------        
        void process_messages()
        {
            GraphElem g_l[2];
             
            // access local window and process data
            MPI_Win_flush_all(nwin_);
            
            // incoming data sizes
            MPI_Neighbor_alltoall(scounts_.data(), 1, MPI_GRAPH_TYPE, 
                    rcounts_.data(), 1, MPI_GRAPH_TYPE, g_comm_);
            
            // indegree == outdegree in our case
            for (int k = 0; k < indegree_; k++)
            {
                const GraphElem index = nghosts_target_indices_[k];
                const int start = prcounts_[k];
                const int end = rcounts_[k];
                
                for (int i = start; i < end; i+=3)
                {
                    g_l[0] = winbuf_[index + i];
                    g_l[1] = winbuf_[index + i + 1];
                    int status = (int)(winbuf_[index + i + 2]);

                    // REQUEST: may result in a match
                    if (status == MATE_REQUEST_TAG) 
                    {
                        bool matched = false;

                        // check if y is already matched
                        if (!is_matched(g_l[0]))
                        { 
                            // deactivate edge
                            deactivate_edge(g_l[0], g_l[1]);

                            if (mate_[g_->global_to_local(g_l[0])] == g_l[1])
                            {
                                M_.emplace_back(g_l[0], g_l[1], 0.0);

                                D_.push_back(g_l[0]);
                                D_.push_back(g_l[1]);

                                matched = true;
                            }
                        } 

                        // send REJECT if matching not possible
                        if (!matched)
                        {
                            // deactivate edge
                            deactivate_edge(g_l[0], g_l[1]);

                            const int x_owner = g_->get_owner(g_l[1]);
                            GraphElem data[2] = {g_l[1], g_l[0]};
                            Put(MATE_REJECT_TAG, x_owner, data);
                        }                
                    } 
                    else if (status == MATE_REJECT_TAG)
                    {
                        deactivate_edge(g_l[0], g_l[1]);

                        // recalculate mate[x]
                        if (mate_[g_->global_to_local(g_l[0])] == g_l[1])
                            find_mate(g_l[0]);
                    }
                    else // INVALID: deactivate x -- v
                        deactivate_edge(g_l[0], g_l[1]);
                }
                
                // retain past recv counts
                prcounts_[k] = rcounts_[k];
            }
        }
        
        // maximal weight matching main 
        void maxematch_nrm()
        {           
            const GraphElem lnv = g_->get_lnv();

            /* Phase #1: Part #1 -- Process locally owned vertices */
            for (GraphElem i = 0; i < lnv; i++)
                find_mate(g_->local_to_global(i));

            /* Phase #1: Part #2 -- Handle remotely owned vertices */
            while(1)
            {
                process_messages();
                do_matching();

                // exit criteria
                // check if all cross edges have been processed
                GraphElem count = std::accumulate(nghosts_in_target_.begin(), 
                        nghosts_in_target_.end(), 0);

                //std::cout << "[" << rank_ << "] count: " << count << std::endl;
                MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_GRAPH_TYPE, MPI_SUM, comm_);

                if (count == 0)
                    break;
            } // end of while(D_)
        }

        // locally process matched vertices
        // ignore ghost vertices
        void do_matching()
        {
            while (!D_.empty())
            {
                GraphElem v = D_.back();
                D_.pop_back();
                const int v_owner = g_->get_owner(v);
                
                if (v_owner == rank_) // check neighbors of v
                    process_neighbors(v);
            }
        }

    private:
        Graph* g_;
        std::vector<GraphElem> D_, mate_;
        std::vector<EdgeTuple> M_;
        
        // count of ghost vertices not owned by me
        // and counters
        std::unordered_map<int, GraphElem> 
            pindex_;                    // index of neighbor processes

        // intermediate communication buffers
        GraphElem* sendbuf_;
        
        std::vector<int> sources_, targets_;
        std::vector<GraphElem> scounts_, rcounts_, prcounts_;
        int indegree_, outdegree_;
    
        // RMA
        MPI_Win nwin_; // window to store data
        GraphElem* winbuf_;
        std::vector<GraphElem> rdispls_, // target displacement
            nghosts_in_target_,          // ghost vertices in target rank
            nghosts_target_indices_;     // indices of data 
        
        int rank_, size_;
        MPI_Comm comm_;
        MPI_Comm g_comm_; // neighborhood comm
};

#endif
