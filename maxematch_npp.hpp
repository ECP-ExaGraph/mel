/* Maximal edge matching using graph communicator and send/recv */

#pragma once
#ifndef MAXEMATCHNPP_HPP
#define MAXEMATCHNPP_HPP

#include "graph.hpp"

#include <numeric>
#include <cstring>
#include <cassert>

#define MATE_REQUEST_TAG        101 // mate[x] = y, if mate[y] = x, then match (x/y in different processes)
#define MATE_REJECT_TAG         102 // reject vertex mate, requires to update mate
#define MATE_INVALID_TAG        103 // invalidate edge

class MaxEdgeMatchNPP
{
    public:
        MaxEdgeMatchNPP(Graph* g): 
            g_(g), D_(0), M_(0), 
            sources_(0), targets_(0), 
            sendbuf_(0), recvbuf_(0),
            nghosts_in_target_(0), ghost_per_node_(0),
            nghosts_target_indices_(0), sendbuf_ctr_(0),
            sendreq_ctr_(0), indegree_(-1), outdegree_(-1)
        {
            const GraphElem lnv = g_->get_lnv();
            comm_ = g_->get_comm();

            MPI_Comm_size(comm_, &size_);
            MPI_Comm_rank(comm_, &rank_);
             
            // create communication graph, 
            // update in|out degrees
            create_graph_topo();
                       
            for (int i = 0; i < outdegree_; i++)
            {
                nghosts_in_target_.insert({targets_[i], 0});
                nghosts_target_indices_.insert({targets_[i], 0});
                sendbuf_ctr_.insert({targets_[i], 0});
            }
             
            ghost_per_node_.resize(lnv, 0);
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
                        nghosts_in_target_[target] += 1;
                        ghost_per_node_[i] += 1; // per vertex cross edges
                    }
                }                

                tot_ghosts += ghost_per_node_[i];
            }
            
            // initialize input buffer
            
            // sends a pair of vertices with tag,
            // can send at most 2 messages along a
            // cross edge
           sendbuf_ = new GraphElem[tot_ghosts*2*2];
           sendreq_ = new MPI_Request[tot_ghosts*2];
                      
            // prefix sum for calculating 
            // indices of outgoing buffers
            GraphElem disp = 0;
            for (int t = 0; t < outdegree_; t++)
            {
                nghosts_target_indices_[targets_[t]] = disp;
                disp += (nghosts_in_target_[targets_[t]]*2*2); 
            }
        }

        ~MaxEdgeMatchNPP() {}

        void clear()
        {
            M_.clear();
            mate_.clear();
           
            sources_.clear();
            targets_.clear();

            nghosts_in_target_.clear();
            ghost_per_node_.clear();
            nghosts_target_indices_.clear();

            sendbuf_ctr_.clear();
            recvbuf_.clear();
            
            delete []sendbuf_;
            delete []sendreq_;

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
            maxematch_npp();
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

            g_->edge_range(lx, e0, e1);

            for (GraphElem e = e0; e < e1; e++)
            {
                EdgeActive& edge = g_->get_active_edge(e);
                if (edge.edge_.tail_ == y && edge.active_)
                {
                    edge.active_ = false;

                    if (y_owner != rank_)
                        ghost_per_node_[lx] -= 1;
                    
                    break;
                }
            }
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

                    // FIXME TODO maintain a separate function
                    // for Isend like base P2P version, aka TaggedIsend
                    const GraphElem index = nghosts_target_indices_[y_owner] 
                        + sendbuf_ctr_[y_owner];

                    sendbuf_[index] = y;
                    sendbuf_[index + 1] = x;

                    // FIXME if ranks are reordered, then 
                    // this may fail
                    MPI_Isend(&sendbuf_[index], 2, MPI_GRAPH_TYPE, 
                            y_owner, MATE_REQUEST_TAG, g_comm_, 
                            &sendreq_[sendreq_ctr_]);

                    MPI_Request_free(&sendreq_[sendreq_ctr_]);
                    sendbuf_ctr_[y_owner] += 2;
                    sendreq_ctr_++;
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
                            ghost_per_node_[lx] -= 1;
                            
                            const GraphElem index = nghosts_target_indices_[z_owner] + sendbuf_ctr_[z_owner];
                            
                            sendbuf_[index] = z;
                            sendbuf_[index + 1] = x;
                    
                            MPI_Isend(&sendbuf_[index], 2, MPI_GRAPH_TYPE, 
                                    z_owner, MATE_INVALID_TAG, g_comm_, 
                                    &sendreq_[sendreq_ctr_]);

                            MPI_Request_free(&sendreq_[sendreq_ctr_]);
                            sendbuf_ctr_[z_owner] += 2;
                            sendreq_ctr_++;
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
                            ghost_per_node_[lv] -= 1; 

                            const GraphElem index = nghosts_target_indices_[x_owner] + sendbuf_ctr_[x_owner];

                            sendbuf_[index] = x;
                            sendbuf_[index + 1] = v;
                            
                            MPI_Isend(&sendbuf_[index], 2, MPI_GRAPH_TYPE, 
                                    x_owner, MATE_REJECT_TAG, g_comm_, 
                                    &sendreq_[sendreq_ctr_]);

                            MPI_Request_free(&sendreq_[sendreq_ctr_]);
                            sendbuf_ctr_[x_owner] += 2;
                            sendreq_ctr_++;
                        }
                    }
                }
            }
        }
                
        // remote operations, needs ncomm before progressing
        // ----------------------------------------------------        
        void process_messages()
        {
            GraphElem g_l[2];
            MPI_Status status;
            int flag = -1;           
            
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, g_comm_, 
                    &flag, &status);

            // post recv
            if (flag)
            { 
                MPI_Recv(g_l, 2, MPI_GRAPH_TYPE, status.MPI_SOURCE, 
                        status.MPI_TAG, g_comm_, &status);            
            }
            else
                return;
                
            // REQUEST: may result in a match
            if (status.MPI_TAG == MATE_REQUEST_TAG) 
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
                    const GraphElem index = nghosts_target_indices_[x_owner] + sendbuf_ctr_[x_owner];

                    sendbuf_[index] = g_l[1];
                    sendbuf_[index + 1] = g_l[0];
                    
                    MPI_Isend(&sendbuf_[index], 2, MPI_GRAPH_TYPE, 
                            x_owner, MATE_REJECT_TAG, g_comm_, 
                            &sendreq_[sendreq_ctr_]);

                    MPI_Request_free(&sendreq_[sendreq_ctr_]);
                    sendbuf_ctr_[x_owner] += 2;
                    sendreq_ctr_++;
                }                
            } 
            else if (status.MPI_TAG == MATE_REJECT_TAG)
            {
                deactivate_edge(g_l[0], g_l[1]);

                // recalculate mate[x]
                if (mate_[g_->global_to_local(g_l[0])] == g_l[1])
                    find_mate(g_l[0]);
            }
            else // INVALID: deactivate x -- v
                deactivate_edge(g_l[0], g_l[1]);
        }
        
        // maximal weight matching main 
        void maxematch_npp()
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
                GraphElem count = std::accumulate(ghost_per_node_.begin(), 
                        ghost_per_node_.end(), 0);
                //std::cout << "[" << rank_ << "] count: " << count << std::endl;

                if (count == 0)
                    break;
            } // end of while(D_)
            
            MPI_Waitall(sendreq_ctr_, sendreq_, MPI_STATUSES_IGNORE);
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
        std::vector<GraphElem> ghost_per_node_; 
        std::unordered_map<int, GraphElem> 
            nghosts_in_target_, // ghost vertices in target rank
            nghosts_target_indices_, // indices of data
            sendbuf_ctr_;

        // intermediate communication buffers
        GraphElem* sendbuf_;
        MPI_Request* sendreq_;
        GraphElem sendreq_ctr_;
        std::vector<GraphElem> recvbuf_;
        
        std::vector<int> sources_, targets_;
        int indegree_, outdegree_;
        
        int rank_, size_;
        MPI_Comm comm_;
        MPI_Comm g_comm_; // neighborhood comm
};

#endif
