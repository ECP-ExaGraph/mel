
/* Maximal edge matching using MPI Send/Recv */

#pragma once
#ifndef MAXEMATCHP2P_HPP
#define MAXEMATCHP2P_HPP

#include "graph.hpp"

#include <numeric>
#include <utility>
#include <cstring>

#include <omp.h>

// MPI message tags
#define MATE_REQUEST_TAG        1 // mate[x] = y, if mate[y] = x, then match (x/y in different processes)
#define MATE_REJECT_TAG         2 // reject vertex mate, requires to update mate
#define MATE_INVALID_TAG        3 // invalidate edge

class MaxEdgeMatchP2P
{
    public:
        MaxEdgeMatchP2P(Graph* g): 
            g_(g), D_(0), M_(0), 
            sbuf_ctr_(0), tot_ghosts_(0)
        {
            comm_ = g_->get_comm();
            MPI_Comm_size(comm_, &size_);
            MPI_Comm_rank(comm_, &rank_);

            // initialize mate_
            const GraphElem lnv = g_->get_lnv();
            mate_ = new GraphElem[lnv];
            std::fill(mate_, mate_ + lnv, -1);
                       
            // populate counter that tracks number
            // of ghosts not owned by me
            ghost_count_.resize(lnv, 0);

            for (GraphElem i = 0; i < lnv; i++)
            {
                GraphElem e0, e1;
                g_->edge_range(i, e0, e1);

                for (GraphElem e = e0; e < e1; e++)
                {
                    Edge const& edge = g_->get_edge(e);
                    
                    // all edge is active in the beginning,
                    // so no need to check edge.active_
                    if (g_->get_owner(edge.tail_) != rank_)
                        ghost_count_[i] += 1;
                }

                tot_ghosts_ += ghost_count_[i];
            }
            
            // each vertex can send at most
            // 2 messages along a cross edge
            sbuf_ = new GraphElem[tot_ghosts_*2*2]; // sending a pair
            sreq_ = new MPI_Request[tot_ghosts_*2];
        }

        ~MaxEdgeMatchP2P() {}

        void clear()
        {
            D_.clear();
            M_.clear();
            ghost_count_.clear();
            
            delete []mate_;
            delete []sbuf_;
            delete []sreq_;
        }
       
        // TODO FIXME not expecting a, b to
        // be large, if large then following
        // absolute tolerance test will fail:
        // http://realtimecollisiondetection.net/blog/?p=89
        inline bool is_same(GraphWeight a, GraphWeight b) 
        { return std::abs(a - b) <= std::numeric_limits<GraphWeight>::epsilon(); }

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
            MPI_Gatherv(mate_, lnv, MPI_GRAPH_TYPE, mate_global, m_rcounts, m_rdispls, 
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
            
            MPI_Barrier(comm_);

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
            
            MPI_Barrier(comm_);

            // clear buffers
            delete []M_global;
            delete []M_buf;
            delete []rcounts;
            delete []rdispls;
        }
                
        /* Maximal edge matching */
        std::vector<EdgeTuple> const& operator()()
        {
            maxematch_p2p();
            return M_;
        }

        // MPI_Isend
        void TaggedIsend(int tag, int target, GraphElem data[2])
        {
            // copy into persistent buffer before messaging 
            memcpy(&sbuf_[sbuf_ctr_], data, 2*sizeof(GraphElem));

            MPI_Isend(&sbuf_[sbuf_ctr_], 2, MPI_GRAPH_TYPE, 
                    target, tag, comm_, &sreq_[(sbuf_ctr_/2)]);

            MPI_Request_free(&sreq_[(sbuf_ctr_/2)]);

            sbuf_ctr_ += 2;
        }

        // search v in M_ (local list
        // of matched vertices)
        bool is_matched(const GraphElem v)
        {
#ifdef M_LIST_SORTED
#else
            auto found = std::find_if(M_.begin(), M_.end(), 
                    [&](EdgeTuple const& et) 
                    { return ((et.ij_[0] == v) || (et.ij_[1] == v)); });
            if (found == std::end(M_))
                return false;
#endif

            return true;
        }

        // process matched vertices
        // in my local queue, ignore
        // ghost vertices
        void do_matching()
        {
            while(!D_.empty())                               
            {    
                GraphElem v = D_.back();
                D_.pop_back();
                const int v_owner = g_->get_owner(v);

                if (v_owner == rank_)
                    process_neighbors(v);
            }
        }
        
        // process messages
        void process_messages()
        {
            // probe incoming msg
            MPI_Status status;
            int flag = -1;
            GraphElem g_l[2];
                           
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, 
                    &flag, &status);

            // post recv
            if (flag)
            { 
                MPI_Recv(g_l, 2, MPI_GRAPH_TYPE, status.MPI_SOURCE, 
                        status.MPI_TAG, comm_, MPI_STATUS_IGNORE);            
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

                    GraphElem pair[2] = {g_l[1], g_l[0]}; // {g,l}
                    TaggedIsend(MATE_REJECT_TAG, status.MPI_SOURCE, pair);
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
        
        // maximal edge matching using MPI P2P
        void maxematch_p2p()
        {
            // Part #1 -- initialize candidate mate
            GraphElem lnv = g_->get_lnv();
#if defined(USE_MPI_ONLY)
            for (GraphElem i = 0; i < lnv; i++)
                find_mate(g_->local_to_global(i));
#else
            // local computation (to be offloaded)
            std::vector<Edge> max_edges(lnv);
#ifdef OMP_TARGET_OFFLOAD
	    int ndevs = omp_get_num_devices();
	    int to_offload = (ndevs > 0);
            max_edges.resize(lnv);
            M_.resize(lnv);
            D_.resize(lnv*2);
            Edge* max_edges_ptr = max_edges.data();
            EdgeTuple* M_ptr = M_.data();
            GraphElem* D_ptr = D_.data();
            GraphElem* ghost_count_ptr = ghost_count_.data();
#pragma omp target parallel for if (to_offload) \
            map(to:g_) \
	    map(from:mate_[0:lnv], M_ptr[0:lnv], D_ptr[0:lnv*2], ghost_count_ptr[0:lnv], max_edges_ptr[0:lnv]) 
            for (GraphElem i = 0; i < lnv; i++)
            {
                GraphElem e0, e1;
                const GraphElem x = g_->local_to_global(i);
                const GraphElem lx = g_->global_to_local(x);
                g_->edge_range(lx, e0, e1);
                for (GraphElem e = e0; e < e1; e++)
                {
                    EdgeActive& edge = g_->get_active_edge(e);
                    if (edge.active_)
                    {
                        if (edge.edge_.weight_ > max_edges_ptr[i].weight_)
                            max_edges_ptr[i] = edge.edge_;
                        // break tie using vertex index
                        if (is_same(edge.edge_.weight_, max_edges_ptr[i].weight_))
                            if (edge.edge_.tail_ > max_edges_ptr[i].tail_)
                                max_edges_ptr[i] = edge.edge_;
                    }
                }
                #pragma omp atomic write
                mate_[lx] = max_edges_ptr[i].tail_;
                const GraphElem y = mate_[lx];
                // initiate matching request
                if (y != -1)
                {
                    // check if y can be matched
                    const int y_owner = g_->get_owner(y);
                    if (y_owner == rank_)
                    {
                        GraphElem mate_y;
                        #pragma omp atomic read
                        mate_y = mate_[g_->global_to_local(y)]; 
                        if (mate_y == x)
                        {
                            D_ptr[i    ] = x;
                            D_ptr[i + 1] = y;
                            EdgeTuple et(x, y, max_edges_ptr[i].weight_); 
                            M_ptr[i] = et;
			    // mark y<->x inactive, because its matched
                            deactivate_edge_device(y, x, ghost_count_ptr);
                            deactivate_edge_device(x, y, ghost_count_ptr);
                        }
                    }
                }
                else // invalidate all neigboring vertices 
                {
                    for (GraphElem e = e0; e < e1; e++)
                    {
                        EdgeActive& edge = g_->get_active_edge(e);
                        if (edge.active_)
                        {
                            const GraphElem z = edge.edge_.tail_;
                            const int z_owner = g_->get_owner(z);
                            if (z_owner == rank_)
                            {
                                edge.active_ = false;
                                deactivate_edge_device(z, x, ghost_count_ptr); // invalidate x -- z
                            }
                        }
                    }
                }
            }            
#else
#pragma omp declare reduction(merge : std::vector<GraphElem> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction(merge : std::vector<EdgeTuple> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp parallel for default(shared) reduction(merge: D_, M_) schedule(static)
	    for (GraphElem i = 0; i < lnv; i++)
            {
                GraphElem e0, e1;
                const GraphElem x = g_->local_to_global(i);
                const GraphElem lx = g_->global_to_local(x);
                g_->edge_range(lx, e0, e1);
                for (GraphElem e = e0; e < e1; e++)
                {
                    EdgeActive& edge = g_->get_active_edge(e);
                    if (edge.active_)
                    {
                        if (edge.edge_.weight_ > max_edges[i].weight_)
                            max_edges[i] = edge.edge_;

                        // break tie using vertex index
                        if (is_same(edge.edge_.weight_, max_edges[i].weight_))
                            if (edge.edge_.tail_ > max_edges[i].tail_)
                                max_edges[i] = edge.edge_;
                    }
                }
                #pragma omp atomic write
                mate_[lx] = max_edges[i].tail_;
                const GraphElem y = mate_[lx];
                
                // initiate matching request
                if (y != -1)
                {
                    // check if y can be matched
                    const int y_owner = g_->get_owner(y);
                    if (y_owner == rank_)
                    {
                        GraphElem mate_y;
                        #pragma omp atomic read
                        mate_y = mate_[g_->global_to_local(y)]; 
                        if (mate_y == x)
                        {
                            D_.push_back(x);
                            D_.push_back(y);
                            M_.emplace_back(x, y, max_edges[i].weight_);

                            // mark y<->x inactive, because its matched
                            deactivate_edge(y, x);
                            deactivate_edge(x, y);
                        }
                    }
                }
                else // invalidate all neigboring vertices 
                {
                    for (GraphElem e = e0; e < e1; e++)
                    {
                        EdgeActive& edge = g_->get_active_edge(e);
                        if (edge.active_)
                        {
                            const GraphElem z = edge.edge_.tail_;
                            const int z_owner = g_->get_owner(z);
                            if (z_owner == rank_)
                            {
                                edge.active_ = false;
                                deactivate_edge(z, x); // invalidate x -- z
                            }
                        }
                    }
                }
            }
#endif
            // OMP region ends
            // communication
            for (GraphElem i = 0; i < lnv; i++)
            {
                const GraphElem x = g_->local_to_global(i);
                const GraphElem lx = g_->global_to_local(x);
                const GraphElem y = mate_[lx] = max_edges[i].tail_;
                GraphElem e0, e1;
                g_->edge_range(lx, e0, e1);
                // initiate matching request
                if (y != -1)
                {
                    // check if y can be matched
                    const int y_owner = g_->get_owner(y);
                    if (y_owner != rank_) // ghost, send REQUEST
                    {
                        deactivate_edge(x, y);
                        GraphElem y_x[2] = {y, x};
                        TaggedIsend(MATE_REQUEST_TAG, y_owner, y_x);  
                    }
                }
                else // invalidate all neigboring vertices 
                {
                    for (GraphElem e = e0; e < e1; e++)
                    {
                        EdgeActive& edge = g_->get_active_edge(e);
                        if (edge.active_)
                        {
                            const GraphElem z = edge.edge_.tail_;
                            const int z_owner = g_->get_owner(z);
                            if (z_owner != rank_) // ghost, send INVALID
                            {
                                edge.active_ = false;
                                ghost_count_[lx] -= 1;
                                GraphElem z_x[2] = {z, x};
                                TaggedIsend(MATE_INVALID_TAG, z_owner, z_x); // invalidate x -- z
                            }
                        }
                    }
                }
            }
#endif            
            // Part 2 -- complete nb synch sends
            while(1)
            {                         
                process_messages();
                do_matching();
 
                // check if all cross edges have been processed
                GraphElem count = std::accumulate(ghost_count_.begin(), ghost_count_.end(), 0);
                if (count == 0)
                    break;
                //std::cout << "count: " << count << std::endl;
            }
        } 

        // expecting v to be local index
        // require global_to_local
        // before passing local computation
        void compute_mate(const GraphElem v, Edge& max_edge)
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
                    {
#if defined(USE_MPI_ONLY)
#else
                        #pragma omp atomic update
#endif
                        ghost_count_[lx] -= 1;
                    }
                    break;
                }
            }
        }

#if defined(OMP_TARGET_OFFLOAD)   
        inline void deactivate_edge_device(GraphElem x, GraphElem y, GraphElem *ghost_count_ptr)
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
                    {
                        #pragma omp atomic update
                        ghost_count_ptr[lx] -= 1;
                    }
                }
            }
        }
#endif
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

                        // mark y<->x inactive, because its matched
                        deactivate_edge(y, x);
                        deactivate_edge(x, y);
                    }
                }
                else // ghost, send REQUEST
                {
                    deactivate_edge(x, y);

                    GraphElem y_x[2] = {y, x};
                    TaggedIsend(MATE_REQUEST_TAG, y_owner, y_x);  
                }
            }
            else // invalidate all neigboring vertices 
            {
                GraphElem e0, e1;
                g_->edge_range(lx, e0, e1);

                for (GraphElem e = e0; e < e1; e++)
                {
                    EdgeActive& edge = g_->get_active_edge(e);
                    const GraphElem z = edge.edge_.tail_;

                    if (edge.active_)
                    {
                        // invalidate x -- z
                        edge.active_ = false;
                        
                        const int z_owner = g_->get_owner(z);
                        
                        if (z_owner == rank_)
                            deactivate_edge(z, x);
                        else // ghost, send INVALID
                        {
                            ghost_count_[lx] -= 1;
                            
                            GraphElem z_x[2] = {z, x};
                            TaggedIsend(MATE_INVALID_TAG, z_owner, z_x);  
                        }
                    }
                }
            }
        }

        // process neighbors of matched vertices
        void process_neighbors(GraphElem v)
        {
            GraphElem e0, e1;
            const GraphElem lv = g_->global_to_local(v);
            g_->edge_range(lv, e0, e1);

            // find unmatched vertices in v's neighborhood
            for (GraphElem e = e0; e < e1; e++)
            {
                EdgeActive& edge = g_->get_active_edge(e);
                
                if (edge.active_)
                {
                    const GraphElem x = edge.edge_.tail_;

                    // v is matched with one and only one of its neighbor,
                    // recalculate new mate for all x in N(v) whose candidate 
                    // mate is v
                    // if mate[v] != x then mate[x] != v, for cases where
                    // mate[x] == v, mate[x] needs to be recalculated
                    
                    // eagerly send INVALID, no need to wait, as 
                    // mate[x] cannot be v, because v is matched 
                    // already
                    if (mate_[g_->global_to_local(v)] != x)
                    {
                        // invalidate v - x, because v is already matched, 
                        // and not with x
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
                        else // send INVALID to invalidate x-v and recompute mate[x] if necessary
                        {      
                            ghost_count_[lv] -= 1; 
                            
                            GraphElem x_v[2] = {x, v};
                            TaggedIsend(MATE_REJECT_TAG, x_owner, x_v);
                        }
                    }
                }   
            }
        }

    private:
        Graph* g_;
        std::vector<GraphElem> D_;
        std::vector<EdgeTuple> M_;
        // count of ghost vertices not owned by me
        std::vector<GraphElem> ghost_count_;
        GraphElem tot_ghosts_;

        // used to store input data
        // and MPI requests for sends
        GraphElem* sbuf_;
        MPI_Request* sreq_;
        GraphElem sbuf_ctr_;

        int rank_, size_;
        MPI_Comm comm_;

        GraphElem* mate_;
};

#endif
