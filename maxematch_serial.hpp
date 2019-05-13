
/* Serial implementation of maximal edge matching.
 * This version does not uses the EdgeMatch object. 
 */

#pragma once
#ifndef MAXEMATCHSER_HPP
#define MAXEMATCHSER_HPP

#include "graph.hpp"
#include <numeric>

class MaxEdgeMatch
{
    public:
        MaxEdgeMatch(Graph* g): 
            g_(g), D_(0), M_(0), 
            mate_(0)
        {}
        ~MaxEdgeMatch() {}

        void clear()
        {
            D_.clear();
            M_.clear();
            mate_.clear();
        }

        void print_M() const
        {
            std::cout << "Matched vertices: " << std::endl;
            for (GraphElem i = 0; i < M_.size(); i++)
                std::cout << M_[i].ij_[0] << " ---- " << M_[i].ij_[1] << std::endl;
        }
         
        // if mate[mate[v]] == v then
        // we're good
        void check_results()
        {
            bool success = true;
            for (GraphElem i = 0; i < M_.size(); i++)
            {                
                if ((mate_[mate_[M_[i].ij_[0]]] != M_[i].ij_[0])
                        || (mate_[mate_[M_[i].ij_[1]]] != M_[i].ij_[1]))
                {
                    std::cout << "Validation FAILED." << std::endl; 
                    std::cout << "mate_[mate_[" << M_[i].ij_[0] << "]] != " << M_[i].ij_[0] << std::endl;
                    std::cout << "mate_[mate_[" << M_[i].ij_[1] << "]] != " << M_[i].ij_[1] << std::endl;
                    success = false;

                }
            }
            if (success)
                std::cout << "Validation SUCCESS." << std::endl;
        }
        
        // TODO FIXME not expecting a, b to
        // be large, if large then following
        // absolute tolerance test will fail:
        // http://realtimecollisiondetection.net/blog/?p=89
        bool is_same(GraphWeight a, GraphWeight b) 
        { return std::abs(a - b) <= std::numeric_limits<GraphWeight>::epsilon(); }
 
        // expecting v to be local index
        // may require global_to_local
        // before passing v
        void heaviest_edge_unmatched(GraphElem v, Edge& max_edge, GraphElem x = -1)
        {
            GraphElem e0, e1;
            g_->edge_range(v, e0, e1);

            for (GraphElem e = e0; e < e1; e++)
            {
                EdgeActive& edge = g_->get_active_edge(e);
                if (edge.active_)
                {
                    if (edge.edge_.tail_ == x)
                        continue;

                    if ((mate_[edge.edge_.tail_] == -1) 
                            || (mate_[mate_[edge.edge_.tail_]] 
                                != edge.edge_.tail_))
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
        }

        void inactivate_edge(GraphElem x, GraphElem y)
        {
            GraphElem e0, e1;
            g_->edge_range(x, e0, e1);

            for (GraphElem e = e0; e < e1; e++)
            {
                EdgeActive& edge = g_->get_active_edge(e);
                if (edge.edge_.tail_ == y)
                {
                    edge.active_ = false;
                    break;
                }
            }
        }

        // maximal edge matching
        std::vector<EdgeTuple> const& operator()()
        {
            maxematch();
            return M_;
        }

        // maximal edge matching
        void maxematch()
        {
            // initializations
            GraphElem lnv = g_->get_lnv();
            mate_.resize(lnv);
            std::fill(mate_.begin(), mate_.end(), -1);

            /* Phase #1 */
            // part 1: compute max edge for every vertex
            for (GraphElem v = 0; v < lnv; v++)
            {
                Edge max_edge;
                heaviest_edge_unmatched(v, max_edge);
                GraphElem u = mate_[v] = max_edge.tail_; // v's mate
                
                // is mate[u] == v?
                if (mate_[u] == v) // matched
                {
                    D_.push_back(u);
                    D_.push_back(v);
                    M_.emplace_back(u, v, max_edge.weight_);
                        
                    inactivate_edge(v, u);
                    inactivate_edge(u, v);
                }
            }

            /* Phase #2 */
            unsigned int remote_count = 0;
            while(1)
            {     
                // exit criteria
                if (D_.size() == 0) 
                    break;
#if defined(DEBUG) && DEBUG > 0
                std::cout << "Start of iteration #" << remote_count << 
                    " : Size of D = " << D_.size() << std::endl;
#endif
                if (D_.size() > 0)
                {
                    GraphElem v = D_.back();
                    D_.pop_back();
                    update_mate(v);
                }

                remote_count++;
            } // end of while(D_)
        }

        // check if mate[x] = v and mate[v] != x
        // if yes, compute mate[x]
        void update_mate(GraphElem v)
        {
            GraphElem e0, e1;
            g_->edge_range(v, e0, e1);

            for (GraphElem e = e0; e < e1; e++)
            {
                Edge const& edge = g_->get_edge(e);
                GraphElem const& x = edge.tail_;

                // check if vertex is already matched
                auto result = std::find_if(M_.begin(), M_.end(), 
                        [&](EdgeTuple const& et) 
                        { return (((et.ij_[0] == v) || (et.ij_[1] == v)) && 
                                ((et.ij_[0] == x) || (et.ij_[1] == x))); });
                
                //  mate[x] == v and (v,x) not in M
                if ((mate_[x] == v) && (result == std::end(M_)))
                {
                    Edge x_max_edge;
                    heaviest_edge_unmatched(x, x_max_edge, v);
                    GraphElem y = mate_[x] = x_max_edge.tail_;

                    if (y == -1) // if x has no neighbor other than v
                        continue;

                    if (mate_[y] == x) // matched
                    {
                        D_.push_back(x);
                        D_.push_back(y);
                        M_.emplace_back(x, y, x_max_edge.weight_);

                        inactivate_edge(x, y);
                    }
                }
            }
        }

    private:
        Graph* g_;
        
        std::vector<GraphElem> D_;
        std::vector<EdgeTuple> M_;
        std::vector<GraphElem> mate_;
};

#endif
