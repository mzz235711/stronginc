#ifndef CPP_STRONGSIMULATION_INC_H_
#define CPP_STRONGSIMULATION_INC_H_
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <unordered_map>
#include <unordered_set>
#include "cpp/core/graphapi.h"
#include "cpp/core/strongr.h"
#include "dual_incremental.h"
#include "cpp/core/ball_view.h"
#include "cpp/utils/util.h"
#include "dualsimulation.h"
class StrongInc{
public:
    StrongInc();

    ~StrongInc();
    template<class T>
    std::unordered_set<T> diff(const std::unordered_set<T> &a, const std::unordered_set<T>& b) {
		std::unordered_set<T> ret;
		for (auto ele : a) {
			if (b.find(ele) == b.end()) {
				ret.insert(ele);
			}
		}
		return ret;
	}

    template<class T>
    std::unordered_set<T> intersection(const std::unordered_set<T> &a, const std::unordered_set<T>& b) {
		std::unordered_set<T> ret;
		for (auto ele : a) {
			if (b.find(ele) != b.end()) {
				ret.insert(ele);
			}
		}
		return ret;
	}

	void find_affected_center_area(Graph &dgraph,std::unordered_set<std::pair<VertexID,VertexID>> &add_edges,
                                                 std::unordered_set<std::pair<VertexID,VertexID>> &rm_edges,
                                                 int d_hop,
                                                 std::unordered_set<VertexID> &result);

    void find_affected_center_area(Graph &dgraph,std::unordered_set<std::pair<VertexID,VertexID>> &add_edges,
                                                 std::unordered_set<std::pair<VertexID,VertexID>> &rm_edges,
                                                 int d_hop,
                                                 std::unordered_set<VertexID> &result,
												 int flag);

      /**
       calculate diameter for qgraph
       */
      int cal_diameter_qgraph(Graph &qgraph);

      /**
       judge node w is a valid node in sim set
       */
      bool valid_sim_w(Graph &qgraph,std::unordered_map<VertexID, std::unordered_set<VertexID>> &sim,VertexID w);

      void find_node_connectivity_nodes(Ball_View &ball_view,std::unordered_set<VertexID> &v_set, VertexID w);

      void rename_sim(Ball_View &ball_view,Graph &qgraph,
                               std::unordered_map<VertexID, std::unordered_set<VertexID>> &sim);

      void ss_counter_initialization(Ball_View &ball_view,Graph &qgraph,
                                     std::unordered_map<VertexID, std::vector<int>> &sim_counter_pre,
                                     std::unordered_map<VertexID, std::vector<int>> &sim_counter_post,
                                     std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w);

      void dual_filter_match(Ball_View &refined_ball_view, Graph &qgraph,
                              std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w,VertexID w,int d_Q);


      void push_phase(Ball_View &ball,Graph &qgraph,VertexID w,int d_Q,
                                  std::set<std::pair<VertexID,VertexID>> &filter_set,
                                  std::unordered_map<VertexID, std::vector<int>> &sim_counter_pre,
                                  std::unordered_map<VertexID, std::vector<int>> &sim_counter_post,
                                  std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w);

      void update_counter(Ball_View &ball,Graph &qgraph,VertexID u,VertexID v,
                          std::unordered_map<VertexID, std::vector<int>> &sim_counter_pre,
                          std::unordered_map<VertexID, std::vector<int>> &sim_counter_post);

      void decremental_refine(Ball_View &ball_view,Graph &qgraph,
                          std::set<std::pair<VertexID,VertexID>> &filter_set,
                          std::unordered_map<VertexID, std::vector<int>> &sim_counter_pre,
                          std::unordered_map<VertexID, std::vector<int>> &sim_counter_post,
                           std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w);

      void extract_max_pg(Ball_View &ball_view,Graph &dgraph,Graph &qgraph,VertexID w,
                                std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w);

      void print_ball_info(Graph &qgraph,std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w,VertexID w);

      void  recalculate_incrementl_dual(Graph &dgraph, Graph &qgraph,
                                      std::unordered_map<VertexID,std::unordered_set<VertexID>> &dsim,
                                      std::unordered_set<std::pair<VertexID,VertexID>> &add_edges,
                                      std::unordered_set<std::pair<VertexID,VertexID>> &rm_edges);

									  
	void cal_culculate_inc_dhop_nodes_add(Graph &dgraph, int d_Q, std::unordered_set<VertexID> &ball_node,
									  std::vector<int> &dis, std::unordered_set<std::pair<VertexID,VertexID>> &add_edges);
	
	void strong_simulation_inc_only_add(Graph &dgraph, Graph &qgraph,
                               std::unordered_map<VertexID,std::unordered_set<VertexID>> &dsim,
                               std::vector<StrongR> &strong_r,
                               std::unordered_set<std::pair<VertexID,VertexID>> &add_edges,
							   int flag,
							   std::unordered_map<VertexID,std::unordered_set<VertexID>> &whole_ball_nodes,
							   std::unordered_map<VertexID,std::vector<int>> &whole_dist);
	
    void strong_simulation_inc(Graph &dgraph, Graph &qgraph,
                               std::unordered_map<VertexID,std::unordered_set<VertexID>> &dsim,
                               std::vector<StrongR> &strong_r,
                               std::unordered_set<std::pair<VertexID,VertexID>> &add_edges,
                               std::unordered_set<std::pair<VertexID,VertexID>> &rm_edges);
									  
	void strong_simulation_inc(Graph &dgraph, Graph &qgraph,
                               std::unordered_map<VertexID,std::unordered_set<VertexID>> &dsim,
                               std::vector<StrongR> &strong_r,
                               std::unordered_set<std::pair<VertexID,VertexID>> &add_edges,
                               std::unordered_set<std::pair<VertexID,VertexID>> &rm_edges,
							   int flag);
};
#endif //CPP_STRONGSIMULATION_INC_H_
