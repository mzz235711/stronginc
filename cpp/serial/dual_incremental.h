#ifndef CPP_DUALSIMULATION_INC_H_
#define CPP_DUALSIMULATION_INC_H_
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <unordered_map>
#include <unordered_set>
#include "cpp/core/graphapi.h"
#include "cpp/core/ball_view.h"
 /**
  * algorithm for dualsimulation incremental
  */
class DualInc {
 public:
  DualInc();
  ~DualInc();
 /**
  * a node set difference another node set
  */
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

  void propagate_add(Graph &dgraph,Graph &qgraph,
                       std::set<std::pair<VertexID,VertexID>> &candidate_node,
                       std::vector<std::unordered_set<VertexID>> &aff_node,
                       std::vector<std::unordered_set<VertexID>> &dsim,
                       std::set<std::pair<VertexID,VertexID>> &already_matched);

  void update_pre_dec_counter(GraphView &graph_view,Graph &qgraph,VertexID u,VertexID v,
                          std::vector<std::vector<int>> &sim_counter_pre,
                          std::vector<std::vector<int>> &sim_counter_post);

  void propagate_remove(GraphView &graph_view,Graph &qgraph,
                          std::vector<std::unordered_set<VertexID>> &aff_node,
                          std::set<std::pair<VertexID,VertexID>> &filter_set,
                          std::vector<std::vector<int>> &sim_counter_pre,
                          std::vector<std::vector<int>> &sim_counter_post,
                          std::set<std::pair<VertexID,VertexID>> &already_matched);

  void incremental_addedges(Graph &dgraph,Graph &qgraph,
                           std::vector<std::unordered_set<VertexID>> &dsim,
                           std::set<std::pair<VertexID,VertexID>> &add_edges);

  void update_counter(Graph &dgraph,Graph &qgraph,VertexID u,VertexID v,
                          std::vector<std::vector<int>> &sim_counter_pre,
                          std::vector<std::vector<int>> &sim_counter_post);

  void decremental_rmove(Graph &dgraph,Graph &qgraph,
                          std::set<std::pair<VertexID,VertexID>> &filter_set,
                          std::vector<std::unordered_set<VertexID>> &dsim,
                          std::vector<std::vector<int>> &sim_counter_pre,
                          std::vector<std::vector<int>> &sim_counter_post);

  void incremental_removeedgs(Graph &dgraph,Graph &qgraph,
                           std::vector<std::unordered_set<VertexID>> &dsim,
                           std::set<std::pair<VertexID,VertexID>> &rm_edges);

  void dual_incremental(Graph &dgraph,Graph &qgraph,
                          std::vector<std::unordered_set<VertexID>> &dsim,
                          std::set<std::pair<VertexID,VertexID>> &add_edges,std::set<std::pair<VertexID,VertexID>> &rm_edges);
 private:
   int n;
  //Graph qgraph;
  //Graph dgraph;
  //GraphLoader testgraph_loader_;
};

#endif //CPP_DUALSIMULATION_INC_H_