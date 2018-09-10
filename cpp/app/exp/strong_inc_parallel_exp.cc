#include <gflags/gflags.h>
#include <glog/logging.h>
#include<boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
//#include "cpp/core/graphapi.h"
#include "cpp/serial/dualsimulation.h"
#include "cpp/serial/dual_incremental.h"
#include "cpp/serial/strongsimulation.h"
#include "cpp/serial/strong_incremental.h"
#include "cpp/utils/util.h"
#include "cpp/io/io_local.h"
#include "cpp/core/global.h"
#include "cpp/core/strongr.h"
#include "cpp/core/view.h"
#include "cpp/utils/generate.h"
#include "cpp/utils/time.h"
#include<iostream>
#include <fstream>
#include<ctime>

#include<boost/filesystem.hpp>

class Strong_Inc_Parallel_Exp {
 public:
  Strong_Inc_Parallel_Exp() {}

  void get_query_vfile(int index, std::string &query_vfile) {
    query_vfile = this->query_name + std::to_string(index) + ".v";   
  }

  void get_query_efile(int index, std::string &query_efile) {
    query_efile = this->query_name + std::to_string(index) + ".e";
  }

  void PrintInfo(std::string &query_vfile, std::string &query_efile, int extend) {
    std::string add_efile = this->add_name + std::to_string(extend) + ".e";
    std::string add_vfile = this->add_name + std::to_string(extend) + ".v";
    std::string rm_efile = this->rm_name + std::to_string(extend) + ".e";
    LOG(INFO) << "Finish computation.";
    LOG(INFO) << "===============================================";
    LOG(INFO) << "vfile: " + this->vfile;
    LOG(INFO) << "efile: " + this->efile;
    LOG(INFO) << "query_vfile: " + query_vfile;
    LOG(INFO) << "query_efile: " + query_efile;
    LOG(INFO) << "add_efile: " + add_efile;
    LOG(INFO) << "add_vfile: " + add_vfile;
    LOG(INFO) << "remove_efile: " + rm_efile;
    LOG(INFO) << "total time: " + std::to_string(get_timer(WORKER_TIMER));
    LOG(INFO) << "load graph timer: " + std::to_string(get_timer(LOAD_TIMER));
    LOG(INFO) << "strong simulation time: " + std::to_string(get_timer(EVALUATION_TIMER));
    LOG(INFO) << "incremental strong simulation time: "
                    + std::to_string(get_timer(INCREMENTAL_TIMER));
    LOG(INFO) << "===============================================";
  }

  bool Verify(std::vector<StrongR>& result1, std::vector<StrongR>& result2) {
    std::unordered_set<VertexID> global_center, inc_center;
    for (auto &ball : result1) {
      global_center.insert(ball.center());
    }
    for (auto &ball : result2) {
      inc_center.insert(ball.center());
    }
    if (global_center.size() != inc_center.size()) {
      return false;
    }
     return true;
  }

  void Run() {
    Graph dgraph, fragmentgraph;
    GraphLoader dgraph_loader;
    dgraph_loader.LoadGraph(dgraph, this->vfile, this->efile);
    Fragment fragment(fragmentgraph, this->vfile, this->efile);
    for (int index = 1; index <= this->query_num; index++) {
      std::string query_vfile;
      std::string query_efile;
      get_query_vfile(index, query_vfile);
      get_query_efile(index, query_efile);
      Graph qgraph;
      qgraph_loader.LoadGraph(qgraph, query_vfile, query_efile);
      Dual_Parallel dual_parallel;
      std::unordered_map<VertexID, std::unordered_set<VertexID>> partial_sim;
      dual_parallel.dual_parallel(fragment, fragmentgraph, qgraph, partial_sim);
      StrongparallelInc strong_parallel_dir;
      std::vector<StrongR> partial_result =
           strong_parallel_dir.strong_parallel(fragment, fragmentgraph, qgraph);
      for (int extend = 1; extend <= this->extend_num; extend++) {
        Graph inc_dgraph;
        Graph inc_fragmentgraph;
//        dgraph_loader.LoadGraph(inc_dgraph, this->vfile, this->efile);
        std::unordered_set<std::pair<VertexID, VertexID>> add_edges, rm_edges;
        std::vector<std::pair<VertexID, Vertex>> add_vertices;
        Load_bunch_edges(add_vertices, add_edges, add_name, extend);
        Load_bunch_edges(rm_edges, rm_name, extend);
        for (auto &v : add_vertices) {
          dgraph.AddVertex(Vertex(v.first, v.second));
        }
        for (auto &e : add_edges) {
          dgraph.AddEdge(Edge(e.first, e.second, 1));
        }
        for (auto &e : rm_edges) {
          dgraph.RemoveEdge(Edge(e.first, e.second, 1));
        }
        std::unordered_map<VertexID, std::unordereed_set<VertexID>> inc_parallel_dual;
        std::vector<StrongR> inc_parallel_strong;
        for (auto u : qgraph.GetAllVerticesID()) {
          for (auto &v : partial_sim[u]) {
            inc_parallel_dual[u].insert(v);
          }
        }
        for (auto &ball : partial_strong) {
          inc_parallel_strong.push_back(ball);
        }
        StrongparallelInc strong_parallel_inc;
        std::vector<StrongR> partial_result =
          strong_parallel_inc.strong_parallel_inc(inc_fragment, inc_fragmentgraph,
          qgraph, inc_parallel_dual, inc_parallel_strong, add_edges, rm_edges);
        if (fid == 0) {
          std::vector<std::vector<StrongR>> tmp_vec(get_num_workers());
          tmp_vec[fid] = partial_result;
          masterGather(tmp_vec);
          std::vector<StrongR> global_result;
          for (int i = 0; i < get_num_workers(); i++) {
            for (auto &ball : tmp_vec[i]) {
              global_result.push_back(ball);
            }
          }
          std::unordered_set<VertexID> global_center, inc_center;
          StrongSim strongsim;
          std::vector<StrongR> direct_strong = strongsim.strong_simulation_sim(dgraph, qgraph);
          if (!Verify(global_result, direct_strong)) {
            LOG(INFO) << "WARNING! WE ARE UNDER ATTACK!";
          }
          for (auto &e : add_edges) {
            dgraph.RemoveEdge(Edge(e.first, e.second, 1));
          }
        }
      }  
    }
  }
 private:
  std::string vfile = FLAGS_vfile;
  std::string efile = FLAGS_efile;
  std::string rfile = FLAGS_rfile;
  std::string query_name = FLAGS_query_file;
  int query_num = FLAGS_query_num;
  std::string add_name = FLAGS_base_add_file;
  std::string rm_name = FLAGS_base_remove_file;
  int extend_num = FLAGS_extend_num;
}

int main(int argc, char **argv) {
  google::SetUsageMessage("Usage: strong_inc_parallel_exp [gflags_opt]");
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::ShutDownCommandLindFlags();
  google::InitGoogleLogging("");
  google::ShutdownGoogleLogging();
  init_workers();
  int rank = get_worker_id();
  Strong_Inc_Parallel_Exp strong_inc_parallel_exp;
  strong_inc_parallel_exp.Run();
  worker_finalize();
}