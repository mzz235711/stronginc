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
#include<iostream>
#include <fstream>
#include<ctime>

#include<boost/filesystem.hpp>

class Strong_Inc_Exp {
 public:
  Strong_Inc_Exp(){}

  void get_query_vfile(int index, std::string &query_vfile) {
    query_vfile = this->query_name + std::to_string(index) + ".v";  
  }

  void get_query_efile(int index, std::string &query_efile) {
    query_efile = this->query_name + std::to_string(index) + ".e";
  }

  void PrintInfo(std::string &query_vfile, std::string &query_efile, int extend) {
    std::string add_efile = this->add_name + std::to_string(extend) + ".e"
    std::string rm_efile = this->rm_name + std::to_string(extend) + ".e";
    LOG(INFO) << "Finish computation.";
    LOG(INFO) << "==============================================="
    LOG(INFO) << "vfile: " + this->graph_vfile;
    LOG(INFO) << "efile: " + this->graph_efile;
    LOG(INFO) << "query_vfile" + query_vfile;
    LOG(INFO) << "query_efile" + query_efile;
    LOG(INFO) << "total time: " + get_timer(WORKER_TIMER);
    LOG(INFO) << "load graph timer: " + get_timer(LOAD_TIMER);
    LOG(INFO) << "dual simulation time: " + get_timer(EVALUATION_TIMER);
    LOG(INFO) << "incremental dual simulation time: " + get_timer(INCREMENTAL_TIMER);
    LOG(INFO) << "==============================================="
  }

  void Run() {
    init_timer();
    StrongInc strong_inc;
    StrongSim strong_sim;
    DualSim dual_sim;
    Graph dgraph;
    GraphLoader dgraph_loader;
    start_timer(WORKER_TIMER);
    start_timer(LOAD_TIMER);  
    dgraph_loader.LoadGraph(dgraph, this->vfile, this->efile);
    stop_timer(LOAD_TIMER);
    for (int index = 0; index < this->query_num; index++) {
      std::string query_vfile;
      std::string query_efile;
      get_query_vfile(index, query_vfile);
      get_query_efile(index, query_efile);
      Graph qgraph;
      GraphLoader qgraph_loader;
      qgraph_loader.LoadGraph(qgraph, query_vfile, query_efile);
      std::unordered_map<VertexID, std::unordered_set<VertexID>> sim;
      std::vector<StrongR> strongsimr = strongsim.strong_simulation_sim(dgraph, qgraph);
      dual_sim.dual_simulation(dgraph, qgraph, sim, false);
      for (int extend = 0; extend < this->extend_num; extend++) {
        std::set<std::pair<VertexID, VertexID>> add_edges, rm_edges;
        Load_bunch_edges(add_edgesm, base_add_file, extend);
        Load_bunch_edges(rm_edges, extend, j);
        std::vector<StrongR> tmp_r ;
        std::unordered_map<VertexID, std::unordered_set<VertexID>> tmp_sim;
        for (auto &ball : strongsimr) {
          tmp_r.push_bakc(ball);
        }
        for (auto &u : qgraph.GetAllVerticesID()) {
          tmp_sim[u] = std::unordered_set<VertexID>;
          for (auto &v : sim[u]) {
            tmp_sim[u].insert(v);
          }
        }
        StrongSim strong_sim_dir;
        for (auto &e : add_edges) {
          dgraph.AddEdge(Edge(e.first, e.second, 1));
        }
        for (auto &e : rm_edges) {
          dgraph.RemoveEdge(Edge(e.first, e.second, 1));
        }
        start_timer(EVALUATION_TIMER);
        std::vector<StrongR> simr_dir = strong_sim_dir.strong_simulation_sim(dgraph, qgraph);
        stop_timer(EVALUATION_TIMER);
        start_timer(INCREMENTAL_TIMER);
        stronginc.strong_simulation_inc(dgraph, qgraph, tmp_sim, tmp_r, add_edges, rm_edges);
        stop_timer(INCREMENTAL_TIMER);
        PrintInfo(query_vfile, query_efile, extend);
        reset_timer(EVALUATION_TIMER);
        reset_timer(INCREMENTAL_TIMER);
        for (auto &e : rm_edges) {
          dgraph.AddEdge(Edge(e.first, e.second, 1));
        }
        for (auto &e : add_edges) {
          dgraph.RemoveEdge(Edge(e.first, e.second, 1));
        }
        
      }  
    }

  }

 private:
  std::string vfile = FLAGS_vfile;
  std::string efile = FLAGS_efile;
  std::string viewfile = FLAGS_viewfile;
  std::string query_name = FLAGS_query_file;
  int query_num = FLAGS_query_num;
  std::string add_rm_dir = FLAGS_add_rm_dir;
};

int main(int argc, char **argv) {
  google::SetUsageMessage("Usage: test [gflags_opt]");
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::ShutDownCommandLineFlags();
  google::InitGoogleLogging("");
  google::ShutDownGoogleLogging();
  Strong_Inc_Exp strong_inc_exp();
  strong_inc_exp.Run();
  return 0;
}