#include <gflags/gflags.h>
#include <glog/logging.h>
#include<boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
//#include "cpp/core/graphapi.h"
#include "cpp/serial/dualsimulation.h"
#include "cpp/serial/dual_incremental.h"
#include "cpp/serial/strongsimulation.h"
#include "cpp/serial/strong_incremental.h"
#include "cpp/parallel/dualsimulation_parallel.h"
#include "cpp/parallel/dual_parallel_inc.h"
#include "cpp/parallel/strongparallel_incremental.h"
#include "cpp/utils/util.h"
#include "cpp/io/io_local.h"
#include "cpp/core/Fragment.h"
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

  void PrintInfo(std::string &query_vfile, std::string &query_efile,
                 int extend, std::vector<double> inc_timers,
                 std::vector<double> dir_timers) {
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
    LOG(INFO) << "parallel strong simulation time: ";
    for (int i = 0; i < dir_timers.size(); i++) {
      LOG(INFO) << "    --frag" << i << ": " << std::to_string(dir_timers[i]);
    }
    LOG(INFO) << "parallel incremental strong simulation time: ";
    for (int i = 0; i < inc_timers.size(); i++) {
      LOG(INFO) << "    --frag" << i << ": " << std::to_string(inc_timers[i]);
    }
                    
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
    init_timers();
    int fid = get_worker_id();
//    Graph dgraph, fragmentgraph;
      Graph fragmentgraph;
//    GraphLoader dgraph_loader;
//    dgraph_loader.LoadGraph(dgraph, this->vfile, this->efile);
    Fragment fragment(fragmentgraph, this->vfile, this->efile, this->rfile);
    for (int index = 1; index <= this->query_num; index++) {
      std::string query_vfile;
      std::string query_efile;
      get_query_vfile(index, query_vfile);
      get_query_efile(index, query_efile);
      Graph qgraph;
      GraphLoader qgraph_loader;
      qgraph_loader.LoadGraph(qgraph, query_vfile, query_efile);
      Dual_Parallel dual_parallel;
      std::unordered_map<VertexID, std::unordered_set<VertexID>> partial_sim;
      dual_parallel.dual_paraller(fragment, fragmentgraph, qgraph, partial_sim);
      StrongparallelInc strong_parallel_dir;
      std::vector<StrongR> partial_strong =
           strong_parallel_dir.strong_parallel(fragment, fragmentgraph, qgraph);
      LOG(INFO) << "xxxxxxxxxxxxx";
      for (int extend = 1; extend <= this->extend_num; extend++) {
//        Graph inc_dgraph;
//        Graph inc_fragmentgraph;
//        dgraph_loader.LoadGraph(inc_dgraph, this->vfile, this->efile);
        std::unordered_set<std::pair<VertexID, VertexID>> add_edges, rm_edges;
        std::vector<std::pair<VertexID, VertexLabel>> add_vertices;
        Load_bunch_edges(add_vertices, add_edges, add_name, extend);
        Load_bunch_edges(rm_edges, rm_name, extend);
        fragment.update_by_add_vertices(fragmentgraph, add_vertices);
        fragment.update_by_add_edges(fragmentgraph, add_edges, true);
        fragment.update_by_remove_edges(fragmentgraph, rm_edges, true);
        std::unordered_map<VertexID, std::unordered_set<VertexID>> inc_parallel_dual;
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
        worker_barrier();
        start_timer(INCREMENTAL_TIMER);
        std::vector<StrongR> inc_result =
          strong_parallel_inc.strong_parallel_inc(fragment, fragmentgraph,
          qgraph, inc_parallel_dual, inc_parallel_strong, add_edges, rm_edges);
        stop_timer(INCREMENTAL_TIMER);
        worker_barrier();

        start_timer(EVALUATION_TIMER);
        std::vector<StrongR> dir_result =
           strong_parallel_dir.strong_parallel(fragment, fragmentgraph, qgraph);
        stop_timer(EVALUATION_TIMER);
        worker_barrier();
        if (fid == 0) {
          std::vector<std::vector<StrongR>> tmp_vec(get_num_workers());
          std::vector<double> inc_timers(get_num_workers());
          tmp_vec[fid] = inc_result;
          masterGather(tmp_vec);
          inc_timers[fid] = get_timer(INCREMENTAL_TIMER);
          masterGather(inc_timers);
          std::vector<StrongR> global_inc_result;
          for (int i = 0; i < get_num_workers(); i++) {
            for (auto &ball : tmp_vec[i]) {
              global_inc_result.push_back(ball);
            }
          }

          std::vector<std::vector<StrongR>> dir_vec(get_num_workers());
          std::vector<double> dir_timers(get_num_workers());
          dir_vec[fid] = dir_result;
          masterGather(dir_vec);
          dir_timers[fid] = get_timer(EVALUATION_TIMER);
          masterGather(dir_timers);
          std::vector<StrongR> global_dir_result;
          for (int i = 0; i < get_num_workers(); i++) {
            for (auto &ball : dir_vec[i]) {
              global_dir_result.push_back(ball);
            }
          }
          if (!Verify(global_inc_result, global_dir_result)) {
            LOG(INFO) << "WARNING! WE ARE UNDER ATTACK!";
          }
          PrintInfo(query_vfile, query_efile, extend, dir_timers, inc_timers);
        } else {
          slaveGather(inc_result);
          double tmp_timer = get_timer(INCREMENTAL_TIMER);
          slaveGather(tmp_timer);
          slaveGather(dir_result);
          tmp_timer = get_timer(EVALUATION_TIMER);
          slaveGather(tmp_timer);
        }
        fragment.update_by_add_edges(fragmentgraph, rm_edges, true);
        fragment.update_by_remove_edges(fragmentgraph, add_edges, true);
        reset_timer(EVALUATION_TIMER);
        reset_timer(INCREMENTAL_TIMER);  
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
};

int main(int argc, char **argv) {
  google::SetUsageMessage("Usage: strong_inc_parallel_exp [gflags_opt]");
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::ShutDownCommandLineFlags();
  google::InitGoogleLogging("");
  google::ShutdownGoogleLogging();
  init_workers();
  int rank = get_worker_id();
  Strong_Inc_Parallel_Exp strong_inc_parallel_exp;
  strong_inc_parallel_exp.Run();
  worker_finalize();
}
