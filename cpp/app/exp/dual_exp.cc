#include <gflags/gflags.h>
#include <glog/logging.h>
#include<boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include "cpp/serial/dualsimulation.h"
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

class Dual_Exp {
 public:
  Dual_Exp() {}

  void get_query_vfile(int index, std::string &query_vfile) {
    query_vfile = this->query_name + std::to_string(index) + ".v";
  }

  void get_query_efile(int index, std::string &query_efile) {
    query_efile = this->query_name + std::to_string(index) + ".e";
  }

  void PrintInfo(std::string &query_vfile, std::string &query_efile) {
    LOG(INFO) << "Finish computation.";
    LOG(INFO) << "==============================================="
    LOG(INFO) << "vfile: " + this->graph_vfile;
    LOG(INFO) << "efile: " + this->graph_efile;
    LOG(INFO) << "query_vfile" + query_vfile;
    LOG(INFO) << "query_efile" + query_efile;
    LOG(INFO) << "==============================================="
    LOG(INFO) << "total time: " + get_timer(WORKER_TIMER);
    LOG(INFO) << "load graph timer: " + get_timer(LOAD_TIMER);
    LOG(INFO) << "computing time: " + get_timer(EVALUATION_TIMER);
    LOG(INFO) << "==============================================="
  }
  
  void run() {
    init_timers();  
    Graph dgraph;
    GraphLoader dgraph_loader;
    start_timer(WORKER_TIMER);
    start_timer(LOAD_TIMER);
    dgraph_loader.LoadGraph(dgraph, vfile, efile);
    stop_timer(LOAD_TIMER);
    DualSim dualsim;
    for(int index = 0; index < query_num; index++) {
      Graph qgraph;
      GraphLoader dgraph_loader;
      std::string query_vfile, query_efile;
      get_query_vfile(index, query_vfile);
      get_query_efile(index, query_efile);
      qgraph_loader.LoadGraph(qgraph, query_vfile, query_efile);
      std::unorerder_map<VertexID, std::unordered_set<VertexID>> sim;
      start_timer(EVALUATION_TIMER);
      dualsim.dual_simulation(dgraph, qgraph, sim, false);
      stop_timer(EVALUATION_TIMER);
      stop_timer(WORKER_TIMER);
      PrintInfo(query_vfile, query_efile);
      reset_timer(EVALUATION_TIMER);
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

int main (int argc, char **argv) {
  google::SetUsageMessage("Usage: dual_exp [gflags_opt]");
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::ShutDownCommandLineFlags();
  google::InitGoogleLogging("exp for dual simualtion");
  google::ShutdownGoogleLogging();
  Dual_Exp dual_exp();
  dual_exp.run();
  return 0;
}