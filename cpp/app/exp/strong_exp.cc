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

class Strong_Exp {
 public:
  Strong_Exp(){}

  void get_query_vfile(int index, std::string &query_vfile) {
    query_vfile = this->query_name + std::to_string(index) + ".v";  
  }

  void get_query_efile(int index, std::string &query_efile) {
    query_efile = this->query_name + std::to_string(index) + ".e";
  }

  void PrintInfo(std::string &query_vfile, std::string &query_efile) {
    LOG(INFO) << "Finish computation.";
    LOG(INFO) << "===============================================";
    LOG(INFO) << "vfile: " + this->vfile;
    LOG(INFO) << "efile: " + this->efile;
    LOG(INFO) << "query_vfile" + query_vfile;
    LOG(INFO) << "query_efile" + query_efile;
    LOG(INFO) << "===============================================";
    LOG(INFO) << "total time: " + std::to_string(get_timer(WORKER_TIMER));
    LOG(INFO) << "load graph timer: " + std::to_string(get_timer(LOAD_TIMER));
    LOG(INFO) << "computing time: " + std::to_string(get_timer(EVALUATION_TIMER));
    LOG(INFO) << "===============================================";
  }

  void Run() {
    StrongSim strongs;
    Graph dgraph;
    GraphLoader dgraph_loader;
    start_timer(WORKER_TIMER);
    start_timer(LOAD_TIMER);
    dgraph_loader.LoadGraph(dgraph, this->vfile, this->efile);
    reset_timer(LOAD_TIMER);
    for (int index = 0; index < this->query_num; index++) {
      Graph qgraph;
      GraphLoader qgraph_loader;
      std::string query_vfile;
      get_query_vfile(index, query_vfile);
      std::string query_efile;
      get_query_efile(index, query_efile);
      qgraph_loader.LoadGraph(qgraph, query_vfile, query_efile);
      start_timer(EVALUATION_TIMER);
      strongs.strong_simulation_sim(dgraph, qgraph);  
      stop_timer(EVALUATION_TIMER);
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
//  std::string add_rm_dir = FLAGS_add_rm_dir;
};

int main(int argc, char **argv) {
  google::SetUsageMessage("Ussage: strong_exp [gflags_opt]");
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::ShutDownCommandLineFlags();
  google::InitGoogleLogging("exp for strong simulation");
  google::ShutdownGoogleLogging();
  Strong_Exp strong_exp;
  strong_exp.Run();
  return 0;
}
