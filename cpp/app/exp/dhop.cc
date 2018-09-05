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

class DHop {
 public:
  DHop() {}

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
    LOG(INFO) << "query_vfilei: " + query_vfile;
    LOG(INFO) << "query_efile: " + query_efile;
    LOG(INFO) << "add_efile: " + add_efile;
    LOG(INFO) << "add_vfile: " + add_vfile;
    LOG(INFO) << "remove_efile: " + rm_efile;
    LOG(INFO) << "total time: " + std::to_string(get_timer(WORKER_TIMER));
    LOG(INFO) << "load graph timer: " + std::to_string(get_timer(LOAD_TIMER));
    LOG(INFO) << "dual simulation time: " + std::to_string(get_timer(EVALUATION_TIMER));
    LOG(INFO) << "incremental dual simulation time: "
                    + std::to_string(get_timer(INCREMENTAL_TIMER));
    LOG(INFO) << "===============================================";
  }

void dhop_add_origin(Graph& dgraph,int d_Q,
                    std::unordered_set<VertexID> &result, std::vector<int> &dis,
                    std::unordered_set<std::pair<VertexID,VertexID>> &add_edges){
    int dgraph_num_vertices = dgraph.GetNumVertices();
    dis.resize(dgraph_num_vertices,INT_MAX);
    VertexID base_id = 0;
    VertexID inc_id = 0;
    for(auto e : add_edges){
        if(dis[e.first]>=d_Q && dis[e.second]>=d_Q){
            continue;
        }else if(dis[e.first]-2 >= dis[e.second]){
            base_id = e.second; 
            inc_id = e.first;
        }else if(dis[e.second]-2 >= dis[e.first]){
            base_id = e.first;
            inc_id = e.second;
        }else{
            continue;
        }   
        std::queue<VertexID> q;
        dis[inc_id] = dis[base_id] + 1;
                // if(dis[inc_id] <= d_Q){
            result.insert(inc_id);
            q.push(inc_id);
        // }
        while(!q.empty()){
            VertexID root = q.front();
            q.pop(); 
            if(dis[root] == d_Q){
                break ;
            }   
            for(auto v : dgraph.GetChildrenID(root)){
                if(dis[v] > dis[root]+1){
                    q.push(v);
                    result.insert(v);
                    dis[v] = dis[root] + 1;
                }
            }
            for (auto v : dgraph.GetParentsID(root)){
                                if(dis[v] > dis[root]+1){
                                        q.push(v);
                                        result.insert(v);
                                        dis[v] = dis[root] + 1;
                                }       
            }                   
        }   
    }   
}

  void dhop_direct(Graph &dgraph, VertexID vid, int d_hop,
                   std::unordered_set<VertexID> &result, std::vector<int> &dis) {
    int tvnum = dgraph.GetNumVertices();
    dis.resize(tvnum, max);
    std::queue<VertexID> source;
    source.push(vid);
    dis[vid] = 0;
    std::vector<VertexID> visit(tvnum, false);
    while (!source.empty()) {
      VertexID u = source.front();
      source.pop();
      if (visit[u]) continue;
      visit[u] = true;
      for (auto v : dgraph.GetChildrenID(u)) {
        if (dis[v] == max) {
          dis[v] = dis[u] + 1;
          if (dis[v] < d_hop) {
            source.push(v);
          }
         result.insert(v);
        }
      }
      for (auto v : dgraph.GetParentsID(u)) {
        if (dis[v] == max) {
          dis[v] = dis[u] + 1;
          if (dis[v] < d_hop) {
            source.push(v);
          }
          result.insert(v);
        }
      }
    }
  }

  
  void dhop_add(Graph& dgraph, int d_Q, std::unordered_set<VertexID> &result,
                std::vector<int> &dis,
                std::unordered_set<std::pair<VertexID,VertexID>> &add_edges) {  
    std::priority_queue<std::pair<int, VertexID>> source;
    int dgraph_num_vertices = dgraph.GetNumVertices();
    dis.resize(dgraph_num_vertices,INT_MAX);
    for (auto &e : add_edges) {
      if ((dis[e.first] > dis[e.second] + 1) && dis[e.second] < d_Q) {
        dis[e.first] = dis[e.second] + 1;
      result.insert(e.first);
        if (dis[e.first] < d_Q) {
          source.push(make_pair(-dis[e.first], e.first));
        }
      } else if ((dis[e.second] > dis[e.first] + 1) && dis[e.first] < d_Q) {
        dis[e.second] = dis[e.first] + 1;
        result.insert(e.second);
        if (dis[e.second] < d_Q) {
          source.push(make_pair(-dis[e.second], e.second));
        }
      }
    } 
    int tvnum = dgraph.GetNumVertices();
    std::vector<bool> visit(tvnum, false);
//    if (source.size() != 0) {
//      LOG(INFO) << "source: " << source.size();
//    }
    while (!source.empty()) {
      VertexID u = source.top().second;
      int dist = -source.top().first;
      source.pop();
      if(visit[u])  continue;
      visit[u] = true;
      for (auto v : dgraph.GetChildrenID(u)) {
        if (dis[v] > dis[u] + 1) {
          dis[v] = dis[u] + 1;
          result.insert(v);
          if (dis[v] < d_Q) {
            source.push(make_pair(-dis[v], v));
          }
        }
      }
      for (auto v : dgraph.GetParentsID(u)) {
        if (dis[v] > dis[u] + 1) {
          dis[v] = dis[u] + 1;
          result.insert(v);
          if (dis[v] < d_Q) {
            source.push(make_pair(-dis[v], v));
          }
        }
      }
    }
  }

  bool Verify(std::vector<int> &dis1, std::vector<int> &dis2) {
    for (int i = 0; i < dis1.size();i++) {
      if (dis1[i] != dis2[i]) {
        LOG(INFO) << i << "\t" << dis1[i] << "\t" << dis2[i];
        return false;
      }
    }
    return true;
  }
  void Run() {
    init_timers();
    DualSim dualsim;
    Graph dgraph;
    GraphLoader dgraph_loader;
    start_timer(WORKER_TIMER);
    start_timer(LOAD_TIMER);
    dgraph_loader.LoadGraph(dgraph, this->vfile, this->efile);
    stop_timer(LOAD_TIMER);
    for (int index = 1; index <= query_num; index++) {
      Graph qgraph;
      GraphLoader qgraph_loader;
      std::string query_vfile;
      get_query_vfile(index, query_vfile);
      std::string query_efile;
      get_query_efile(index,  query_efile);
      std::unordered_map<VertexID, std::unordered_set<VertexID>> origin_sim;
      bool initialization = false;
      qgraph_loader.LoadGraph(qgraph, query_vfile, query_efile);
      int d_Q = cal_diameter_qgraph(qgraph);
      LOG(INFO) << "d_Q: " << d_Q;
      dualsim.dual_simulation(dgraph, qgraph, origin_sim, initialization);
      std::vector<int> dis;
      std::unordered_set<VertexID> match;
      for (auto u : qgraph.GetAllVerticesID()) {
        for (auto &v : origin_sim[u]) {
          match.insert(v);
        }
      }
      LOG(INFO) << "match size: " << match.size();
      std::unordered_map<VertexID, std::vector<int>> dis_origin;
      std::unordered_map<VertexID, std::unordered_set<VertexID>> ball;
      for (auto &w : match) {
        dhop_direct(dgraph, w, d_Q, ball[w], dis_origin[w]);
      }
      for (int extend = 1; extend <= extend_num; extend++) {
        std::unordered_set<std::pair<VertexID, VertexID>> add_edges, rm_edges;
        std::vector<std::pair<VertexID, VertexLabel>> add_vertices;
        Load_bunch_edges(add_vertices, add_edges, this->add_name, extend);
//        Load_bunch_edges(rm_edges, this->rm_name, extend);
        for (auto &v : add_vertices) {
          dgraph.AddVertex(Vertex(v.first, v.second));
        }
        for (auto &e : add_edges) {
          dgraph.AddEdge(Edge(e.first, e.second, 1));
        }
        dgraph.RebuildGraphProperties();
        double DIR_TIME = 0;
        double INC_TIME = 0;
        double INC_ORIGIN_TIME = 0;
        for (auto &w : match) {
          std::unordered_set<VertexID> ball_node1, ball_node2, ball_node3;
          std::vector<int> dis1, dis2, dis3;
          dis2.resize(dgraph.GetNumVertices(), max);
          dis3.resize(dgraph.GetNumVertices(), max);
          for (int i = 0; i < dis_origin[w].size(); i++) {
            dis2[i] = dis_origin[w][i];
            dis3[i] = dis_origin[w][i];
          }
          ball_node2 = ball[w];
          ball_node3 = ball[w];
          double starttime = get_current_time();
          dhop_direct(dgraph, w, d_Q, ball_node1, dis1);
          double stoptime = get_current_time();
          DIR_TIME += (stoptime - starttime); 
          starttime = get_current_time();
          dhop_add(dgraph, d_Q, ball_node2, dis2, add_edges);
          stoptime = get_current_time();
          INC_TIME += (stoptime - starttime);
          starttime = get_current_time();
          dhop_add_origin(dgraph, d_Q, ball_node3, dis3, add_edges);
          stoptime = get_current_time();
          INC_ORIGIN_TIME += (stoptime - starttime);
          if (!Verify(dis1 ,dis2)) {
            LOG(INFO) << "WARNING! We are under attack! Verify1";
          }
          if (!Verify(dis1 ,dis3)) {
            LOG(INFO) << "WARNING! We are under attack! Verify2";
          }
        }
        LOG(INFO) << "Direct Time: " << DIR_TIME << " Incremental Time: " << INC_TIME
                  << " Incremental Origin Time: " << INC_ORIGIN_TIME;
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
  std::string add_name = FLAGS_base_add_file;
  std::string rm_name = FLAGS_base_remove_file;
  int extend_num = FLAGS_extend_num;
  const int max = std::numeric_limits<int>::max() - 10;
};

int main(int argc, char **argv) {
  google::SetUsageMessage("Usage: dual_exp [gflags_opt]");
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::ShutDownCommandLineFlags();
  google::InitGoogleLogging("exp for dual incremental simulation");
  google::ShutdownGoogleLogging();
  DHop dual_inc_exp;
  dual_inc_exp.Run();
  return 0;  
}
