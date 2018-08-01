#include "dualsimulation.h"

DualSim::DualSim(){}

DualSim::~DualSim(){}

void DualSim::dual_sim_initialization(Graph &dgraph, Graph &qgraph,
                               std::vector<std::unordered_set<VertexID>> &sim,
                               bool &initialized_sim,
                               std::vector<std::unordered_set<VertexID>> &remove_pred,
                               std::vector<std::unordered_set<VertexID>> &remove_succ){
       // LOG(INFO)<<"dual_sim_initialzation"<<std::endl;
//        std::unordered_set<VertexID> pred_dgraph_vertices,succ_dgraph_vertices;
         auto dvnum = dgraph.GetNumVertices();
         LOG(INFO) <<dvnum;
         std::vector<bool> pred_dgraph_vertices(dvnum, false);
         std::vector<bool> succ_dgraph_vertices(dvnum, false);
        //initial pred_dgraph_vertices ,succ_dgraph_vertices
        for (auto v : dgraph.GetAllVerticesID()) {
            if (dgraph.GetOutDegree(v)!=0){
            //    pred_dgraph_vertices.insert(v);
               pred_dgraph_vertices[v] = true;
            }
            if (dgraph.GetInDegree(v)!=0){
            //    succ_dgraph_vertices.insert(v);
              succ_dgraph_vertices[v] = true;
            }
        }
       //initial sim unordered_set
       for(auto u : qgraph.GetAllVerticesID()){
            sim[u] = std::unordered_set<VertexID>();
            remove_pred[u] = std::unordered_set<VertexID>();
            remove_succ[u] = std::unordered_set<VertexID>();
            if (initialized_sim==false){
                if (qgraph.GetOutDegree(u) == 0 && qgraph.GetInDegree(u) ==0){
                    for (auto v : dgraph.GetAllVerticesID()){
                        if (qgraph.GetVertexLabel(u)==dgraph.GetVertexLabel(v)){
                            sim[u].insert(v);
                        }
                    }
                }
                else if (qgraph.GetOutDegree(u) != 0 && qgraph.GetInDegree(u) ==0){
                    for (auto v : dgraph.GetAllVerticesID()){
                        if (qgraph.GetVertexLabel(u)==dgraph.GetVertexLabel(v) && dgraph.GetOutDegree(v) != 0){
                            sim[u].insert(v);
                        }
                    }
                }
                else if (qgraph.GetOutDegree(u) == 0 && qgraph.GetInDegree(u) !=0){
                    for (auto v : dgraph.GetAllVerticesID()){
                        if (qgraph.GetVertexLabel(u)==dgraph.GetVertexLabel(v) && dgraph.GetInDegree(v) != 0){
                            sim[u].insert(v);
                        }
                    }
                }
                else {
                    for (auto v : dgraph.GetAllVerticesID()){
                        if (qgraph.GetVertexLabel(u)==dgraph.GetVertexLabel(v) && dgraph.GetOutDegree(v) != 0 && dgraph.GetInDegree(v) != 0){
                            sim[u].insert(v);
                         }
                      }
                 }
             }

        //inital remove_pred,remove_succ
        //    std::unordered_set<VertexID> pt1,pt2;
            std::vector<VertexID> pt1(dvnum, true);
            std::vector<VertexID> pt2(dvnum, true);
            for (auto w : sim[u]){
                for(auto v : dgraph.GetParentsID(w)){
            //        pt1.insert(v);
                   pt1[v] = false;
                }
                for (auto v : dgraph.GetChildrenID(w)){
            //        pt2.insert(v);
                   pt2[v] = false;
                }
            }

//            remove_pred[u] = diff(pred_dgraph_vertices,pt1);
//            remove_succ[u] = diff(succ_dgraph_vertices,pt2);
            for (auto v : dgraph.GetAllVerticesID()) {
//              if (pred_dgraph_vertices[v]) {
              if (pred_dgraph_vertices[v] & pt1[v]) {
//                if (pt1.find(v) != pt1.end())  LOG(INFO) << v;
                remove_pred[u].insert(v);
              }
//              if (succ_dgraph_vertices[v]) {
              if (succ_dgraph_vertices[v] & pt2[v]) {
//                if (pt2.find(v) != pt2.end()) LOG(INFO) << v;
                remove_succ[u].insert(v);
              } 
            }
       }
    }


void DualSim::reunordered_map_data_id(std::unordered_map<VertexID, VertexID> &t_f,Graph &dgraph){
        VertexID fid = 0;
        for (auto v : dgraph.GetAllVerticesID()){
            t_f[v] = fid;
            fid += 1;
        }
    }

void DualSim::dual_counter_initialization(Graph &dgraph, Graph &qgraph,
                                     std::vector<std::vector<int>> &sim_counter_post,
                                     std::vector<std::vector<int>> &sim_counter_pre,
                                     std::vector<std::unordered_set<VertexID>> &sim,
                                     std::unordered_map<VertexID, VertexID> &t_f){
        for (auto w : dgraph.GetAllVerticesID()){
            sim_counter_post[w] = std::vector<int>(qgraph.GetNumVertices(), 0);
            sim_counter_pre[w] = std::vector<int>(qgraph.GetNumVertices(), 0);
            for (auto u : qgraph.GetAllVerticesID()){
                int len_des=0,len_pre=0;
                for (auto des_w : dgraph.GetChildrenID(w)){
                     if(sim[u].find(des_w) != sim[u].end()){
                        len_des += 1;
                     }
                }
                for (auto pre_w : dgraph.GetParentsID(w)){
                    if (sim[u].find(pre_w) != sim[u].end()){
                        len_pre+=1;
                    }
                }
                sim_counter_post[w][u] = len_des;
                sim_counter_pre[w][u] = len_pre;
            }
        }
    }

VertexID DualSim::find_nonempty_remove(Graph &qgraph,
                           std::vector<std::unordered_set<VertexID>> &remove_pred,
                           std::vector<std::unordered_set<VertexID>> &remove_succ){
        for (auto u : qgraph.GetAllVerticesID())
        {
            if(remove_pred[u].size() !=0){
                return u;
            }
            if (remove_succ[u].size() !=0){
                return u;
            }
        }
        return -1;
    }

void DualSim::update_sim_counter(Graph &dgraph,
                            std::vector<std::vector<int>> &sim_counter_post,
                            std::vector<std::vector<int>> &sim_counter_pre,
                            VertexID u,VertexID w,
                            std::unordered_map<VertexID, VertexID> &t_f){
    for (auto wp : dgraph.GetParentsID(w)){
        if (sim_counter_post[wp][u] > 0){
            --sim_counter_post[wp][u];
        }
    }
    for (auto ws : dgraph.GetChildrenID(w)){
        if (sim_counter_pre[ws][u] > 0){
            --sim_counter_pre[ws][u];
        }

    }
    }

bool DualSim::match_check(Graph &qgraph,std::vector<std::unordered_set<VertexID>> &sim){
        for (auto u : qgraph.GetAllVerticesID()){
            if (sim[u].size() == 0){
                return false;
            }
        }
        return true;
    }

bool DualSim::dual_sim_output(Graph &qgraph,std::vector<std::unordered_set<VertexID>> &sim){
        if (match_check(qgraph, sim) == false){
           for(auto u : qgraph.GetAllVerticesID()){
                sim[u].clear();
           }
        }

    }


void DualSim::dual_sim_refinement(Graph &dgraph, Graph &qgraph,
                           std::vector<std::unordered_set<VertexID>> &sim,
                           std::vector<std::unordered_set<VertexID>> &remove_pred,
                           std::vector<std::unordered_set<VertexID>> &remove_succ){
        //LOG(INFO)<<"dual_sim_refinement"<<std::endl;
        auto dvnum = dgraph.GetNumVertices();
        std::vector<std::vector<int>> sim_counter_post(dvnum);
        std::vector<std::vector<int>> sim_counter_pre(dvnum);
        std::unordered_map<VertexID, VertexID> t_f;
//        reunordered_map_data_id(t_f, dgraph);

        dual_counter_initialization(dgraph, qgraph, sim_counter_post, sim_counter_pre, sim,t_f);
        int len_pre=0,len_succ=0;
        VertexID u= find_nonempty_remove(qgraph, remove_pred, remove_succ);
        while(u!=-1){
            if(remove_pred[u].size() !=0){
                for (auto u_p : qgraph.GetParentsID(u)){
                    for (auto w_pred : remove_pred[u]){
                        if (sim[u_p].find(w_pred)!=sim[u_p].end()){
                            sim[u_p].erase(w_pred);
                            len_pre++;
                            //LOG(INFO)<<u_p<<' '<<w_pred<<std::endl;
                            update_sim_counter(dgraph,sim_counter_post, sim_counter_pre, u_p, w_pred, t_f);
                            for (auto w_pp : dgraph.GetParentsID(w_pred)){
                                if (sim_counter_post[w_pp][u_p] == 0){
                                    remove_pred[u_p].insert(w_pp);
                                }
                            }
                            for(auto w_ps : dgraph.GetChildrenID(w_pred)){
                                if (sim_counter_pre[w_ps][u_p] == 0){
                                    remove_succ[u_p].insert(w_ps);
                                }
                            }
                        }
                    }
                }
                remove_pred[u].clear();
            }

            if(remove_succ[u].size() != 0){
                for (auto u_s : qgraph.GetChildrenID(u)){
                    for (auto w_succ : remove_succ[u]){
                        if (sim[u_s].find(w_succ) != sim[u_s].end()){
                            sim[u_s].erase(w_succ);
                            len_succ++;
                            //LOG(INFO)<<u_s<<' '<<w_succ<<std::endl;
                            update_sim_counter(dgraph,sim_counter_post, sim_counter_pre, u_s, w_succ, t_f);
                            for (auto w_sp : dgraph.GetParentsID(w_succ)){
                                if (sim_counter_post[w_sp][u_s] == 0){
                                    remove_pred[u_s].insert(w_sp);
                                }
                            }
                            for (auto w_ss : dgraph.GetChildrenID(w_succ)){
                                if (sim_counter_pre[w_ss][u_s] == 0){
                                    remove_succ[u_s].insert(w_ss);
                                }
                            }
                        }
                    }

                }
                remove_succ[u].clear();
            }

            u = find_nonempty_remove(qgraph, remove_pred, remove_succ);
        }
    }

void DualSim::dual_simulation(Graph &dgraph, Graph &qgraph, std::vector<std::unordered_set<VertexID>> &sim, bool &initialized_sim){
        //LOG(INFO)<<"begin dual"<<std::endl;
//        std::unordered_map<VertexID, std::unordered_set<VertexID>> remove_pred,remove_succ;
        auto dvnum = dgraph.GetNumVertices();
        std::vector<std::unordered_set<VertexID>> remove_pred(dvnum);
        std::vector<std::unordered_set<VertexID>> remove_succ(dvnum);
        dual_sim_initialization(dgraph, qgraph, sim, initialized_sim, remove_pred, remove_succ);
        dual_sim_refinement(dgraph, qgraph, sim, remove_pred, remove_succ);
        dual_sim_output(qgraph, sim);
  }

using namespace std;
