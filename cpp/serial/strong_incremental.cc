#include "strong_incremental.h"
#include "cpp/utils/time.h"
StrongInc::StrongInc(){}

StrongInc::~StrongInc(){}

void StrongInc::find_affected_center_area(Graph &dgraph,std::unordered_set<std::pair<VertexID,VertexID>> &add_edges,
                                                 std::unordered_set<std::pair<VertexID,VertexID>> &rm_edges,
                                                 int d_hop,
                                                 std::unordered_set<VertexID> &result){
    std::unordered_set<VertexID> incedges_node;
    for(auto e:add_edges){
        if(rm_edges.find(e) == rm_edges.end()){
			incedges_node.insert(e.first);
			incedges_node.insert(e.second);
        }
    }
    for(auto e:rm_edges){
        if(add_edges.find(e) == add_edges.end()){
			incedges_node.insert(e.first);
			incedges_node.insert(e.second);
        }
    }
    dgraph.find_hop_nodes(incedges_node,d_hop,result);
}

void StrongInc::find_affected_center_area(Graph &dgraph,std::unordered_set<std::pair<VertexID,VertexID>> &add_edges,
                                                 std::unordered_set<std::pair<VertexID,VertexID>> &rm_edges,
                                                 int d_hop,
                                                 std::unordered_set<VertexID> &result,
												 int flag){
	std::unordered_set<VertexID> incedges_node;
	if(flag == 0){
		for(auto e : rm_edges){
            if(add_edges.find(e) == add_edges.end()){
				// if(e.first < dgraph.GetNumVertices()){  // it's necessary when first remove edges.
					incedges_node.insert(e.first);
				// }
				// if(e.second < dgraph.GetNumVertices()){ // it's necessary when first remove edges.
					incedges_node.insert(e.second);
				// }
            }
        }		
	}
	else if (flag == 2){
		for(auto e:add_edges){
			if(rm_edges.find(e) == rm_edges.end()){
				incedges_node.insert(e.first);
				incedges_node.insert(e.second);
			}
        }
	}
    dgraph.find_hop_nodes(incedges_node,d_hop,result);
}

int StrongInc::cal_diameter_qgraph(Graph &qgraph){
          int temp_dia = 0;
          int max_dia = qgraph.GetNumVertices()-1;
          for(auto u : qgraph.GetAllVerticesID()){
              std::vector<int> dis;
              qgraph.shortest_distance(u,dis);
              for (int i=0; i<qgraph.GetNumVertices(); i++){
                if (dis[i] <= max_dia && temp_dia < dis[i]){
                    temp_dia = dis[i];
                }
              }
          }
          return temp_dia;
      }


bool StrongInc::valid_sim_w(Graph &qgraph,std::unordered_map<VertexID, std::unordered_set<VertexID>> &sim,VertexID w){
          for(auto u : qgraph.GetAllVerticesID()){
              if(sim[u].size()==0){
                  return false;
              }
           }
           int uid = -1;
           for(auto u : qgraph.GetAllVerticesID()){
               if (sim[u].find(w) != sim[u].end()){
                   uid = u;
                   return true;
               }
           }
           if (uid == -1){
           return false;
           }
      }


void StrongInc::find_node_connectivity_nodes(Ball_View &ball_view,std::unordered_set<VertexID> &v_set,VertexID w){
    std::queue<VertexID> q;
    v_set.clear();
    v_set.insert(w);
    q.push(w);
    while(!q.empty()){
        VertexID root = q.front();
        q.pop();
        for(auto v :ball_view.GetParentsID(root)){
            if(v_set.find(v) == v_set.end()){
                v_set.insert(v);
                q.push(v);
            }
        }
        for(auto v :ball_view.GetChildrenID(root)){
            if(v_set.find(v) == v_set.end()){
                v_set.insert(v);
                q.push(v);
            }
        }
    }
}

void StrongInc::rename_sim(Ball_View &ball_view,Graph &qgraph,
                               std::unordered_map<VertexID, std::unordered_set<VertexID>> &sim){
              //std::cout<<w<<std::endl;
       for(auto u : qgraph.GetAllVerticesID()){
           std::unordered_set<VertexID> tmp_set;
           for(auto v:ball_view.GetAllVerticesID()){
               if (sim[u].find(v) != sim[u].end()){
                 tmp_set.insert(v);
               }
           }
           sim[u].clear();
           sim[u]=tmp_set;
       }
     }

void StrongInc::ss_counter_initialization(Ball_View &ball_view,Graph &qgraph,
                                     std::unordered_map<VertexID, std::vector<int>> &sim_counter_pre,
                                     std::unordered_map<VertexID, std::vector<int>> &sim_counter_post,
                                     std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w){
        for (auto w : ball_view.GetAllVerticesID()){
            sim_counter_post[w] = std::vector<int>(qgraph.GetNumVertices(), 0);
            sim_counter_pre[w] = std::vector<int>(qgraph.GetNumVertices(), 0);
            for (auto u : qgraph.GetAllVerticesID()){
                int len_des=0,len_pre=0;
                for (auto des_w : ball_view.GetChildrenID(w)){
                     if(S_w[u].find(des_w) != S_w[u].end()){
                        len_des += 1;
                     }
                }
                for (auto pre_w : ball_view.GetParentsID(w)){
                    if (S_w[u].find(pre_w) != S_w[u].end()){
                        len_pre+=1;
                    }
                }
                sim_counter_post[w][u] = len_des;
                sim_counter_pre[w][u] = len_pre;
            }
        }
 }

 void StrongInc::dual_filter_match(Ball_View &refined_ball, Graph &qgraph,
                      std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w,VertexID w,int d_Q){
        std::set<std::pair<VertexID,VertexID> > filter_set;
        std::unordered_map<VertexID, std::vector<int> > sim_counter_pre,sim_counter_post;
        ss_counter_initialization(refined_ball,qgraph, sim_counter_pre,sim_counter_post,S_w);
        push_phase(refined_ball,qgraph,w,d_Q, filter_set, sim_counter_pre, sim_counter_post,S_w);
        decremental_refine(refined_ball,qgraph, filter_set,sim_counter_pre,sim_counter_post,S_w);
   }

void StrongInc::push_phase(Ball_View &ball,Graph &qgraph,VertexID w,int d_Q,
                          std::set<std::pair<VertexID,VertexID>> &filter_set,
                          std::unordered_map<VertexID, std::vector<int>> &sim_counter_pre,
                          std::unordered_map<VertexID, std::vector<int>> &sim_counter_post,
                          std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w){

         for (auto u :qgraph.GetAllVerticesID()){
            for (auto v : S_w[u]){
                for (auto u_s : qgraph.GetChildrenID(u)){
                    if (sim_counter_post[v][u_s]==0){
                       filter_set.insert(std::pair<VertexID,VertexID>(u,v));
                        break;
                    }
                }
                for(auto u_p : qgraph.GetParentsID(u)){
                    if(sim_counter_pre[v][u_p]==0){
                        filter_set.insert(std::pair<VertexID,VertexID>(u,v));
                       break;
                    }
                }
            }
        }
 }


 void StrongInc::update_counter(Ball_View &ball,Graph &qgraph,VertexID u,VertexID v,
                          std::unordered_map<VertexID, std::vector<int>> &sim_counter_pre,
                          std::unordered_map<VertexID, std::vector<int>> &sim_counter_post){
        for(auto vp : ball.GetParentsID(v)){
            if (sim_counter_post.find(vp)!=sim_counter_post.end()){
                if(sim_counter_post[vp][u]>0){
                    sim_counter_post[vp][u]-=1;
                }
            }
        }
        for (auto vs : ball.GetChildrenID(v)){
            if (sim_counter_pre.find(vs)!=sim_counter_pre.end()){
                if(sim_counter_pre[vs][u]>0){
                    sim_counter_pre[vs][u]-=1;
                }
            }
        }
    }

void StrongInc::decremental_refine(Ball_View &ball_view,Graph &qgraph,
                          std::set<std::pair<VertexID,VertexID>> &filter_set,
                          std::unordered_map<VertexID, std::vector<int>> &sim_counter_pre,
                          std::unordered_map<VertexID, std::vector<int>> &sim_counter_post,
                          std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w){
        while(!filter_set.empty()){
            std::pair<VertexID,VertexID> pmatch = *filter_set.begin();
            VertexID u = pmatch.first;
            VertexID v = pmatch.second;
            filter_set.erase(filter_set.begin());
            S_w[u].erase(v);
            update_counter(ball_view,qgraph,u,v,sim_counter_pre,sim_counter_post);
            for (auto u_p : qgraph.GetParentsID(u)){
                for (auto v_p : ball_view.GetParentsID(v)){
                    if (S_w[u_p].find(v_p)!=S_w[u_p].end()){
                        if(sim_counter_post[v_p][u]==0){
                            filter_set.insert(std::pair<VertexID,VertexID>(u_p,v_p));
                        }
                    }
                }
            }
            for(auto u_s :qgraph.GetChildrenID(u)){
                for(auto v_s: ball_view.GetChildrenID(v)){
                    if (S_w[u_s].find(v_s) != S_w[u_s].end()){
                        if(sim_counter_pre[v_s][u] == 0){
                            filter_set.insert(std::pair<VertexID,VertexID>(u_s,v_s));
                        }
                    }
                }

            }
        }
    }

void StrongInc::extract_max_pg(Ball_View &ball_view,Graph &dgraph,Graph &qgraph,VertexID w,
                                std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w){
    if(!valid_sim_w(qgraph,S_w,w)){
        for (auto u : qgraph.GetAllVerticesID()){
            S_w[u].clear();
        }
    }

    std::unordered_set<VertexID> vertex_match_set;
    for (auto u : qgraph.GetAllVerticesID()){
            for(auto v: S_w[u]){
             vertex_match_set.insert(v);
            }
     }

    std::unordered_set<Edge> edge_match_set;
    for(auto e: qgraph.GetAllEdges()){
        VertexID sourceid=e.src();
        VertexID targetid=e.dst();
        for (auto sim_v1 : S_w[sourceid]){
            for(auto sim_v2 : S_w[targetid]){
                 if (ball_view.ExistEdge(sim_v1,sim_v2)){
                     edge_match_set.insert(Edge(sim_v1,sim_v2,1));
                 }

             }
        }
   }
   Ball_View pg_view(vertex_match_set,edge_match_set);
   std::unordered_set<VertexID> vertex_match_set1;
   find_node_connectivity_nodes(pg_view,vertex_match_set1,w);


   for(auto u : qgraph.GetAllVerticesID()){
       std::unordered_set<VertexID> tmp_set;
       for(auto v:vertex_match_set1){
           if (S_w[u].find(v) != S_w[u].end()){
             tmp_set.insert(v);
           }
       }
       S_w[u].clear();
       S_w[u]=tmp_set;
       }
   }



void StrongInc::print_ball_info(Graph &qgraph,std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w,VertexID w){

              std::unordered_map<VertexID,std::set<VertexID>> printset;
              for(auto u :qgraph.GetAllVerticesID()){
                  printset[u]=std::set<VertexID>();
                  for(auto v:S_w[u]){
                      printset[u].insert(v);
                  }

              }
              for(auto u :qgraph.GetAllVerticesID()){
               std::cout<<u;
               for(auto v:printset[u]){
                   std::cout<<' '<<v;
               }
              std::cout<<std::endl;
              }
}

void  StrongInc::recalculate_incrementl_dual(Graph &dgraph, Graph &qgraph,
                                      std::unordered_map<VertexID,std::unordered_set<VertexID>> &dsim,
                                      std::unordered_set<std::pair<VertexID,VertexID>> &add_edges,
                                      std::unordered_set<std::pair<VertexID,VertexID>> &rm_edges){
          DualInc dualinc;
//          clock_t start1 = clock();
//          for (auto e:add_edges){
//             dgraph.AddEdge(Edge(e.first,e.second,1));
//          }
//          clock_t end1 = clock();
		  
//	      clock_t s1 = clock();
          dualinc.incremental_addedges(dgraph,qgraph,dsim,add_edges);
//          clock_t e1 = clock();
          
//          clock_t start2 = clock();
//          for(auto e :rm_edges){
//              dgraph.RemoveEdge(Edge(e.first,e.second,1));
//          }
//		  std::cout<<"Inc.........After add edges : "<<std::endl;
//          clock_t end2 = clock();
//		  dgraph.printGraphInfo();

//          clock_t s2 = clock();
          dualinc.incremental_removeedgs(dgraph,qgraph,dsim,rm_edges);
//          clock_t e2 = clock();

//          std::fstream outfile("time_info_add_rm.txt",std::ios::app);
//          outfile<<(float)(end2-start2+end1-start1)/CLOCKS_PER_SEC<<std::endl;
//          outfile.close();          
  
//          std::fstream outfile2("whole_time_info_add_rm.txt",std::ios::app);
//          outfile2<<"inc strong sim.......dualinc.incremental_addedges and incremental_removeedgs....time = "<<(float)(e2-s2+e1-s1)/CLOCKS_PER_SEC<<"s"<<std::endl;
//          outfile2.close();

   }

void StrongInc::cal_culculate_inc_dhop_nodes_add(Graph &dgraph, int d_Q,
                std::unordered_set<VertexID> &ball_node,
                std::vector<int> &dis,
                std::unordered_set<std::pair<VertexID,VertexID>> &add_edges){
    std::priority_queue<std::pair<int, VertexID>> source;
    for (auto &e : add_edges) {
      if ((dis[e.first] > dis[e.second] + 1) && dis[e.second] < d_Q) {
        dis[e.first] = dis[e.second] + 1;
        ball_node.insert(e.first);
        if (dis[e.first] < d_Q) {
          source.push(make_pair(-dis[e.first], e.first));
        }
      } else if ((dis[e.second] > dis[e.first] + 1) && dis[e.first] < d_Q) {
        dis[e.second] = dis[e.first] + 1;
        ball_node.insert(e.second);
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
          ball_node.insert(v);
          if (dis[v] < d_Q) {
            source.push(make_pair(-dis[v], v));
          }
        }
      }
      for (auto v : dgraph.GetParentsID(u)) {
        if (dis[v] > dis[u] + 1) {
          dis[v] = dis[u] + 1;
          ball_node.insert(v);
          if (dis[v] < d_Q) {
            source.push(make_pair(-dis[v], v));
          }
        }
      }
    }

}

void StrongInc::strong_simulation_inc_only_add(Graph &dgraph, Graph &qgraph,
                                      std::unordered_map<VertexID,std::unordered_set<VertexID>> &dsim,
                                      std::vector<StrongR> &strong_r,
                                      std::unordered_set<std::pair<VertexID,VertexID>> &add_edges,
									  int flag,
									  std::unordered_map<VertexID,std::unordered_set<VertexID>> &whole_ball_nodes,
									  std::unordered_map<VertexID,std::vector<int>> &whole_dist){
    /**
	*calculate qgaraph diameter
	*/
	std::vector<StrongR> max_result;
	clock_t s1 = clock();
	int d_Q = cal_diameter_qgraph(qgraph);
	clock_t e1 = clock();
		
	std::unordered_map<VertexID, std::unordered_set<VertexID>> global_sim;
		
	std::unordered_set<VertexID> affected_center_nodes;
	std::unordered_set<VertexID> affected_center_nodes_rm;
	std::unordered_set<VertexID> affected_center_nodes_add;
	DualInc dualinc;
	clock_t s4,e4,s5,e5,start2,end2;
	if(flag == 1){ // only add.
		start2 = clock();
		int max = dgraph.GetNumVertices() - 1;
		for(auto e : add_edges){
			if(e.first > max){
				max = e.first;
			}
			if(e.second > max){
				max = e.second;
			}
		}
		int tmp = dgraph.GetNumVertices();
		for(int i=0; i<=(max-tmp); i++){
			dgraph.AddVertex(Vertex(i, (i+123456789)%MAX_LABEL)); // save it.
		}
		for(auto e : add_edges){
			dgraph.AddEdge(Edge(e.first,e.second,1));
		}
		dgraph.RebuildGraphProperties();
		// dgraph.printGraphInfo();
		end2 = clock();
		std::fstream outfile3("time_info_only_add.txt",std::ios::app);
		outfile3<<(float)(end2-start2)/CLOCKS_PER_SEC<<std::endl;
		outfile3.close();
		
		s4 = clock();
		dualinc.incremental_addedges(dgraph,qgraph,dsim,add_edges);
		e4 = clock();
		
		s5 = clock();
		std::unordered_set<std::pair<VertexID,VertexID>> rm_edges;
		find_affected_center_area(dgraph,add_edges,rm_edges,d_Q,affected_center_nodes,2); //在更新后的图上，找增边的影响区域
		e5 = clock();
	}
	
    std::unordered_set<VertexID> max_dual_set;
    for(auto u : qgraph.GetAllVerticesID()){
		for(auto v : dsim[u]){
			max_dual_set.insert(v);
        }
    }
    clock_t s7 = clock();
    affected_center_nodes = intersection(affected_center_nodes,max_dual_set);
    clock_t e7 = clock();
	
    int i = 0;
    clock_t stime = clock();
    for (auto w : max_dual_set) {
		/**
		* calculate ball for center w if w if a valid center
		*/
		// if (valid_sim_w(qgraph,dsim,w)){
		if (affected_center_nodes.find(w) == affected_center_nodes.end()){
			for(auto strong_ball:strong_r){
				if(strong_ball.center()==w){
					max_result.push_back(strong_ball);
				}
            }
			continue;
        }
		clock_t start = clock();
		
		/**
		*find d_hop_nodes for w in dgraph
		*/
		
		std::unordered_set<VertexID> ball_node;
		int tmp_flag = 0;
		if(whole_ball_nodes.find(w) == whole_ball_nodes.end()){
			dgraph.find_hop_nodes(w,d_Q,ball_node);
			tmp_flag = 1;
		}
		else{
			cal_culculate_inc_dhop_nodes_add(dgraph, d_Q, whole_ball_nodes[w], whole_dist[w], add_edges);
		}
		
		std::unordered_set<VertexID> ball_filter_node;
        std::unordered_set<Edge> ball_filter_edge;
        std::unordered_map<VertexID, std::unordered_set<VertexID>> S_w;
		if(tmp_flag == 1){
			for(auto u : qgraph.GetAllVerticesID()){
				for (auto v : dsim[u]){
					if(ball_node.find(v) != ball_node.end()){
						S_w[u].insert(v);
						ball_filter_node.insert(v);
					}
				}
			}
		}
		else{
			for(auto u : qgraph.GetAllVerticesID()){
				for (auto v : dsim[u]){
					if(whole_ball_nodes[w].find(v) != whole_ball_nodes[w].end()){
						S_w[u].insert(v);
						ball_filter_node.insert(v);
					}
				}
			}
		}
		
        for(auto e : qgraph.GetAllEdges()){
            VertexID sourceid=e.src();
            VertexID targetid=e.dst();
            for (auto sim_v1 : S_w[sourceid]){
                for(auto sim_v2 : S_w[targetid]){
                    if (dgraph.ExistEdge(sim_v1,sim_v2)){
                        ball_filter_edge.insert(Edge(sim_v1,sim_v2,1));
                    }
                }
            }
        }
		/*
		std::unordered_set<VertexID> ball_node;
		int tmp_size=0;
        if(whole_ball_nodes.find(w)!=whole_ball_nodes.end()){
			tmp_size=whole_ball_nodes[w].size();
            cal_culculate_inc_dhop_nodes_add(dgraph,d_Q,whole_ball_nodes[w],whole_dist[w],add_edges);
            for(auto v:whole_ball_nodes[w]){
                ball_node.insert(v);
            }
        }else{
            dgraph.find_hop_nodes(w,d_Q,ball_node);
        }

        std::unordered_set<VertexID> ball_filter_node;
        std::unordered_set<Edge> ball_filter_edge;
        std::unordered_map<VertexID, std::unordered_set<VertexID>> S_w;
        for(auto u : qgraph.GetAllVerticesID()){
            for (auto v : dsim[u]){
				if(ball_node.find(v) != ball_node.end()){
					S_w[u].insert(v);
					ball_filter_node.insert(v);
				}	
			}
        }
		if(tmp_size==ball_filter_node.size()){
            for(auto strong_ball:strong_r){
                if(strong_ball.center()==w){
                    max_result.push_back(strong_ball);
                }
            }
            continue;				     
		}
        for(auto e: qgraph.GetAllEdges()){
            VertexID sourceid=e.src();
            VertexID targetid=e.dst();
            for (auto sim_v1 : S_w[sourceid]){
                for(auto sim_v2 : S_w[targetid]){
                    if (dgraph.ExistEdge(sim_v1,sim_v2)){
                        ball_filter_edge.insert(Edge(sim_v1,sim_v2,1));
                    }
                }
            }
        }
		*/
		
        Ball_View ball_view(ball_filter_node,ball_filter_edge);

        std::unordered_set<VertexID> refined_ball_vertex;
        std::unordered_set<Edge> refinded_ball_edge;
        find_node_connectivity_nodes(ball_view,refined_ball_vertex,w);
        for(auto e :ball_filter_edge){
            if(refined_ball_vertex.find(e.src()) != refined_ball_vertex.end() && refined_ball_vertex.find(e.dst())!=refined_ball_vertex.end()){
                refinded_ball_edge.insert(e);
            }
        }
        Ball_View refined_ball_view(refined_ball_vertex,refinded_ball_edge);
        rename_sim(refined_ball_view,qgraph,S_w);
        dual_filter_match(refined_ball_view, qgraph,S_w,w,d_Q);

        extract_max_pg(refined_ball_view,dgraph,qgraph, w,S_w);

        max_result.emplace_back(w,S_w);
        // print_ball_info(qgraph,S_w,w);
        // break;
	    clock_t end = clock();
        std::fstream out_("calculate_ball_info_inc.txt",std::ios::app);
        out_<<"calculate one ball time "<<(float)(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
        out_.close();
        // std::cout<<"calculate one ball time "<<(float)(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
    }
    clock_t etime = clock();
    strong_r.clear();
    for(auto strong_ball :max_result){
        strong_r.push_back(strong_ball);
    }
	
	std::fstream tmp_outfile("whole_time_info_only_add.txt",std::ios::app);
	tmp_outfile<<"inc strong sim...calculate d_Q = "<<(float)(e1-s1)/CLOCKS_PER_SEC<<"s"<<std::endl;
	tmp_outfile<<"inc strong sim....dualinc = "<<(float)(e4-s4)/CLOCKS_PER_SEC<<"s"<<std::endl;
	tmp_outfile<<"inc strong sim......find_affected_center_area = "<<(float)(e5-s5)/CLOCKS_PER_SEC<<"s"<<std::endl;
	tmp_outfile<<"inc strong sim........find affected_center_nodes by intersection = "<<(float)(e7-s7)/CLOCKS_PER_SEC<<"s"<<std::endl;
	tmp_outfile<<"inc strong sim..........calculate all ball =  "<<(float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
	tmp_outfile<<"strongr.size = "<<strong_r.size()<<std::endl;
	tmp_outfile<<"------------------------------------------------------------------------------------------"<<std::endl;
	tmp_outfile.close();

    // std::cout<<"inc strong "<< (float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
}

void StrongInc::strong_simulation_inc(Graph &dgraph, Graph &qgraph,
                std::unordered_map<VertexID,std::unordered_set<VertexID>> &dsim,
                std::vector<StrongR> &strong_r,
                std::unordered_set<std::pair<VertexID,VertexID>> &add_edges,
                std::unordered_set<std::pair<VertexID,VertexID>> &rm_edges,
                std::unordered_map<VertexID, std::unordered_set<VertexID>> &whole_ball_nodes,
                std::unordered_map<VertexID, std::vector<int>> &whole_dist){
          /**
           *calculate qgaraph diameter
          */
          std::vector<StrongR> max_result;
          int d_Q = cal_diameter_qgraph(qgraph);
          std::unordered_map<VertexID, std::unordered_set<VertexID>> global_sim;
          double start = get_current_time();
//          recalculate_incrementl_dual(dgraph,qgraph,dsim,add_edges,rm_edges);
          DualInc dualinc;
          dualinc.incremental_addedges(dgraph, qgraph, dsim, add_edges);
          dualinc.incremental_removeedgs(dgraph, qgraph, dsim, rm_edges);
          double end = get_current_time();
          double time = end - start;
          LOG(INFO) << "incremental dual time: " << time;

          std::unordered_set<VertexID> affected_center_nodes;
          start = get_current_time();
          find_affected_center_area(dgraph,add_edges,rm_edges,d_Q,affected_center_nodes);
          end = get_current_time();
          time = end - start;
          LOG(INFO) << "find affected center area time: " << time;

          start = get_current_time();
          std::unordered_set<VertexID> max_dual_set;
          for(auto u:qgraph.GetAllVerticesID()){
              for(auto v :dsim[u]){
                  max_dual_set.insert(v);
              }
          }
          end = get_current_time();
          LOG(INFO) << "max dual set time: " << end - start;

          start = get_current_time();
          affected_center_nodes = intersection(affected_center_nodes,max_dual_set);
          end = get_current_time();
          LOG(INFO) << "intersection: " << end - start;

          int i=0;

          start = get_current_time();
          vector<double> parttime(11, 0.0);
          for (auto w : max_dual_set) {
              /**
               * calculate ball for center w if w if a valid center
               */
//              if (valid_sim_w(qgraph,dsim,w)){
              if (affected_center_nodes.find(w) == affected_center_nodes.end()){
                double partstart = get_current_time();
                 for(auto strong_ball:strong_r){
                     if(strong_ball.center()==w){
                          max_result.push_back(strong_ball);
                     }
                 }
                double partend = get_current_time();
                parttime[0] = partend - partstart;
                  continue;
              }
              /**
               *find d_hop_nodes for w in dgraph
               */
              std::unordered_set<VertexID> ball_node;
              bool already_ball_node;
              double partstart = get_current_time();
              if (whole_ball_nodes.find(w) == whole_ball_nodes.end()) {
                dgraph.find_hop_nodes(w, d_Q, ball_node);
                already_ball_node = false; 
              } else {
                cal_culculate_inc_dhop_nodes_add(dgraph, d_Q, whole_ball_nodes[w],
                       whole_dist[w], add_edges);    
                already_ball_node = true;
              }
              double partend = get_current_time();
              parttime[1] += (partend - partstart);
              
              std::unordered_set<VertexID> ball_filter_node;
              std::unordered_set<Edge> ball_filter_edge;
              std::unordered_map<VertexID, std::unordered_set<VertexID>> S_w;

              partstart = get_current_time();
              if (!already_ball_node) {
                for(auto u : qgraph.GetAllVerticesID()){
                    for (auto v : dsim[u]){
                        if(ball_node.find(v) != ball_node.end()){
                            S_w[u].insert(v);
                            ball_filter_node.insert(v);
                        }
                    }
                }
              } else {
                for (auto u : qgraph.GetAllVerticesID()) {
                  for (auto &v : dsim[u]) {
                    if (whole_ball_nodes[w].find(v) != whole_ball_nodes[w].end()) {
                      S_w[u].insert(v);
                      ball_filter_node.insert(v);
                    }
                  }
                }
              }
              partend = get_current_time();
              parttime[2] += (partend - partstart);
 
              partstart = get_current_time();
              for(auto e: qgraph.GetAllEdges()){
                  VertexID sourceid=e.src();
                  VertexID targetid=e.dst();
                  for (auto sim_v1 : S_w[sourceid]){
                      for(auto sim_v2 : S_w[targetid]){
                          if (dgraph.ExistEdge(sim_v1,sim_v2)){
                              ball_filter_edge.insert(Edge(sim_v1,sim_v2,1));
                          }
                      }
                  }
              }
              partend = get_current_time();
              parttime[3] += (partend - partstart);

              partstart = get_current_time();
              Ball_View ball_view(ball_filter_node,ball_filter_edge);
              partend = get_current_time();
              parttime[4] += (partend - partstart);

              std::unordered_set<VertexID> refined_ball_vertex;
              std::unordered_set<Edge> refinded_ball_edge;
              find_node_connectivity_nodes(ball_view,refined_ball_vertex,w);
           
              partstart = get_current_time();
              for(auto e :ball_filter_edge){
                   if(refined_ball_vertex.find(e.src()) != refined_ball_vertex.end() && refined_ball_vertex.find(e.dst())!=refined_ball_vertex.end()){
                       refinded_ball_edge.insert(e);
                   }
              }
              partend = get_current_time();
              parttime[5] += (partend - partstart);

              partstart = get_current_time();
              Ball_View refined_ball_view(refined_ball_vertex,refinded_ball_edge);
              partend = get_current_time();
              parttime[6] += (partend - partstart);
             
              partstart = get_current_time();
              rename_sim(refined_ball_view,qgraph,S_w);
              partend = get_current_time();
              parttime[7] += (partend - partstart);

              partstart = get_current_time();
              dual_filter_match(refined_ball_view, qgraph,S_w,w,d_Q);
              partend = get_current_time();
              parttime[8] += (partend - partstart);

              partstart = get_current_time();
              extract_max_pg(refined_ball_view,dgraph,qgraph, w,S_w);
              partend = get_current_time();
              parttime[9] += (partend - partstart);

              partstart = get_current_time();
              max_result.emplace_back(w,S_w);
              partend = get_current_time();
              parttime[10] += (partend - partstart);

//              print_ball_info(qgraph,S_w,w);
//              break;
             // std::cout<<"calculate one ball time "<<(float)(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
//              }
          }
          end  = get_current_time();
          time = end - start;
          LOG(INFO) << "main calculate time: " << time; 
          LOG(INFO) << "    --max result time: " << parttime[0];
          LOG(INFO) << "    --dhop time: " << parttime[1];
          LOG(INFO) << "    --ball filter node time: " << parttime[2];
          LOG(INFO) << "    --ball filter edge time: " << parttime[3];
          LOG(INFO) << "    --ball view time: " << parttime[4];
          LOG(INFO) << "    --refined ball edge time: " << parttime[5];
          LOG(INFO) << "    --refined ball view time: " << parttime[6];
          LOG(INFO) << "    --rename sim time: " << parttime[7];
          LOG(INFO) << "    --dual filter match time: " << parttime[8];
          LOG(INFO) << "    --extract max pg time: " << parttime[9];
          LOG(INFO) << "    --max result time: " << parttime[10];
          LOG(INFO) << "==============================================" ;
          strong_r.clear();
          for(auto strong_ball :max_result){
              strong_r.push_back(strong_ball);
          }
         // std::cout<<"inc strong "<< (float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
//          return max_result;
      }
   
void StrongInc::strong_simulation_inc(Graph &dgraph, Graph &qgraph,
                                      std::unordered_map<VertexID,std::unordered_set<VertexID>> &dsim,
                                      std::vector<StrongR> &strong_r,
                                      std::unordered_set<std::pair<VertexID,VertexID>> &add_edges,
                                      std::unordered_set<std::pair<VertexID,VertexID>> &rm_edges,
									  int flag){
    /**
	*calculate qgaraph diameter
	*/
	std::vector<StrongR> max_result;
	clock_t s1 = clock();
	int d_Q = cal_diameter_qgraph(qgraph);
	clock_t e1 = clock();
		
	std::unordered_map<VertexID, std::unordered_set<VertexID>> global_sim;
	
//     recalculate_incrementl_dual(dgraph,qgraph,dsim,add_edges,rm_edges);
	
	std::unordered_set<VertexID> affected_center_nodes;
	std::unordered_set<VertexID> affected_center_nodes_rm;
	std::unordered_set<VertexID> affected_center_nodes_add;
	DualInc dualinc;
	clock_t s2,e2,s3,e3,s4,e4,s5,e5,s6,e6,start1,end1,start2,end2;
	
	if(flag == 2){ //first remove, then add.
		// first remove edges
//		s2 = clock();
		find_affected_center_area(dgraph,add_edges,rm_edges,d_Q,affected_center_nodes_rm,0); //在原图上，找减边的影响区域
//		e2 = clock();
		
//		start1 = clock();
		for(auto e : rm_edges){
			dgraph.RemoveEdge(Edge(e.first,e.second,1));
		}
		dgraph.RebuildGraphProperties();
		// dgraph.printGraphInfo();
//		end1 = clock();
		
//		s3 = clock();
		dualinc.incremental_removeedgs(dgraph,qgraph,dsim,rm_edges);
//		e3 = clock();
		
		// then add edges
//		start2 = clock();
		int max = dgraph.GetNumVertices() - 1;
		for(auto e : add_edges){
			if(e.first > max){
				max = e.first;
			}
			if(e.second > max){
				max = e.second;
			}
		}
		int tmp = dgraph.GetNumVertices();
		for(int i=0; i<=(max-tmp); i++){
			dgraph.AddVertex(Vertex(i, (i+123456789)%MAX_LABEL)); // save it.
		}
		for(auto e : add_edges){
			dgraph.AddEdge(Edge(e.first,e.second,1));
		}
		dgraph.RebuildGraphProperties();
		// dgraph.printGraphInfo();
//		end2 = clock();
//		std::fstream outfile3("time_info_rm_and_add.txt",std::ios::app);
//		outfile3<<(float)(end2-start2+end1-start1)/CLOCKS_PER_SEC<<std::endl;
//		outfile3.close();
		
//		s4 = clock();
		dualinc.incremental_addedges(dgraph,qgraph,dsim,add_edges);
//		e4 = clock();
//		s5 = clock();
		find_affected_center_area(dgraph,add_edges,rm_edges,d_Q,affected_center_nodes_add,2); //在更新之后的图上，找增边的影响区域
//		e5 = clock();

		// union
//		s6 = clock();
		affected_center_nodes = unions(affected_center_nodes_rm,affected_center_nodes_add);
//		e6 = clock();
	}
	else if(flag == 0){ // only remove.
		s2 = clock();
		find_affected_center_area(dgraph,add_edges,rm_edges,d_Q,affected_center_nodes,0); //在原图上，找减边的影响区域
		e2 = clock();
		
		start1 = clock();
		for(auto e : rm_edges){
			dgraph.RemoveEdge(Edge(e.first,e.second,1));
		}
		dgraph.RebuildGraphProperties();
		// dgraph.printGraphInfo();
		end1 = clock();
		std::fstream outfile3("time_info_only_rm.txt",std::ios::app);
		outfile3<<(float)(end1-start1)/CLOCKS_PER_SEC<<std::endl;
		outfile3.close();
		
		s3 = clock();
		dualinc.incremental_removeedgs(dgraph,qgraph,dsim,rm_edges);
		e3 = clock();	
	}
	else if(flag == 1){ // only add.
		start2 = clock();
		int max = dgraph.GetNumVertices() - 1;
		for(auto e : add_edges){
			if(e.first > max){
				max = e.first;
			}
			if(e.second > max){
				max = e.second;
			}
		}
		int tmp = dgraph.GetNumVertices();
		for(int i=0; i<=(max-tmp); i++){
			dgraph.AddVertex(Vertex(i, (i+123456789)%MAX_LABEL)); // save it.
		}
		for(auto e : add_edges){
			dgraph.AddEdge(Edge(e.first,e.second,1));
		}
		dgraph.RebuildGraphProperties();
		// dgraph.printGraphInfo();
		end2 = clock();
		std::fstream outfile3("time_info_only_add.txt",std::ios::app);
		outfile3<<(float)(end2-start2)/CLOCKS_PER_SEC<<std::endl;
		outfile3.close();
		
		s4 = clock();
		dualinc.incremental_addedges(dgraph,qgraph,dsim,add_edges);
		e4 = clock();
		
		s5 = clock();
		find_affected_center_area(dgraph,add_edges,rm_edges,d_Q,affected_center_nodes,2); //在更新后的图上，找增边的影响区域
		e5 = clock();
	}
	else if(flag == 3){ // first add, then remove.
		// first add.
		start2 = clock();
		int max = dgraph.GetNumVertices() - 1;
		for(auto e : add_edges){
			if(e.first > max){
				max = e.first;
			}
			if(e.second > max){
				max = e.second;
			}
		}
		int tmp = dgraph.GetNumVertices();
		for(int i=0; i<=(max-tmp); i++){
			dgraph.AddVertex(Vertex(i, (i+123456789)%MAX_LABEL)); // save it.
		}
		for(auto e : add_edges){
			dgraph.AddEdge(Edge(e.first,e.second,1));
		}
		dgraph.RebuildGraphProperties();
		// dgraph.printGraphInfo();
		end2 = clock();
		
		s4 = clock();
		dualinc.incremental_addedges(dgraph,qgraph,dsim,add_edges);
		e4 = clock();
		s5 = clock();
		find_affected_center_area(dgraph,add_edges,rm_edges,d_Q,affected_center_nodes_add,2); //在更新之后的图上，找增边的影响区域
		e5 = clock();
		
		// then remove edges
		s2 = clock();
		find_affected_center_area(dgraph,add_edges,rm_edges,d_Q,affected_center_nodes_rm,0); //找减边的影响区域
		e2 = clock();
		
		start1 = clock();
		for(auto e : rm_edges){
			dgraph.RemoveEdge(Edge(e.first,e.second,1));
		}
		dgraph.RebuildGraphProperties();
		// dgraph.printGraphInfo();
		end1 = clock();
		
		s3 = clock();
		dualinc.incremental_removeedgs(dgraph,qgraph,dsim,rm_edges);
		e3 = clock();
		
		// union
		s6 = clock();
		affected_center_nodes = unions(affected_center_nodes_rm,affected_center_nodes_add);
		e6 = clock();
		
		std::fstream outfile3("time_info_add_and_rm.txt",std::ios::app);
		outfile3<<(float)(end2-start2+end1-start1)/CLOCKS_PER_SEC<<std::endl;
		outfile3.close();
	}
	
    std::unordered_set<VertexID> max_dual_set;
    for(auto u : qgraph.GetAllVerticesID()){
		for(auto v : dsim[u]){
			max_dual_set.insert(v);
        }
    }
    clock_t s7 = clock();
    affected_center_nodes = intersection(affected_center_nodes,max_dual_set);
    clock_t e7 = clock();
	
    int i = 0;
    clock_t stime = clock();
    for (auto w : max_dual_set) {
		/**
		* calculate ball for center w if w if a valid center
		*/
		// if (valid_sim_w(qgraph,dsim,w)){
		if (affected_center_nodes.find(w) == affected_center_nodes.end()){
			for(auto strong_ball:strong_r){
				if(strong_ball.center()==w){
					max_result.push_back(strong_ball);
				}
            }
			continue;
        }
		clock_t start = clock();
		/**
		*find d_hop_nodes for w in dgraph
		*/
		std::unordered_set<VertexID> ball_node;
		dgraph.find_hop_nodes(w,d_Q,ball_node);
		
		std::unordered_set<VertexID> ball_filter_node;
        std::unordered_set<Edge> ball_filter_edge;
        std::unordered_map<VertexID, std::unordered_set<VertexID>> S_w;
        for(auto u : qgraph.GetAllVerticesID()){
            for (auto v : dsim[u]){
                if(ball_node.find(v) != ball_node.end()){
                    S_w[u].insert(v);
                    ball_filter_node.insert(v);
                }
            }
        }
        for(auto e : qgraph.GetAllEdges()){
            VertexID sourceid=e.src();
            VertexID targetid=e.dst();
            for (auto sim_v1 : S_w[sourceid]){
                for(auto sim_v2 : S_w[targetid]){
                    if (dgraph.ExistEdge(sim_v1,sim_v2)){
                        ball_filter_edge.insert(Edge(sim_v1,sim_v2,1));
                    }
                }
            }
        }

        Ball_View ball_view(ball_filter_node,ball_filter_edge);

        std::unordered_set<VertexID> refined_ball_vertex;
        std::unordered_set<Edge> refinded_ball_edge;
        find_node_connectivity_nodes(ball_view,refined_ball_vertex,w);
        for(auto e :ball_filter_edge){
            if(refined_ball_vertex.find(e.src()) != refined_ball_vertex.end() && refined_ball_vertex.find(e.dst())!=refined_ball_vertex.end()){
                refinded_ball_edge.insert(e);
            }
        }
        Ball_View refined_ball_view(refined_ball_vertex,refinded_ball_edge);
        rename_sim(refined_ball_view,qgraph,S_w);
        dual_filter_match(refined_ball_view, qgraph,S_w,w,d_Q);

        extract_max_pg(refined_ball_view,dgraph,qgraph, w,S_w);

        max_result.emplace_back(w,S_w);
        // print_ball_info(qgraph,S_w,w);
        // break;
	    clock_t end = clock();
        std::fstream out_("calculate_ball_info_inc.txt",std::ios::app);
        out_<<"calculate one ball time "<<(float)(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
        out_.close();
        // std::cout<<"calculate one ball time "<<(float)(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
    }
    clock_t etime = clock();
    strong_r.clear();
    for(auto strong_ball :max_result){
        strong_r.push_back(strong_ball);
    }

	if(flag == 2){ // first remove, then add.
		std::fstream tmp_outfile("whole_time_info_rm_and_add.txt",std::ios::app);
		tmp_outfile<<"inc strong sim...calculate d_Q = "<<(float)(e1-s1)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim....dualinc = "<<(float)(e4-s4+e3-s3)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim......find_affected_center_area = "<<(float)(e5-s5+e2-s2)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim........find affected_center_nodes by unions and intersection = "<<(float)(e7-s7+e6-s6)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim..........calculate all ball =  "<<(float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
		tmp_outfile<<"strongr.size = "<<strong_r.size()<<std::endl;
		tmp_outfile<<"------------------------------------------------------------------------------------------"<<std::endl;
		tmp_outfile.close();
	}
	else if(flag == 0){ // only remove.
		std::fstream tmp_outfile("whole_time_info_only_rm.txt",std::ios::app);
		tmp_outfile<<"inc strong sim...calculate d_Q = "<<(float)(e1-s1)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim....dualinc = "<<(float)(e3-s3)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim......find_affected_center_area = "<<(float)(e2-s2)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim........find affected_center_nodes by intersection = "<<(float)(e7-s7)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim..........calculate all ball =  "<<(float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
		tmp_outfile<<"strongr.size = "<<strong_r.size()<<std::endl;
		tmp_outfile<<"------------------------------------------------------------------------------------------"<<std::endl;
		tmp_outfile.close();
	}
	else if(flag == 1){ // only add.
		std::fstream tmp_outfile("whole_time_info_only_add.txt",std::ios::app);
		tmp_outfile<<"inc strong sim...calculate d_Q = "<<(float)(e1-s1)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim....dualinc = "<<(float)(e4-s4)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim......find_affected_center_area = "<<(float)(e5-s5)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim........find affected_center_nodes by intersection = "<<(float)(e7-s7)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim..........calculate all ball =  "<<(float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
		tmp_outfile<<"strongr.size = "<<strong_r.size()<<std::endl;
		tmp_outfile<<"------------------------------------------------------------------------------------------"<<std::endl;
		tmp_outfile.close();
	}
	else if(flag == 3){ // first add, then remove.
		std::fstream tmp_outfile("whole_time_info_add_and_rm.txt",std::ios::app);
		tmp_outfile<<"inc strong sim...calculate d_Q = "<<(float)(e1-s1)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim....dualinc = "<<(float)(e4-s4+e3-s3)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim......find_affected_center_area = "<<(float)(e5-s5+e2-s2)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim........find affected_center_nodes by unions and intersection = "<<(float)(e7-s7+e6-s6)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"inc strong sim..........calculate all ball =  "<<(float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
		tmp_outfile<<"strongr.size = "<<strong_r.size()<<std::endl;
		tmp_outfile<<"------------------------------------------------------------------------------------------"<<std::endl;
		tmp_outfile.close();
	}
    // std::cout<<"inc strong "<< (float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
}


