#include "strongsimulation.h"
#include "cpp/utils/time.h"

StrongSim::StrongSim(){}

StrongSim::~StrongSim(){}

int StrongSim::cal_diameter_qgraph(Graph &qgraph){
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


bool StrongSim::valid_sim_w(Graph &qgraph,std::unordered_map<VertexID, std::unordered_set<VertexID>> &sim,VertexID w){
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


void StrongSim::find_node_connectivity_nodes(Ball_View &ball_view,std::unordered_set<VertexID> &v_set,VertexID w){
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

void StrongSim::rename_sim(Ball_View &ball_view,Graph &qgraph,
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

void StrongSim::ss_counter_initialization(Ball_View &ball_view,Graph &qgraph,
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

 void StrongSim::dual_filter_match(Ball_View &refined_ball, Graph &qgraph,
                      std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w,VertexID w,int d_Q){
        std::set<std::pair<VertexID,VertexID> > filter_set;
        std::unordered_map<VertexID, std::vector<int> > sim_counter_pre,sim_counter_post;
        ss_counter_initialization(refined_ball,qgraph, sim_counter_pre,sim_counter_post,S_w);
        push_phase(refined_ball,qgraph,w,d_Q, filter_set, sim_counter_pre, sim_counter_post,S_w);
        decremental_refine(refined_ball,qgraph, filter_set,sim_counter_pre,sim_counter_post,S_w);
   }

void StrongSim::push_phase(Ball_View &ball,Graph &qgraph,VertexID w,int d_Q,
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


 void StrongSim::update_counter(Ball_View &ball,Graph &qgraph,VertexID u,VertexID v,
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

void StrongSim::decremental_refine(Ball_View &ball_view,Graph &qgraph,
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

void StrongSim::extract_max_pg(Ball_View &ball_view,Graph &dgraph,Graph &qgraph,VertexID w,
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
//                     edge_match_set.insert(Edge(sim_v1,sim_v2,1));
                   edge_match_set.insert(Edge(sim_v1, sim_v2));
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

void StrongSim::print_ball_info(Graph &qgraph,std::unordered_map<VertexID, std::unordered_set<VertexID>> &S_w,VertexID w){
      std::unordered_map<VertexID,std::set<VertexID>> printset;
      for(auto u :qgraph.GetAllVerticesID()){
          printset[u]=std::set<VertexID>();
          for(auto v:S_w[u]){
              printset[u].insert(v);
          }
      }
      for(auto u :qgraph.GetAllVerticesID()){
          std::cout<<w<<' '<<u;
          for(auto v:printset[u]){
             std::cout<<' '<<v;
          }
          std::cout<<std::endl;
      }
}

std::vector<StrongR> StrongSim::strong_simulation_sim(Graph &dgraph, Graph &qgraph,
        std::unordered_map<VertexID, std::unordered_set<VertexID>> &whole_ball_nodes,
        std::unordered_map<VertexID, std::vector<int>> &whole_dist){
    /**
     *calculate qgaraph diameter
    */
    std::vector<StrongR> max_result;
        double start = get_current_time();
	int d_Q = cal_diameter_qgraph(qgraph);
        double end = get_current_time();
        double time = end-start;
        LOG(INFO) << "cal diameter qgraph time: " << time;
	
	std::unordered_map<VertexID, std::unordered_set<VertexID>> global_sim;
	/**
	*calculate dual simulation for dgraph
	*/
	DualSim dualsim;
	bool inital_sim = false;
	
        start = get_current_time();
	dualsim.dual_simulation(dgraph,qgraph,global_sim,inital_sim);
        end = get_current_time();
        time = end - start;
        LOG(INFO) << "dual simulation time: " << time;
	
	int i=0;
        int count = 0;
        start = get_current_time();
        vector<double> parttime(10, 0.0);
	for (auto w : dgraph.GetAllVerticesID()) {
	  /**
	  * calculate ball for center w if w if a valid center
	  */
	  if (valid_sim_w(qgraph,global_sim,w)){
	    clock_t start =clock();
	    /**
	    *find d_hop_nodes for w in dgraph
	    */
            count++;
            double partstart = get_current_time();
            cal_culculate_directed_dhop_nodes(dgraph, w, d_Q, whole_ball_nodes[w],
                 whole_dist[w]);
            double partend = get_current_time();
            parttime[0] += partend - partstart;
//	    std::unordered_set<VertexID> ball_node;
//	    dgraph.find_hop_nodes(w,d_Q,ball_node);
			
	    std::unordered_set<VertexID> ball_filter_node;
	    std::unordered_set<Edge> ball_filter_edge;
            std::unordered_map<VertexID, std::unordered_set<VertexID>> S_w;
            partstart = get_current_time();
            for(auto u : qgraph.GetAllVerticesID()){
	      for (auto v : global_sim[u]){
	        if(whole_ball_nodes[w].find(v) != whole_ball_nodes[w].end()){
		  S_w[u].insert(v);
                  ball_filter_node.insert(v);
                }
              }
            }
            partend = get_current_time();
            parttime[1] += (partend - partstart);

            partstart = get_current_time();
            for(auto e: qgraph.GetAllEdges()){
	      VertexID sourceid=e.src();
              VertexID targetid=e.dst();
              for (auto sim_v1 : S_w[sourceid]){
		for(auto sim_v2 : S_w[targetid]){
                  if (dgraph.ExistEdge(sim_v1,sim_v2)){
//                    ball_filter_edge.insert(Edge(sim_v1,sim_v2,1));
                    ball_filter_edge.insert(Edge(sim_v1, sim_v2));
                  }
                }
              }
            }
            partend = get_current_time();
            parttime[2] += (partend - partstart);
			
            partstart = get_current_time();
	    Ball_View ball_view(ball_filter_node,ball_filter_edge);
            partend = get_current_time();
            parttime[3] += (partend - partstart);

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
            parttime[4] += (partend - partstart);

            partstart = get_current_time();
            Ball_View refined_ball_view(refined_ball_vertex,refinded_ball_edge);
            partend = get_current_time();
            parttime[5] += (partend - partstart);

            partstart = get_current_time();
            rename_sim(refined_ball_view,qgraph,S_w);
            partend = get_current_time();
            parttime[6] += (partend - partstart);

            partstart = get_current_time();
            dual_filter_match(refined_ball_view, qgraph,S_w,w,d_Q);
            partend = get_current_time();
            parttime[7] += (partend - partstart);

            partstart = get_current_time();
            extract_max_pg(refined_ball_view,dgraph,qgraph, w,S_w);
            partend = get_current_time();
            parttime[8] += (partend - partstart);

            partstart = get_current_time();
            max_result.emplace_back(w,S_w);
            partend = get_current_time();
            parttime[9] += (partend - partstart);
            //print_ball_info(qgraph,S_w,w);
			
			// std::fstream out3("calculate_ball_info.txt",std::ios::app);
	        // out3<<"calculate one ball time "<<(float)(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
            // out3.close();
            // std::cout<<"calculate one ball time "<<(float)(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
         
        }
    }
    end = get_current_time();
    time = end - start;
    LOG(INFO) << count;
    LOG(INFO) << "main calculation time: " << time;
    LOG(INFO) << "    --dhop time: " << parttime[0];
    LOG(INFO) << "    --ball filter node time: " << parttime[1];
    LOG(INFO) << "    --ball filter edge time: " << parttime[2];
    LOG(INFO) << "    --ball view time: " << parttime[3];
    LOG(INFO) << "    --refined ball edge time: " << parttime[4];
    LOG(INFO) << "    --refined ball view time: " << parttime[5];
    LOG(INFO) << "    --rename sim time: " << parttime[6];
    LOG(INFO) << "    --dual filter match time: " << parttime[7];
    LOG(INFO) << "    --extract max pg time: " << parttime[8];
    LOG(INFO) << "    --max result time: " << parttime[9];
    LOG(INFO) << "=============================================" ;
//	clock_t etime = clock();
    // std::cout<<"all strong "<< (float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;

	// std::fstream tmp_outfile("whole_time_info_rm_and_and.txt",std::ios::app);
	// tmp_outfile<<"strong sim...calculate d_Q = "<<(float)(e1-s1)/CLOCKS_PER_SEC<<"s"<<std::endl;
	// tmp_outfile<<"strong sim......dual_simulation time = "<<(float)(tmpe-tmps)/CLOCKS_PER_SEC<<"s"<<std::endl;
	// tmp_outfile<<"strong sim.........calculate all ball = "<<(float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
	// tmp_outfile<<"strongr.size = "<<max_result.size()<<std::endl;
	// tmp_outfile.close();
	
    return max_result;
}

void StrongSim::cal_culculate_directed_dhop_nodes(Graph &dgraph, VertexID vid, int d_hop, std::unordered_set<VertexID> &result,std::vector<int> &dis){
    int dgraph_num_vertices = dgraph.GetNumVertices();
    std::vector<int> color(dgraph_num_vertices,0);
    dis.resize(dgraph_num_vertices,INT_MAX - 10);
    std::queue<VertexID> q;
    q.push(vid);
    dis[vid] = 0;
    color[vid] = 1;
    result.insert(vid);
    while(!q.empty()){
         VertexID root = q.front();
         q.pop();
         if(dis[root] == d_hop){
            break ;
         }
         for(auto v : dgraph.GetChildrenID(root)){
             if(color[v] == 0){
                 q.push(v);
                 color[v] = 1;
                 dis[v] = dis[root] + 1;
                 result.insert(v);
             }
         }
         for (auto v: dgraph.GetParentsID(root)){
             if(color[v] == 0){
                 q.push(v);
                 color[v] = 1;
                 dis[v] = dis[root] + 1;
                 result.insert(v);
             }
         }
    }

 }

std::vector<StrongR> StrongSim::strong_simulation_sim_only_add(Graph &dgraph,
                     Graph &qgraph, int flag,
		     std::unordered_map<VertexID,std::unordered_set<VertexID>> &whole_ball_nodes,
		     std::unordered_map<VertexID,std::vector<int>> &whole_dist){
    /**
     *calculate qgaraph diameter
    */
    std::vector<StrongR> max_result;
//	clock_t s1 = clock();
	int d_Q = cal_diameter_qgraph(qgraph);
//	clock_t e1 = clock();
	
	std::unordered_map<VertexID, std::unordered_set<VertexID>> global_sim;
	/**
	*calculate dual simulation for dgraph
	*/
	DualSim dualsim;
	bool inital_sim = false;
//	clock_t tmps = clock();
	dualsim.dual_simulation(dgraph,qgraph,global_sim,inital_sim);
//	clock_t tmpe = clock();
	
	int i=0;
	clock_t stime =clock();
	for (auto w : dgraph.GetAllVerticesID()) {
		/**
		* calculate ball for center w if w if a valid center
		*/
		if (valid_sim_w(qgraph,global_sim,w)){
//			clock_t start =clock();
			/**
			*find d_hop_nodes for w in dgraph
			*/
			
			// std::unordered_set<VertexID> ball_node;			
			// dgraph.find_hop_nodes(w,d_Q,ball_node); // replaced by the next line.
			cal_culculate_directed_dhop_nodes(dgraph, w, d_Q, whole_ball_nodes[w], whole_dist[w]);			
			
			std::unordered_set<VertexID> ball_filter_node;
			std::unordered_set<Edge> ball_filter_edge;
            std::unordered_map<VertexID, std::unordered_set<VertexID>> S_w;
            for(auto u : qgraph.GetAllVerticesID()){
				for (auto v : global_sim[u]){
					// if(ball_node.find(v) != ball_node.end()){
					if(whole_ball_nodes[w].find(v) != whole_ball_nodes[w].end()){
						S_w[u].insert(v);
                        ball_filter_node.insert(v);
                    }
                }
            }
            for(auto e: qgraph.GetAllEdges()){
				VertexID sourceid=e.src();
                VertexID targetid=e.dst();
                for (auto sim_v1 : S_w[sourceid]){
					for(auto sim_v2 : S_w[targetid]){
                        if (dgraph.ExistEdge(sim_v1,sim_v2)){
//                            ball_filter_edge.insert(Edge(sim_v1,sim_v2,1));
                          ball_filter_edge.insert(Edge(sim_v1, sim_v2));
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
            //print_ball_info(qgraph,S_w,w);
//            clock_t end = clock();
			
//            std::fstream out3("calculate_ball_info.txt",std::ios::app);
//	        out3<<"calculate one ball time "<<(float)(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
//            out3.close();
            // std::cout<<"calculate one ball time "<<(float)(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
        }
    }
//	clock_t etime=clock();
    // std::cout<<"all strong "<< (float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
	
//	std::fstream tmp_outfile("whole_time_info_only_add.txt",std::ios::app);
//	tmp_outfile<<"strong sim...calculate d_Q = "<<(float)(e1-s1)/CLOCKS_PER_SEC<<"s"<<std::endl;
//	tmp_outfile<<"strong sim......dual_simulation time = "<<(float)(tmpe-tmps)/CLOCKS_PER_SEC<<"s"<<std::endl;
//	tmp_outfile<<"strong sim.........calculate all ball = "<<(float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
//	tmp_outfile<<"strongr.size = "<<max_result.size()<<std::endl;
//	tmp_outfile.close();

    return max_result;
}


std::vector<StrongR> StrongSim::strong_simulation_sim(Graph &dgraph, Graph &qgraph, int flag){
    /**
     *calculate qgaraph diameter
    */
    std::vector<StrongR> max_result;
	clock_t s1 = clock();
	int d_Q = cal_diameter_qgraph(qgraph);
	clock_t e1 = clock();
	
	std::unordered_map<VertexID, std::unordered_set<VertexID>> global_sim;
	/**
	*calculate dual simulation for dgraph
	*/
	DualSim dualsim;
	bool inital_sim = false;
	clock_t tmps = clock();
	dualsim.dual_simulation(dgraph,qgraph,global_sim,inital_sim);
	clock_t tmpe = clock();
	
    // std::cout<<"dual "<<(float)(end1-start1)/CLOCKS_PER_SEC<<std::endl;
	int i=0;
	clock_t stime =clock();
	for (auto w : dgraph.GetAllVerticesID()) {
		/**
		* calculate ball for center w if w if a valid center
		*/
		if (valid_sim_w(qgraph,global_sim,w)){
			clock_t start =clock();
			/**
			*find d_hop_nodes for w in dgraph
			*/
			std::unordered_set<VertexID> ball_node;
			dgraph.find_hop_nodes(w,d_Q,ball_node);
			
			std::unordered_set<VertexID> ball_filter_node;
			std::unordered_set<Edge> ball_filter_edge;
            std::unordered_map<VertexID, std::unordered_set<VertexID>> S_w;
            for(auto u : qgraph.GetAllVerticesID()){
				for (auto v : global_sim[u]){
					if(ball_node.find(v) != ball_node.end()){
						S_w[u].insert(v);
                        ball_filter_node.insert(v);
                    }
                }
            }
            for(auto e: qgraph.GetAllEdges()){
				VertexID sourceid=e.src();
                VertexID targetid=e.dst();
                for (auto sim_v1 : S_w[sourceid]){
					for(auto sim_v2 : S_w[targetid]){
                        if (dgraph.ExistEdge(sim_v1,sim_v2)){
//                            ball_filter_edge.insert(Edge(sim_v1,sim_v2,1));
                          ball_filter_edge.insert(Edge(sim_v1, sim_v2));
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
            //print_ball_info(qgraph,S_w,w);
            clock_t end = clock();
			
            std::fstream out3("calculate_ball_info.txt",std::ios::app);
	        out3<<"calculate one ball time "<<(float)(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
            out3.close();
            // std::cout<<"calculate one ball time "<<(float)(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
        }
    }
	clock_t etime=clock();
    // std::cout<<"all strong "<< (float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
	
	if(flag == 0){
		std::fstream tmp_outfile("whole_time_info_only_rm.txt",std::ios::app);
		tmp_outfile<<"strong sim...calculate d_Q = "<<(float)(e1-s1)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"strong sim......dual_simulation time = "<<(float)(tmpe-tmps)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"strong sim.........calculate all ball = "<<(float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
		tmp_outfile<<"strongr.size = "<<max_result.size()<<std::endl;
		tmp_outfile.close();
	}
	else if(flag == 1){
		std::fstream tmp_outfile("whole_time_info_only_add.txt",std::ios::app);
		tmp_outfile<<"strong sim...calculate d_Q = "<<(float)(e1-s1)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"strong sim......dual_simulation time = "<<(float)(tmpe-tmps)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"strong sim.........calculate all ball = "<<(float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
		tmp_outfile<<"strongr.size = "<<max_result.size()<<std::endl;
		tmp_outfile.close();
	}
	else if(flag == 2){
		std::fstream tmp_outfile("whole_time_info_rm_and_add.txt",std::ios::app);
		tmp_outfile<<"strong sim...calculate d_Q = "<<(float)(e1-s1)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"strong sim......dual_simulation time = "<<(float)(tmpe-tmps)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"strong sim.........calculate all ball = "<<(float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
		tmp_outfile<<"strongr.size = "<<max_result.size()<<std::endl;
		tmp_outfile.close();
	}
	else if(flag == 3){
		std::fstream tmp_outfile("whole_time_info_add_and_rm.txt",std::ios::app);
		tmp_outfile<<"strong sim...calculate d_Q = "<<(float)(e1-s1)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"strong sim......dual_simulation time = "<<(float)(tmpe-tmps)/CLOCKS_PER_SEC<<"s"<<std::endl;
		tmp_outfile<<"strong sim.........calculate all ball = "<<(float)(etime-stime)/CLOCKS_PER_SEC<<std::endl;
		tmp_outfile<<"strongr.size = "<<max_result.size()<<std::endl;
		tmp_outfile.close();
	}
    return max_result;
}


