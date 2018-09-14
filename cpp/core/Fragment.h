#ifndef CPP_FRAGMENT_H_
#define CPP_FRAGMENT_H_
#include "cpp/core/global.h"
#include "cpp/core/vertex.h"
#include "cpp/core/edge.h"
#include "cpp/core/graph.h"
#include "cpp/core/config.h"
#include "cpp/core/FragmentLoader.h"
#include "cpp/utils/MessageBuffer.h"
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <string>
#include <vector>
#include <utility>
//#include <map>
/**
  * a class record each fragment graph and processor communication information
  */
class Fragment{
public :
    Fragment();

    Fragment(Graph &graph, const std::string v_file, const std::string e_file, const std::string r_file);

    Fragment(Graph &graph,const std::vector<Vertex> &_vertices, const std::vector<Edge> &_edges,
              const std::unordered_map<VertexID,int> &_fragTable, int _FID);

    ~Fragment();

    int getVertexFragmentID(VertexID u) const;

    int getLocalID(VertexID gvid);

    int getGlobalID(VertexID localid);

    int getNumVertices() const;

    int getNumEdges() const;

    const std::unordered_set<VertexID>* getInnerVertices() const;

    const std::unordered_set<VertexID>* getOuterVertices() const;

    const bool isBorderVertex(const VertexID gvid);

    const bool isInnerVertex(const VertexID gvid);

    const bool isOuterVertex(const VertexID gvid);

	const std::bitset<NUM_FRAGMENTS> getMsgThroughDest(const VertexID gvid);

	const bool has_vertex(const VertexID gvid) const;

	void add_global_info(const VertexID gvid, const VertexID localid);

	void add_local_info(const VertexID localid, const VertexID gvid);

	void add_outerVertices(const VertexID gvid);

    void update_by_add_vertices(Graph &dgraph,
         std::vector<std::pair<VertexID, VertexLabel>>& add_vertices) {
      int fnum = get_num_workers();
      for (auto &v : add_vertices) {
        auto gid = v.first;
        if (gid % fnum == FID) {
          VertexID lid = dgraph.AddVertex(Vertex(gid, v.second));
          local2global[lid] = gid;
          global2local[gid] = lid;
          innerVertices.insert(gid);
          fragTable[gid] = FID;
          numVertices++;
        }
      }
    }

    void update_by_add_edges(Graph &dgraph,
            std::unordered_set<std::pair<VertexID, VertexID>> &add_edges,
            bool communication_next) {
      std::unordered_set<Vertex> vertices_;
      std::unordered_set<Edge> edges_;
      for (auto &e : add_edges) {
        edges_.insert(Edge(e.first, e.second));
        if (this->isInnerVertex(e.first) && !this->isInnerVertex(e.second)) {
          const int dst = this->getVertexFragmentID(e.second); 
          vertexBuffers.add_message(dst, Vertex(e.first, dgraph.GetVertexLabel(this->getLocalID(e.first)))); 
        } else if (this->isInnerVertex(e.second) && !this->isInnerVertex(e.first)) {
          const int dst = this->getVertexFragmentID(e.first);
          vertexBuffers.add_message(dst, Vertex(e.second, dgraph.GetVertexLabel(this->getLocalID(e.second))));
        } 
      }
      vertexBuffers.sync_messages();

      for (auto &iter : vertexBuffers.get_messages()) {
        if (!this->isOuterVertex(iter.id())) {
          vertices_.insert(iter);
        }
      }
      vertexBuffers.reset_in_messages();
      this->update_fragment_add_edges(dgraph, edges_, vertices_, communication_next);
    }

    template<class T1,class T2>
    void update_fragment_add_edges(Graph &graph,T1 &add_edges,T2 &vertices,bool communication_next=true){
        for(auto ver : vertices){
            VertexID gvid = ver.id();
            if(global2local.find(gvid) == global2local.end()){
                VertexID localid = graph.AddVertex(ver);
                local2global[localid] = gvid;
                global2local[gvid] = localid;
                outerVertices.insert(gvid);
                numVertices++;
            }
        }
        LOG(INFO) << "inner 4.1 ";
        for(auto edge :add_edges){
            VertexID src = edge.src();
            VertexID dst = edge.dst();
            if(communication_next){
                if (fragTable.at(src) == FID || fragTable.at(dst) == FID){
//                    graph.AddEdge(Edge(global2local[src], global2local[dst], edge.attr()));
                     graph.AddEdge(Edge(global2local[src], global2local[dst]));
                     numEdges++;
                     if(fragTable.at(src) != FID && fragTable.at(dst) ==FID){
                        msgThroughDest[dst].set(fragTable.at(src));
                    }else if(fragTable.at(src) == FID && fragTable.at(dst) != FID){
                        msgThroughDest[src].set(fragTable.at(dst));
                    }
                }
            }else{
                if(global2local.find(src)!=global2local.end() && global2local.find(dst) != global2local.end()){
//                    graph.AddEdge(Edge(global2local[src], global2local[dst], edge.attr()));
                  numEdges++;
                  graph.AddEdge(Edge(global2local[src], global2local[dst]));
                }
            }
        }
        LOG(INFO) << "inner 4.2";
//            if(global2local.find(src)!=global2local.end() && global2local.find(dst) != global2local.end()){
//                graph.AddEdge(Edge(global2local[src], global2local[dst], edge.attr()));
//                if(communication_next){
//                    if(fragTable.at(src) != FID && fragTable.at(dst) ==FID){
//                        msgThroughDest[dst].set(fragTable.at(src));
//                    }else if(fragTable.at(src) == FID && fragTable.at(dst) != FID){
//                        msgThroughDest[src].set(fragTable.at(dst));
//                    }
//                }
//            }

       graph.RebuildGraphProperties();
       LOG(INFO) << "inner 4.3";
    }

    void update_by_remove_edges(Graph &graph, std::unordered_set<std::pair<VertexID, VertexID>> &rm_edges, bool communication_next) {
      std::vector<Edge> rm_edges_;
      for (auto &e : rm_edges) {
        rm_edges_.emplace_back(e.first, e.second);
      }
      this->update_fragment_remove_edges(graph, rm_edges_, communication_next);
    }

    template<class T1>
    void update_fragment_remove_edges(Graph &graph,T1 &rm_edges,bool communication_next=true){
        for(auto edge: rm_edges){
            VertexID src = edge.src();
            VertexID dst = edge.dst();
            if(global2local.find(src)!=global2local.end() && global2local.find(dst) != global2local.end()){
//                graph.RemoveEdge(Edge(global2local[edge.src()], global2local[edge.dst()], edge.attr()));
              graph.RemoveEdge(Edge(global2local[src], global2local[dst]));
              numEdges--;
            }
        }
        if(communication_next){
            msgThroughDest.clear();
            for(auto v:innerVertices){
                VertexID local_id = getLocalID(v);
                for(auto u_pre : graph.GetParentsID(local_id)){
                    VertexID global_id = getGlobalID(u_pre);
                    msgThroughDest[v].set(fragTable.at(global_id));
                }
                for(auto u_des :graph.GetChildrenID(local_id)){
                    VertexID global_id = getGlobalID(u_des);
                    msgThroughDest[v].set(fragTable.at(global_id));
                }
            }
        }
        graph.RebuildGraphProperties();
    }
private :
       int numVertices;
       int numEdges;
       std::unordered_map<VertexID,int> fragTable;
       std::unordered_map<VertexID,VertexID> global2local;
       std::unordered_map<VertexID,VertexID> local2global;
       std::unordered_set<VertexID> innerVertices;
       std::unordered_set<VertexID> outerVertices;
       int FID;
       std::unordered_map<VertexID, std::bitset<NUM_FRAGMENTS>> msgThroughDest;//recored meassage should send to process's id
       MessageBuffer<Vertex>  vertexBuffers;
};
#endif //CPP_FRAGMENT_H_
