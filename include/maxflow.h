//
// Created by wangxuedong on 11/6/23.
//

#ifndef SEQGRAPH_MAXFLOW_H
#define SEQGRAPH_MAXFLOW_H
#include "Graph.h"
#include "algorithm"
#include <deque>
#include "util.h"
#include <cstring>
#include <algorithm>
#include <queue>
#include <iostream>
#include "limits"
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <cfloat>
#include <climits>
extern int VERBOSE;
extern bool BREAK_C;
extern bool IS_RNA;
extern const float ZERO;
extern const float INF;
extern const double epsilon;
extern const float edge_epsilon;


class maxFlow {
private:
    seqGraph::Graph* graph;
    uint V;
    std::vector<std::vector<int>> augmentationPaths;
    float ** ConjugateMatrix;
    int sId;
    int tId;
    std::map<uint, std::string> matrixIdx2nodeName;
    std::map<uint, std::string> matrixIdx2nodeIdxStr;
    std::map<int, int> matrixIdx2nodeIdx;
    int startNode;
    int graphInNode;
    int graphOutNode;



public:
    seqGraph::Graph *getGraph() const;

private:
    float** currentMatrix;
    float** originalMatrix;
    std::vector<seqGraph::Vertex*>* originalVertices;
    seqGraph::Graph* originalGraph;
    int N;


public:
    explicit maxFlow(seqGraph::Graph* graph1);
    ~maxFlow();
    std::vector<std::vector<seqGraph::preTransFrag>> preTrans;
    std::vector<std::vector<int>> preTransPath;

    void printM(int i);

    static bool cmpVertex (int i,int j);
    inline int getN() const {
        return N;
    }

    inline float** getMatrix() const {
        return currentMatrix;
    };
//    inline seqGraph::SparseMatrix& getMatrix() const {
//        return this->graph->getConjugateMatrix();
//    }

    std::map<int, std::vector<int>*>* resolvePath(std::map<int, std::vector<int>*>* prevPaths);





    bool vertexLookup(int i, int j);
    int  inDegree(int idx);
    int  outDegree(int idx);
    bool bfs(float **rGraph, int s, int t, int parent[]);
    void dfs(float **rGraph, int s, int t, bool visited[], std::vector<int>& path);
    void fordFulkerson();



    void iterate();

    void calculateMaxFlow();

    void subtractMaxPath(const std::pair<std::vector<int>, float> &max_flow_path_info);

    std::pair<std::vector<int>, float> findMaxPath();

    std::vector<int> findMaxWeightPathForward(float **matrix, int start);

    std::vector<int> findMaxWeightPathBackward(float **matrix, int start);

    std::vector<int> findMaxWeightPathBackward(int start, int end);

    std::vector<std::pair<int, int>> findPredecessors(int v);

    std::vector<int> findMaxWeightPathForward(int start, int end);

    std::vector<std::pair<int, int>> findSuccessors(int v);


    bool is_subsequence(const std::vector<int> &input_path, const std::vector<int> &file_path);

    std::vector<int> rabin_karp(const std::vector<int> &text, const std::vector<int> &pattern);

    bool check_equal(const std::vector<int> &seq1, int start1, const std::vector<int> &seq2, int start2, int length);

    int compute_hash(const std::vector<int> &seq);

    std::vector<int> kmpProcessPattern(const std::vector<int> &pattern);

    std::vector<int> kmpSearch(const std::vector<int> &text, const std::vector<int> &pattern);


    float pushMaxFlow(std::vector<int> &path, std::vector<float> &nodeFlux);

    std::map<uint, std::string> GetmatrixIdx2nodeIdxStr();
};


#endif SEQGRAPH_MAXFLOW_BK_HDR

