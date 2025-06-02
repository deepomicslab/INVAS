//
// Created by wangxuedong on 11/6/23.
//

//
// Created by caronkey on 28/12/2021.
//

#include "../include/maxflow.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

// 使用较大的质数作为基数
const int base = 257;
// 选择一个足够大的数以避免哈希冲突
const int mod = 1e9 + 9;




maxFlow::maxFlow(seqGraph::Graph* graph1) {
    graph = graph1;
    graph->initMatrix();
    graph->updateMaxI();
    ConjugateMatrix = graph->getConjugateMatrix2();

    /// print matrix
    int N = 4 * graph1->getVCount()-4;
    for (int i=0; i< N; i++) {
        if (i == N-4) {
            matrixIdx2nodeName[i] = "s5+";
            matrixIdx2nodeIdxStr[i] = "s5+";
            std::cout<<"s\t";
        } else if (i == N-3) {
            matrixIdx2nodeName[i] = "s3+";
            matrixIdx2nodeIdxStr[i] = "s3+";
            std::cout<<"t\t";
        } else if (i == N-2) {
            matrixIdx2nodeName[i] = "t5+";
            matrixIdx2nodeIdxStr[i] = "t5+";
            std::cout<<"s'\t";
        } else if (i == N-1) {
            matrixIdx2nodeName[i] = "t3+";
            matrixIdx2nodeIdxStr[i] = "t3+";
            std::cout<<"t'\t";
        } else {
            int node_idx = i/4;
            int node_type = i%4;
            std::string node_name = (*graph->getVertices())[node_idx]->getOriginId();
            std::string node_type_str;
            if (node_type == 0) {
                node_type_str = " 5+";
            } else if (node_type == 1) {
                node_type_str = " 3+";
            } else if (node_type == 2) {
                node_type_str = " 5-";
            } else if (node_type == 3) {
                node_type_str = " 3-";
            }
            matrixIdx2nodeName[i] = node_name+node_type_str;
            matrixIdx2nodeIdxStr[i] = (std::to_string(node_idx)+node_type_str);
            matrixIdx2nodeIdx[i] = node_idx;
            std::cout<<node_name<<node_type_str<<"\t";
        }
    }
    V=N;
    sId=N-3;
    tId=N-2;
    std::cout<<" \t";
    for (int i = 1; i < N + 1; i++) {
        std::cout<<i<<"\t";
    }
    std::cout<<std::endl;
    std::cout<< "original matrix: "<<std::endl;
    for (int i=0; i< N; i++) {
        std::cout<<"[";
        for (int j=0; j< N; j++) {
            std::cout<<graph1->getConjugateMatrix2()[i][j]<<",";
        }
        std::cout<<"],"<<std::endl;
    }
}

maxFlow::~maxFlow() {
    if (graph != nullptr)
        free(graph);
}


bool maxFlow::bfs(float **rGraph, int s, int t, int parent[])
{
    bool visited[V];
    memset(visited, 0, sizeof(visited));

    std::queue<int> q;
    q.push(s);
    visited[s] = true;
    parent[s] = -1;

//    std::cout << "Starting BFS from source: " << s << std::endl;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

//        std::cout << "Exploring node: " << u << std::endl;

        for (int v = 0; v < V; v++) {
            if (!visited[v] && rGraph[u][v] > 0) {
//                std::cout << "  Found path to node: " << v << ", capacity: " << rGraph[u][v] << std::endl;
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }

    if (visited[t]) {
        std::vector<int> path;
        for (int v = t; v != s; v = parent[v]) {
            path.push_back(v);
        }
        path.push_back(s);
        reverse(path.begin(), path.end());
        augmentationPaths.push_back(path);

    }

    return visited[t];
}
void maxFlow::calculateMaxFlow()
{
    int u, v;
    int s = sId;
    int t = tId;
    float** rGraph = new float*[V];
    for(int i = 0; i < V; ++i) {
        rGraph[i] = new float[V];
        for(int j = 0; j < V; ++j) {
            rGraph[i][j] = ConjugateMatrix[i][j];
        }
    }

    int parent[V];
    float max_flow = 0; // Initialize max flow

    while (bfs(rGraph, s, t, parent)) {
        float path_flow = FLT_MAX;

        for (v = t; v != s; v = parent[v]) {
            u = parent[v];
//            std::cout<< path_flow << " " << rGraph[u][v] << std::endl;
            path_flow = std::min(path_flow, rGraph[u][v]);
        }

        for (v = t; v != s; v = parent[v]) {
            u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }
        max_flow += path_flow;
    }

    for(int i = 0; i < V; ++i)
        delete [] rGraph[i];
    delete [] rGraph;
}

std::pair<std::vector<int>, float> maxFlow::findMaxPath()
{
    std::set<std::vector<int>> uniquePaths(augmentationPaths.begin(), augmentationPaths.end());
    std::vector<int> max_flow_path;
    float max_min_path_weight = -INF; // 假设-INF为负无穷大

    // 在所有路径中找到最大最小权重的路径
    for (const auto& path : uniquePaths) {
        float min_weight_in_path = INF; // 假设INF为正无穷大
        for (int i = 0; i < path.size() - 1; ++i) {
            min_weight_in_path = std::min(min_weight_in_path, ConjugateMatrix[path[i]][path[i + 1]]);
        }
        // 找到最大的最小权重值
        if (min_weight_in_path > max_min_path_weight) {
            max_flow_path = path;
            max_min_path_weight = min_weight_in_path;
        }
    }

    return {max_flow_path, max_min_path_weight};
}

void maxFlow::subtractMaxPath(const std::pair<std::vector<int>, float>& max_flow_path_info)
{
    std::vector<int> max_flow_path = max_flow_path_info.first;
    float min_weight = max_flow_path_info.second;

    for (int i = 1; i < max_flow_path.size() - 2; ++i) {
        ConjugateMatrix[max_flow_path[i]][max_flow_path[i + 1]] -= min_weight;
        /// minus conjugate edge flow as well
        int con_raw, con_col;
        con_raw = graph->getConjugateIndexMap()[std::make_pair(max_flow_path[i], max_flow_path[i + 1])].first;
        con_col = graph->getConjugateIndexMap()[std::make_pair(max_flow_path[i], max_flow_path[i + 1])].second;
        ConjugateMatrix[con_raw][con_col] -= min_weight;
//        std::cout << "subtracting " << min_weight << " from " << max_flow_path[i] << " to " << max_flow_path[i + 1] << std::endl;
//        std::cout << "subtracting conjugate " << min_weight << " from " << con_raw << " to " << con_col << std::endl;
        if(ConjugateMatrix[max_flow_path[i]][max_flow_path[i + 1]] < ZERO) { // Avoid negative weights due to floating point precision
            ConjugateMatrix[max_flow_path[i]][max_flow_path[i + 1]] = 0;
        }
        if(ConjugateMatrix[con_raw][con_col] < ZERO) { // Avoid negative weights due to floating point precision
            ConjugateMatrix[con_raw][con_col] = 0;
        }
    }
}

std::vector<int> maxFlow::findMaxWeightPathForward(float** matrix, int start) {
    std::vector<int> path;
    int current = start;

    while (current < V - 1) {
        float maxWeight = 0;
        int maxNode = -1;

        for (int i = current + 1; i < V; ++i) {
            std::cout << "考虑节点 " << current << " 到节点 " << i << " 的权重: " << matrix[current][i] << std::endl;
            if (matrix[current][i] > maxWeight) {
                maxWeight = matrix[current][i];
                maxNode = i;
            }
        }

        if (maxNode == -1) { // 没有更多的节点可以前进
            std::cout << "节点 " << current << " 没有更多的前进路径。\n";
            break;
        }

        path.push_back(maxNode);
        std::cout << "选择节点 " << maxNode << " 以最大权重前进。\n";
        current = maxNode;
    }

    return path;

}

std::vector<int> maxFlow::findMaxWeightPathBackward(float** matrix, int start) {
    std::vector<int> path;
    int current = start;

    while (current > 0) {
        float maxWeight = 0;
        int maxNode = -1;

        for (int i = current - 1; i >= 0; --i) {
            std::cout << "考虑节点 " << i << " 到节点 " << current << " 的权重: " << matrix[i][current] << std::endl;
            if (matrix[i][current] > maxWeight) {
                maxWeight = matrix[i][current];
                maxNode = i;
            }
        }

        if (maxNode == -1) { // 没有更多的节点可以后退
            std::cout << "节点 " << current << " 没有更多的后退路径。\n";
            break;
        }

        path.push_back(maxNode);
        std::cout << "选择节点 " << maxNode << " 以最大权重后退。\n";
        current = maxNode;
    }

    // 得到的路径是从start到0的顺序，需要反转来表示从0到start的路径
    std::reverse(path.begin(), path.end());

    return path;
}

std::vector<std::pair<int, int>> maxFlow::findPredecessors(int v) {
    std::vector<std::pair<int, int>> predecessors;
    for (int i = 0; i < V; ++i) {
        if (ConjugateMatrix[i][v] != 0) { // 如果i到v有边，则i是v的前序节点
            predecessors.push_back({i, ConjugateMatrix[i][v]});
        }
    }
    return predecessors;
}

std::vector<std::pair<int, int>> maxFlow::findSuccessors(int v) {
    std::vector<std::pair<int, int>> successors;
    // 遍历邻接矩阵的第v行来找到所有从节点v出发的边
    for (int j = 0; j < V; ++j) {
        if (ConjugateMatrix[v][j] != 0) { // 如果v到j有边，则j是v的后序节点
            successors.push_back({j, ConjugateMatrix[v][j]});
        }
    }
    return successors;
}

std::vector<int> maxFlow::findMaxWeightPathBackward(int start, int end) {
    std::vector<int> path;
    if (start < 0 || start >= V || end < 0 || end >= V) return path; // 无效的顶点检查

    int current = end;
    path.push_back(end);

    // 从终点向起点遍历，每次选择最大权重的前序节点
    while (current != start) {
        std::vector<std::pair<int, int>> preds = findPredecessors(current);
        if (preds.empty()) { // 如果没有前序节点，则无法到达起点
            return std::vector<int>(); // 返回空路径
        }

        // 找到权重最大的前序节点
        int maxWeight = INT_MIN;
        int maxPred = -1;
        for (auto &pred : preds) {
            int u = pred.first;
            int weight = pred.second;
            if (weight > maxWeight) {
                maxWeight = weight;
                maxPred = u;
            }
        }

        // 如果找到了权重更大的前序节点，则更新当前节点并将其添加到路径中
        if (maxPred != -1 && maxWeight != INT_MIN) {
            current = maxPred;
            path.push_back(current);
        } else {
            // 如果没有找到权重更大的前序节点，则说明没有可行路径
            return std::vector<int>(); // 返回空路径
        }
    }

    std::reverse(path.begin(), path.end()); // 反转路径，因为我们是从终点回溯到起点
    return path;
}

std::vector<int> maxFlow::findMaxWeightPathForward(int start, int end) {
    std::vector<int> path;
    if (start < 0 || start >= V || end < 0 || end >= V) return path; // 无效的顶点检查

    // 从起点开始寻找最大权重的后继节点
    int current = start;
    while (current != end) {
        path.push_back(current); // 将当前节点添加到路径中

        std::vector<std::pair<int, int>> successors = findSuccessors(current);
        if (successors.empty()) { // 如果没有后继节点，则无法到达终点
            return std::vector<int>(); // 返回空路径
        }

        // 找到权重最大的后继节点
        int maxWeight = INT_MIN;
        int maxSucc = -1;
        for (const auto &succ : successors) {
            int v = succ.first;
            int weight = succ.second;
            if (weight > maxWeight) {
                maxWeight = weight;
                maxSucc = v;
            }
        }

        // 如果找到了权重更大的后继节点，则更新当前节点
        if (maxSucc != -1 && maxWeight != INT_MIN) {
            current = maxSucc;
        } else {
            // 如果没有找到权重更大的后继节点，则说明没有可行路径
            return std::vector<int>(); // 返回空路径
        }
    }

    path.push_back(end); // 确保终点也在路径中
    return path; // 返回包含起点到终点的路径
}


void maxFlow::iterate()
{
    int i = 0;
//    while (true) {
//        augmentationPaths.clear();
//        calculateMaxFlow();
//        if (augmentationPaths.empty()) {
//            break;
//        }
//        std::pair<std::vector<int>, float> res = findMaxPath();
//        std::cout << "\nAugmentation paths:\n";
//        ///
//        std::vector<int> maxPath = res.first;
//        for (int node : maxPath) {
//            std::cout << node <<": (" <<matrixIdx2nodeIdxStr[node] <<")"<< "-->";
//        }
//        std::cout << "min flow: " << res.second << "\n";
//
//        // Print the matrix after each iteration
//        std::cout << "Max flow matrix after this iteration: " << std::endl;
//        for (int i = 0; i < V; i++) {
//            for (int j = 0; j < V; j++) {
//                std::cout << ConjugateMatrix[i][j] << " ";
//            }
//            std::cout << std::endl;
//        }
//        i++;
//        subtractMaxPath(res);
//        if (i==3){
//            break;
//        }
//    }
//    std::cout<< "final matrix: "<<std::endl;
//    for (int i=0; i< V; i++) {
//        std::cout<<"[";
//        for (int j=0; j< V; j++) {
//            std::cout<<ConjugateMatrix[i][j]<<",";
//        }
//        std::cout<<"],"<<std::endl;
//    }

    fprintf(stderr, "current matrix: \n");
    for (int i=0; i< V; i++) {
        fprintf(stderr, "[");
        for (int j=0; j< V; j++) {
            fprintf(stderr, "%f,", ConjugateMatrix[i][j]);
        }
        fprintf(stderr, "],\n");
    }
    startNode = 4*(graph->maxI+1)-3;
    graphInNode = 4*(graph->getVCount()-2)+1;
    graphOutNode = 4*(graph->getVCount()-2)+2;
    std::vector<int> backwardPath = findMaxWeightPathBackward(graphInNode, startNode);
    if (backwardPath.empty()) {
        std::cout << "No backward path found from " << graphInNode << " to " << startNode << std::endl;
        return;
    }
    std::cout << "Max weight backward path to " << startNode << ": ";
    for (int node : backwardPath) {
        std::cout << node << " ";
    }
    std::cout<<std::endl;

    // 从给定节点向前找路径
    std::vector<int> forwardPath = findMaxWeightPathForward(startNode, graphOutNode);
    if (forwardPath.empty()) {
        std::cout << "No forward path found from " << startNode << " to " << graphOutNode << std::endl;
        return;
    }
    std::cout << "Max weight forward path from " << startNode << ": ";
    for (int node : forwardPath) {
        std::cout << node << " ";
    }
    std::cout << std::endl;

    backwardPath.insert(backwardPath.end(), forwardPath.begin()+1, forwardPath.end());
    std::cout<<"Find path: ";
    for (int node : backwardPath) {
        std::cout << node << " ";
    }
    fprintf(stderr, "Find path: ");
    for (int i=0; i< backwardPath.size(); i++) {
        fprintf(stderr, "%d,", backwardPath[i]);
    }
    fprintf(stderr, "\n");

    std::cout << std::endl;
    if (backwardPath.front()!=graphInNode || backwardPath.back()!=graphOutNode){
        fprintf(stderr, "return because path not found\n");
        fprintf(stderr, "left node covs:\n");
        for (auto &v: *graph->getVertices()){
            fprintf(stderr, "%f\n", v->getWeight()->getCoverage());
        }
        fprintf(stderr, "final matrix:\n");
        for (int i=0; i< V; i++) {
            fprintf(stderr, "[");
            for (int j=0; j< V; j++) {
                fprintf(stderr, "%f,", ConjugateMatrix[i][j]);
            }
            fprintf(stderr, "],\n");
        }
        fprintf(stderr, "final transfrags:\n");
        for (int trans_i=0; trans_i < graph->getTransFragment().size(); trans_i++){
            fprintf(stderr, "%d\t%f\n", trans_i, graph->getTransFragment()[trans_i].aboundance);

        }
        return;
    }

    std::vector<float> nodeFlux;
    float flux = pushMaxFlow(backwardPath, nodeFlux);
    if (flux<=epsilon)
    {
        fprintf(stderr, "return because %f<=%F\n", flux, epsilon);
        fprintf(stderr, "left node covs:\n");
        for (auto &v: *graph->getVertices()){
            fprintf(stderr, "%f\n", v->getWeight()->getCoverage());
        }
        /// print final transfrg

        /// print conjugatematrix
        fprintf(stderr, "final matrix:\n");

        for (int i=0; i< V; i++) {
            fprintf(stderr, "[");
            for (int j=0; j< V; j++) {
                fprintf(stderr, "%f,", ConjugateMatrix[i][j]);
            }
            fprintf(stderr, "],\n");
        }
        fprintf(stderr, "final transfrags:\n");
        for (int trans_i=0; trans_i < graph->getTransFragment().size(); trans_i++){
            fprintf(stderr, "%d\t%f\n", trans_i, graph->getTransFragment()[trans_i].aboundance);

        }
        return;
    }


//    if(flux>epsilon) {

        // store trans path, calculate used covs
    int pre_node_id = -1;
    float t_len=0.0;
    float t_cov_total = 0.0;
    std::vector<seqGraph::preTransFrag> transPath;
    for(int path_i=0; path_i <backwardPath.size(); path_i++){
        auto &path = backwardPath[path_i];
        int node_idx = matrixIdx2nodeIdx[path];

        if (path==backwardPath.front() || path == backwardPath.back() || node_idx == pre_node_id){
            continue;
        }
        auto &vertex = (*graph->getVertices())[node_idx];
        float exon_cov;
        float used_cov;

        if ((path==node_idx*4+2||path==node_idx*4+3)&&vertex->getSignal()=="INV"){
            exon_cov = vertex->getWeight()->getConCoverage();
            used_cov = exon_cov * nodeFlux[path_i];
            vertex->getWeight()->setConCoverage(exon_cov - used_cov);
        }
        else{
            exon_cov = vertex->getWeight()->getCoverage();
            used_cov = exon_cov * nodeFlux[path_i];
            vertex->getWeight()->setCoverage(exon_cov - used_cov);
        }

        auto &vertex_Id_str=vertex->getId();
        std::stringstream id_str(vertex_Id_str);
        std::string segment;

        // 获取 start
        std::getline(id_str, segment, '-');
        int start = std::stoi(segment);

        // 获取 end
        std::getline(id_str, segment, '_');
        int end = std::stoi(segment);

        // 获取 copy
        std::getline(id_str, segment);
        int copy = std::stoi(segment);
        seqGraph::preTransFrag pre_trans_f(node_idx, used_cov, start, end);
        transPath.emplace_back(pre_trans_f);
        pre_node_id = node_idx;
        t_len+= (end-start);
        t_cov_total+=((end-start)*used_cov);
    }
    float t_cov = t_cov_total/t_len;
//    if (t_cov <1){
//        fprintf(stderr, "return because t_cov: %f <1\n", t_cov);
//        fprintf(stderr, "left node covs:\n");
//        for (auto &v: *graph->getVertices()){
//            fprintf(stderr, "%f\n", v->getWeight()->getCoverage());
//        }
//        /// print conjugatematrix
//        std::cout<< "final matrix: "<<std::endl;
//        for (int i=0; i< V; i++) {
//            fprintf(stderr, "[");
//            for (int j=0; j< V; j++) {
//                fprintf(stderr, "%f,", ConjugateMatrix[i][j]);
//            }
//            fprintf(stderr, "],\n");
//        }
//        fprintf(stderr, "final transfrags:\n");
//        for (int trans_i=0; trans_i < graph->getTransFragment().size(); trans_i++){
//            fprintf(stderr, "%d\t%f\n", trans_i, graph->getTransFragment()[trans_i].aboundance);
//
//        }
//        return;
//
//    }
    preTrans.emplace_back(transPath);
    preTransPath.emplace_back(backwardPath);
    graph->updateMaxI();
    for (int trans_i=0; trans_i < graph->getTransFragment().size(); trans_i++){
        fprintf(stderr, "%d\t%f\n", trans_i, graph->getTransFragment()[trans_i].aboundance);

    }



//    }
//    for (auto &pretrans: preTrans) {
//        std::cout<<"pretrans: "<<std::endl;
//        for (auto &pretransfrag: pretrans) {
//            std::cout<<pretransfrag.seg_idx<<" "<<pretransfrag.cov<<" "<<pretransfrag.start<<" "<<pretransfrag.end<<std::endl;
//        }
//    }
    iterate();





}

// KMP算法的预处理函数，生成部分匹配表
std::vector<int> maxFlow::kmpProcessPattern(const std::vector<int>& pattern) {
    int m = pattern.size();
    std::vector<int> lps(m, 0);
    int len = 0; // 表示lps[0...len-1]是当前字符匹配的最长前缀
    int i = 1;
    while (i < m) {
        if (pattern[i] == pattern[len]) {
            len++;
            lps[i] = len;
            i++;
        } else { // 如果不匹配
            if (len != 0) {
                len = lps[len - 1];
            } else { // 如果len为0
                lps[i] = 0;
                i++;
            }
        }
    }
    return lps;
}

// KMP算法的搜索函数
std::vector<int> maxFlow::kmpSearch(const std::vector<int>& text, const std::vector<int>& pattern) {
    std::vector<int> lps = kmpProcessPattern(pattern);
    std::vector<int> occurrences;
    int i = 0; // text的索引
    int j = 0; // pattern的索引
    int n = text.size();
    int m = pattern.size();

    while (i < n) {
        if (pattern[j] == text[i]) {
            j++;
            i++;
        }
        if (j == m) {
            occurrences.push_back(i - j); // 完全匹配，记录pattern开始的索引
            j = lps[j - 1];
        } else if (i < n && pattern[j] != text[i]) {
            if (j != 0) {
                j = lps[j - 1];
            } else {
                i = i + 1;
            }
        }
    }
    return occurrences;
}

float maxFlow::pushMaxFlow(std::vector<int> &path,  std::vector<float> &nodeFlux){
    int n = path.size();
    float watch=ConjugateMatrix[17][24];
    std::vector<float> capacityleft;   // how many transcripts compatible to path enter node
    std::vector<float> capacityright;  // how many transcripts compatible to path exit node
    capacityleft.resize(n);
    capacityright.resize(n);
    std::vector<float> sumleft;        // how many transcripts enter node
    std::vector<float> sumright;       // how many transcripts exit node
    sumleft.resize(n);
    sumright.resize(n);
    std::map<int, bool> istranscript;  //true if transfreg support the current path
    /// init istranscript with false
    for(int i=0;i<graph->getTransFragment().size();i++) istranscript[i]=false;

    std::map<int, int> node2path;
    for(int i=0;i<n;i++) {
        node2path[path[i]]=i;
        nodeFlux.emplace_back(0.0);
    }
    nodeFlux.resize(n);
    // compute capacities and sums for all nodes
    for (int point_idx=0; point_idx < path.size(); point_idx++) {
        auto &point = path[point_idx];
        std::cout<<"For "<< point_idx<<"th node: "<<point<<std::endl;
        if (graph->paths_through_point.find(point) != graph->paths_through_point.end()) {
            // 获取穿过该点的所有路径
            for (const auto& pair : graph->paths_through_point[point]) {
                const auto& file_path = pair.first;
//                if (file_path.front()==graphInNode ){
//                    continue;
//                }
                // 如果穿过该点的路径是输入路径的子序列，则打印路径及其在文件中的索引
                if (kmpSearch(path, file_path).size()>0) {
                    // for path after first real node
                    if (point != file_path.back()) {
                        // transfrag ends after this node
                        capacityright[point_idx] += graph->getTransFragment()[pair.second].aboundance;
                        sumright[point_idx] += graph->getTransFragment()[pair.second].aboundance;
                    }
                    if (point != file_path.front()){
                        // transfrag starts before this node
                        capacityleft[point_idx] += graph->getTransFragment()[pair.second].aboundance;
                        sumleft[point_idx] += graph->getTransFragment()[pair.second].aboundance;
                    }
                    istranscript[pair.second] = true;
                    std::cout << "Path index " << pair.second << " supports the input path and passes through point " << point << std::endl;
                }
                else {
//                    if (point == file_path.back()) {
                        // transfrag ends after this node
                        sumright[point_idx] += graph->getTransFragment()[pair.second].aboundance;
//                    }
//                    if(point == file_path.front()){
                        // transfrag starts before this node
                        sumleft[point_idx] += graph->getTransFragment()[pair.second].aboundance;
//                    }
                    std::cout << "Path index " << pair.second << " does not support the input path and passes through point " << point << std::endl;
                }
            }
        }


    }

    { // DEBUG ONLY
        for(int i=1;i<n-1;i++) {
            fprintf(stderr,"debug Node %d LEFT: capacity=%f total=%f ",path[i],capacityleft[i],sumleft[i]);
            if(sumleft[i]) fprintf(stderr,"perc=%f ",capacityleft[i]/sumleft[i]);
            else fprintf(stderr,"perc=n/a ");
            fprintf(stderr,"RIGHT: capacity=%f total=%f ",capacityright[i],sumright[i]);
            if(sumright[i]) fprintf(stderr,"perc=%f\n",capacityright[i]/sumright[i]);
            else fprintf(stderr,"perc=n/a\n");
        }
        fprintf(stderr,"Used transcripts:");
        for(int i=0;i<graph->getTransFragment().size();i++) if(istranscript[i]) fprintf(stderr," %d(%f)",i,graph->getTransFragment()[i].aboundance);
        fprintf(stderr,"\n");
    }

    // compute flow
    for (int capli = 0; capli < capacityleft.size(); capli++){
        fprintf(stderr, "capacityleft[%d]: %f\n", capli, capacityleft[capli]);
    }
    for (int capri = 0; capri < capacityright.size(); capri++){
        fprintf(stderr, "capacityright[%d]: %f\n", capri, capacityright[capri]);
    }
    for (int sumli = 0; sumli < sumleft.size(); sumli++){
        fprintf(stderr, "sumleft[%d]: %f\n", sumli, sumleft[sumli]);
    }
    for (int sumri = 0; sumri < sumright.size(); sumri++){
        fprintf(stderr, "sumright[%d]: %f\n", sumri, sumright[sumri]);
    }

    float prevflow=capacityleft[1];
    for(int i=1;i<n-1;i++) {
        float percleft=prevflow/sumleft[i];
        float percright=capacityright[i]/sumright[i];
        if(percright>percleft) { // more transfrags leave node
            percright=percleft;
        }
        prevflow=percright*sumright[i];
    }
    if(!prevflow) return(0);

    for(int i=n-2;i>0;i--) {
        fprintf(stderr,"i=%d sumright=%f prevflow=%f\n",i,sumright[i],prevflow);
        nodeFlux[i]=prevflow/sumright[i];
        fprintf(stderr,"nodeflux=%f\n",nodeFlux[i]);
        capacityright[i]=prevflow;
        prevflow=nodeFlux[i]*sumleft[i];
        fprintf(stderr,"i=%d sumright=%f sumleft=%f prevflow=%f capacityright=%f nodeflux=%f\n",i,sumright[i],sumleft[i],prevflow,capacityright[i],nodeFlux[i]);
        capacityleft[i]=prevflow; // I don't use this
    }
    // then reduce support fragment abundances
    // then calculate how many abundances need reduces on each exon
    fprintf(stderr,"strat reduce aboundance here ..\n");
    for(int i=1;i<n-1;i++) if(capacityright[i]){
        if (path[i]== graphInNode && path[i]== graphOutNode)
            continue;
        int nt = graph->paths_through_point[path[i]].size();
        for(int j=0;j<nt;j++) {
            int transfrag_idx = graph->paths_through_point[path[i]][j].second;
            auto &cur_transf = graph->getTransFragment()[transfrag_idx];
            std::vector<int> &transfrag_path = cur_transf.paths;
            float transfrag_aboundance = cur_transf.aboundance;
            float trans_reducer = transfrag_aboundance;
            if (istranscript[transfrag_idx] && trans_reducer && (path[i] == transfrag_path[0] || path[i] == transfrag_path[1])){  // transfrag starts at this node
                if (i%2==0) continue; /// avoid reduce conjugate edge flow
                /// reduce the capacity of edge in conjugate matrix


                if (capacityright[i]> trans_reducer){
                    capacityright[i]-=trans_reducer;
                    int row, col;
//                    if (transfrag_path.front()!=graphInNode && transfrag_path.back()!=graphOutNode){
                        for (int pi=0; pi < transfrag_path.size()-1; pi++){
                            row=transfrag_path[pi];
                            col=transfrag_path[pi+1];
                            auto test0= ConjugateMatrix[row][col];
                            ConjugateMatrix[row][col] -= trans_reducer;
                            if (ConjugateMatrix[row][col] <= edge_epsilon){
                                ConjugateMatrix[row][col] = 0;
                            }
                        }
//                    }


                    if(capacityright[i]<epsilon) {
                        capacityright[i]=0;
                    }

                    // then reduce capacityright of next nodes in path, from current node to last node in current tranfrg(transfrag_path):
                    int idx_start_in_tranf;
                    for (int k=0; k<transfrag_path.size(); k++){
                        if (path[i] == transfrag_path[k]){
                            idx_start_in_tranf = k;
                            break;
                        }
                    }
//                    if (path[i] == transfrag_path[0]){
//                        idx_start_in_tranf = 0;
//                    }else if(path[i] == transfrag_path[1])
//                        idx_start_in_tranf = 1;
                    int n_remain = transfrag_path.size()-1-idx_start_in_tranf;
                    /// if i is the conjugate index(even number), continue;
                    for (int k=i+1; k < i+n_remain; k++){
                        capacityright[k] -= trans_reducer;
                    }
                    cur_transf.aboundance-=trans_reducer;
                    if(cur_transf.aboundance<epsilon) {
                        cur_transf.aboundance=0;
                    }
                }
                else{
                    cur_transf.aboundance-=capacityright[i];
                    if(cur_transf.aboundance<epsilon) {
                        cur_transf.aboundance=0;
                    }
                    int row, col;
//                    if (transfrag_path.front()!=graphInNode && transfrag_path.back()!=graphOutNode) {
                        for (int pi=0; pi < transfrag_path.size()-1; pi++){
                            row=transfrag_path[pi];
                            col=transfrag_path[pi+1];
                            auto test= ConjugateMatrix[row][col];
                            ConjugateMatrix[row][col] -= capacityright[i];
                            if (ConjugateMatrix[row][col] <= edge_epsilon){
                                ConjugateMatrix[row][col] = 0;
                            }
//                        }

                    }

                    capacityright[i]-=trans_reducer;
                    if(capacityright[i]<epsilon) {
                        capacityright[i]=0;
                    }
                    // then reduce capacityright of next nodes in path, from current node to last node in current tranfrg(transfrag_path):
                    int idx_start_in_tranf;
                    // get the index of path[i] in transfrag_path
                    for (int k=0; k<transfrag_path.size(); k++){
                        if (path[i] == transfrag_path[k]){
                            idx_start_in_tranf = k;
                            break;
                        }
                    }
//                    if (path[i] == transfrag_path[0]){
//                        idx_start_in_tranf = 0;
//                    }else if(path[i] == transfrag_path[1])
//                        idx_start_in_tranf = 1;
                    int n_remain = transfrag_path.size()-1-idx_start_in_tranf;
                    /// if i is the conjugate index(even number), continue;
                    for (int k=i+1; k < i+n_remain; k++){
                        capacityright[k] -= trans_reducer;
                    }
                }
            }
        }
    }
    return nodeFlux[1];
}

std::map<uint, std::string> maxFlow::GetmatrixIdx2nodeIdxStr() {
    return matrixIdx2nodeIdxStr;
}



