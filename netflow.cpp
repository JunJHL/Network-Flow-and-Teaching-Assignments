#include "netflow.hpp"
#include <queue>
#include <iostream>
#include <vector>
#include <climits>
#include <algorithm>

#define TOO_FEW_VERTICES "Too few vertices."
#define TOO_FEW_EDGES "Too few edges."
#define EDGE_WEIGHT_ZERO "Detected edge weight of 0."
#define EDGE_BAD_ENDPOINT "Edge interacts with nonexistent vertex."
#define SELF_LOOP "At least one self-loop."
#define MULTI_EDGES "Detected multi-edges."
#define NOT_ONE_SRC "Zero or more than one source."
#define NOT_ONE_SINK "Zero or more than one sink."

bool BFS(unsigned f, unsigned t, std::vector<Edge>& rGraph, std::vector<int>& parent){
    Edge edge;

    std::vector<bool> discovered(rGraph.size());    // num edges discovered
    for(unsigned i = 0; i < discovered.size(); i++) {
        discovered[i] = false;
    }

    std::queue<Edge> q;
    for(unsigned i = 0; i < rGraph.size(); i++) {
        if(rGraph[i].from == f) {
            q.push(rGraph[i]);
            discovered[i] = true;
        }
    }
    parent[f] = -1;

    while(!q.empty()) {
        edge = q.front();
        q.pop();
        for(unsigned i = 0; i < rGraph.size(); i++) {
            if(edge.to == t && !discovered[i] && edge.weight > 0) {
                parent[t] = edge.from;
                return true;
            }
            // a) iterate all possible next edge of current edge
            // b) check if current edge is seen
            // c) check if weight of current edge is > 0
            if(edge.to == rGraph[i].from && !discovered[i] && edge.weight > 0) {    // if current node is the to node in the edge to find the neighbor node of vertex
                discovered[i] = true;
                parent[rGraph[i].from] = (int)edge.from;
                q.push(rGraph[i]);
            }
        }
    }
    return false;
}

std::vector<Edge> solveNetworkFlow(
    const std::vector<Edge>& flowNetwork,
    unsigned numVertices)
{
    unsigned source;
    unsigned sink;
    unsigned src_count = 0;
    unsigned sink_count = 0;
    std::vector<bool> is_source(numVertices, true);
    std::vector<bool> is_sink(numVertices, true);

    std::vector<Edge> max_flow(flowNetwork.size());
    std::vector<Edge> rGraph(flowNetwork.size());
    std::vector<int> parent(numVertices);

    if(numVertices < 2) {
        throw std::runtime_error(TOO_FEW_VERTICES);
    }else if(flowNetwork.size() == 0) {
        throw std::runtime_error(TOO_FEW_EDGES);
    }

    for(unsigned i = 0; i < flowNetwork.size(); i++) {
        if(flowNetwork[i].weight == 0) {
            throw std::runtime_error(EDGE_WEIGHT_ZERO);
        }else if(flowNetwork[i].from >= numVertices || flowNetwork[i].to >= numVertices) {
            throw std::runtime_error(EDGE_BAD_ENDPOINT);
        }else if(flowNetwork[i].from == flowNetwork[i].to) {
            throw std::runtime_error(SELF_LOOP);
        }

        //  set all nodes that is not source/to to false
        is_source[flowNetwork[i].to] = false;
        is_sink[flowNetwork[i].from] = false;

        for(unsigned  j = i+1; j < flowNetwork.size(); j++) {
            if((flowNetwork[i].from == flowNetwork[j].from && flowNetwork[i].to == flowNetwork[j].to) && (flowNetwork[i].weight != flowNetwork[j].weight)) {
                throw std::runtime_error(MULTI_EDGES);
            }
        }
    }

    for(unsigned i = 0; i < is_source.size(); i++) {
        if(is_source[i]) {
            ++src_count; 
            source = i;
        }
    }
    if(src_count != 1) {
        throw std::runtime_error(NOT_ONE_SRC);
    }
    for(unsigned i = 0; i < is_sink.size(); i++) {
        if(is_sink[i]) {
            ++sink_count; 
            sink = i;
        }
    }
    if(sink_count != 1) {
        throw std::runtime_error(NOT_ONE_SINK);
    }

    for(unsigned i = 0; i < flowNetwork.size(); i++) {   // create max_flow graph
        max_flow[i] = flowNetwork[i];
    }

    for(unsigned i = 0; i < flowNetwork.size(); i++) {   // create residual graph
        rGraph[i] = flowNetwork[i];
        rGraph.push_back({flowNetwork[i].to, flowNetwork[i].from, 0});
    }

    // BFS(source, sink, rGraph, parent);
    
    while(BFS(source, sink, rGraph, parent)) {  // while path to sink != false
        unsigned bottle_neck = UINT_MAX;
        int node = (int)sink;
        while(parent[node] != -1) {     // find path location in rGraph
            for(unsigned i = 0; i < rGraph.size(); i++) {
                if(parent[node] == (int)rGraph[i].from && node == (int)rGraph[i].to) {
                    bottle_neck = std::min(bottle_neck, rGraph[i].weight);     // look for bottleneck
                    node = parent[node];
                    break;
                }
            }
        }

        node = (int)sink;
        while(parent[node] != -1) {     // find path location in rGraph
            for(unsigned i = 0; i < rGraph.size(); i++) {
                if(parent[node] == (int)rGraph[i].from && node == (int)rGraph[i].to) {  // update edge weight in rGraph
                    rGraph[i].weight -= bottle_neck;
                    node = parent[node];    // update node location
                }
                if(parent[node] == (int)rGraph[i].to && node == (int)rGraph[i].from) {
                    rGraph[i].weight += bottle_neck;   // since this is flow back, no need to update node location
                }
            }
        }
    }

    for(unsigned i = 0; i < max_flow.size(); i++) {
        for(unsigned j = 0; j < rGraph.size(); j++) {
            if(flowNetwork[i].from == rGraph[j].from && flowNetwork[i].to == rGraph[j].to) {
                max_flow[i].weight = flowNetwork[i].weight - rGraph[i].weight;
                break;
            }
        }
    }

    return max_flow;
}

void assignCourses(
        std::vector<Instructor>& instructors,
        const std::vector<std::string>& courses)
{
    std::vector<Edge> flowNetwork;

    for(unsigned i = 1; i <= instructors.size(); i++) {
        flowNetwork.push_back({0, i, instructors[i-1].maxCourses});
        for(unsigned j = 0; j < instructors[i-1].preferences.size(); j++) {
            unsigned k = (unsigned)(std::find(courses.begin(), courses.end(), instructors[i-1].preferences[j]) - courses.begin() + 1);
            flowNetwork.push_back({i, (unsigned)instructors.size() + k, 1});
            flowNetwork.push_back({(unsigned)instructors.size() + k, (unsigned)(instructors.size() + courses.size()) + 1, 1});
            // std::cout << i << "->" << instructors.size() + k << ": " << courses[k-1] << std::endl;
            // std::cout << instructors.size() + k << "->" << instructors.size() + courses.size() + 1 << std::endl;
        }
    }
    
    for(unsigned i = 0; i < flowNetwork.size(); i++) {
        for(unsigned j = i + 1; j < flowNetwork.size(); j++) {  //  doesn't need to check back because front is checked
            if(flowNetwork[i].from == flowNetwork[j].from && flowNetwork[i].to == flowNetwork[j].to) {
                flowNetwork.erase(flowNetwork.begin() + j);    // remove duplicate edges
                --j;    // stay if erased
            }
        }
    }
    std::vector<Edge> maxFlow = solveNetworkFlow(flowNetwork, (unsigned)(instructors.size() + courses.size() + 2));
    
    for(unsigned i = 0; i < maxFlow.size(); i++) {
        if(maxFlow[i].from == 0) {  // if instructors
            for(unsigned j = 0; j < maxFlow.size(); j++) {
                if(maxFlow[i].to == maxFlow[j].from && maxFlow[j].weight > 0) {
                    instructors[maxFlow[i].to - 1].assignedCourses.push_back(courses[maxFlow[j].to - instructors.size() - 1]);
                }
            }
        }
    }

    for(unsigned i = 0; i < instructors.size(); i++) {
        for(unsigned j = 0; j < instructors.size(); j++) {
            for(unsigned k = 0; k < instructors[i].assignedCourses.size(); k++) {
                auto pos = std::find(instructors[j].assignedCourses.begin(), instructors[j].assignedCourses.end(), instructors[i].assignedCourses[k]);
                if(pos != instructors[j].assignedCourses.end() && i != j) {
                    if(instructors[i].maxCourses > instructors[j].maxCourses) {
                        instructors[i].assignedCourses.erase(instructors[i].assignedCourses.begin() + k);
                    }else {
                        instructors[j].assignedCourses.erase(pos);
                    }
                }
            } 
        }
    }
    
    // for (const Edge& edge : maxFlow)
    //     std::cout << edge.from << " -> " << edge.to
    //         << " (" << edge.weight<< ")\n";

}
