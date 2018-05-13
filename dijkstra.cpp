#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <mpi.h>

#include "Graph.h"

const int MAX_INT = 2147483647;
const int rootID = 0;

std::pair<int, int> getMpiWorkerNodeRanges(int nodesCount, int mpiNodesCount, int mpiNodeId) {
    mpiNodesCount -= 1; // counting only worker nodes
    int fromNode = (nodesCount / mpiNodesCount) * (mpiNodeId - 1);
    int toNode = (nodesCount / mpiNodesCount) * mpiNodeId - 1;
    int restNodes = nodesCount % mpiNodesCount;

    if (mpiNodeId - 1 < restNodes) {
        fromNode += mpiNodeId - 1;
        toNode += mpiNodeId;
    }
    else {
        fromNode += restNodes;
        toNode += restNodes;
    }

    return std::pair<int, int>(fromNode, toNode);
}

void dijkstraMain(const Graph *graph, const std::string& initialNodeName, const std::string& goalNodeName, const int mpiNodesCount) {
    const intVectors weights = graph->getWeights();
    const vector<string> nodes = graph->getNodes();
    auto nodesCount = nodes.size();

    std::vector<int> distances(nodesCount);
    std::vector<int> prevNodes(nodesCount);

    std::set<int> visited;

    for(auto node = 0u; node < nodesCount; ++node) {
        distances[node] = MAX_INT;
        prevNodes[node] = MAX_INT;
    }

    auto indexOf = [&] (auto nodeName) { return std::find(nodes.begin(), nodes.end(), nodeName) - nodes.begin(); };
    auto isVisited = [&] (auto node) { return visited.find(node) != visited.end(); };
    auto isNeighbour = [&] (auto currentNode, auto node) { return weights[currentNode][node] != -1; };

    auto initialNode = static_cast<int>(indexOf(initialNodeName));
    auto currentNode = initialNode;
    auto goalNode = indexOf(goalNodeName);

    int workerNodes = mpiNodesCount - 1;
    int nodesStep = nodesCount / workerNodes;

    std::cout << "Sending initial data to workers..." << std::endl;
    int data[3] = {nodesCount, initialNode, goalNode};
    MPI_Bcast(&data, 3, MPI_INT, rootID, MPI_COMM_WORLD);

    for(auto i=0u; i<nodesCount; ++i)
        MPI_Bcast((int*)&weights[i][0], nodesCount, MPI_INT, rootID, MPI_COMM_WORLD);

    distances[initialNode] = 0;

    while (true) {
        std::cout << "Sending currNode=" << currentNode << " distance=" << distances[currentNode] << std::endl;
        int data[2] = {currentNode, distances[currentNode]};
        MPI_Bcast(&data, 2, MPI_INT, rootID, MPI_COMM_WORLD);

        //Get distances calculated by workers
        for(auto mpiNodeId = 1; mpiNodeId < mpiNodesCount; ++mpiNodeId) {
            const auto nodeRanges = getMpiWorkerNodeRanges(nodesCount, mpiNodesCount, mpiNodeId);
            const auto fromNode = nodeRanges.first;
            const auto toNode = nodeRanges.second;

            MPI_Recv(&distances[fromNode], toNode - fromNode + 1, MPI_INT, mpiNodeId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::cout << "Recv from " << mpiNodeId << std::endl;
        }

        // test for goal
        if (currentNode == goalNode) {
            std::cout << "Goal node found" << std::endl;
            for(auto mpiNodeId=1; mpiNodeId<mpiNodesCount; ++mpiNodeId) {
                const auto nodeRanges = getMpiWorkerNodeRanges(nodesCount, mpiNodesCount, mpiNodeId);
                const auto fromNode = nodeRanges.first;
                const auto toNode = nodeRanges.second;
                MPI_Recv(&prevNodes[fromNode], toNode - fromNode + 1, MPI_INT, mpiNodeId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }


            std::vector<int> stack;
            while(currentNode != initialNode) {
                stack.push_back(currentNode);

                auto prev = prevNodes[currentNode];
                currentNode = prev;
            }


            std::cout << "Total cost: " << distances[goalNode] << std::endl;
            std::cout << "Path: " << nodes[currentNode];

            for(auto it=stack.rbegin(); it != stack.rend(); ++it) {
                auto nextNodeName = nodes[*it];
                std::cout << " -> " << nextNodeName;
            }
            std::cout << std::endl;

            break;
        }

		std::cout << "Visited " << currentNode;
        visited.insert(currentNode);

        auto minCost = MAX_INT; 
        auto nextNode = -1;

        for(auto node = 0u; node<nodesCount; ++node) {
            auto totalCost = distances[node];

            if (!isVisited(node) && totalCost < minCost) {
                minCost = totalCost;
                nextNode = node;
            }
        }
        
        currentNode = nextNode;
		std::cout << "Next currentNode = " << currentNode;

        if (currentNode == -1) {
            std::cout << "Path not found" << std::endl;
            return;
        }
    }
}


void dijkstraNode(int mpiNodeId, int mpiNodesCount) {
    int data[3];
    MPI_Bcast(&data, 3, MPI_INT, 0, MPI_COMM_WORLD);

    int nodesCount = data[0];
    int initialNode = data[1];
    int goalNode = data[2];

    std::cout << "mpiID=" << mpiNodeId << " bcast recv nodesCount=" << nodesCount << " initialNode=" << initialNode << " goalNode=" << goalNode << std::endl;

    std::vector<std::vector<int>> weights(nodesCount);
    std::vector<int> distances(nodesCount, MAX_INT);
    std::vector<int> prevNodes(nodesCount, MAX_INT);
    std::set<int> visited;

    auto isVisited = [&] (auto node) { return visited.find(node) != visited.end(); };
    auto isNeighbour = [&] (auto currentNode, auto node) { return weights[currentNode][node] != -1; };

    for(auto i=0u; i<nodesCount; ++i) {
        weights[i].resize(nodesCount);
        MPI_Bcast((int*)&weights[i][0], nodesCount, MPI_INT, rootID, MPI_COMM_WORLD);
    }

    const auto nodeRanges = getMpiWorkerNodeRanges(nodesCount, mpiNodesCount, mpiNodeId);
    const auto fromNode = nodeRanges.first;
    const auto toNode = nodeRanges.second;
    std::cout << "mpiID=" << mpiNodeId << " from=" << fromNode << " to=" << toNode << std::endl;

    // real work
    while (1) {
        std::cout << "~~~~" << std::endl;
        MPI_Bcast(&data, 2, MPI_INT, rootID, MPI_COMM_WORLD);

        int currentNode = data[0];
        distances[currentNode] = data[1];
        std::cout << "mpiId=" << mpiNodeId << " bcast recv currNode=" << currentNode << " dist=" << distances[currentNode] << std::endl;

        // goal node not found
        if (currentNode == -1)
            return;

        for (auto node=fromNode; node<=toNode; ++node) {
            if (isVisited(node)) {
				std::cout << "Node " << node << " already visited - skipping.";
                continue;
            }

            if (isNeighbour(currentNode, node)) {
                auto nodeDistance = weights[currentNode][node];
                auto totalCostToNode = distances[currentNode] + nodeDistance;

				std::cout << "Node " << node << " is neighbour of " << currentNode << " (distance: " << nodeDistance << ", totalCostToNode: " << totalCostToNode << ")";
                if (totalCostToNode < distances[node]) {
                    distances[node] = totalCostToNode;
                    prevNodes[node] = currentNode;
					std::cout << "New total cost is less than the old, replacing";
                }
            }
        }

        visited.insert(currentNode);
        MPI_Send(&distances[fromNode], toNode - fromNode + 1, MPI_INT, rootID, 0, MPI_COMM_WORLD);

        // goal node found
        if (currentNode == goalNode) {
            MPI_Send(&prevNodes[fromNode], toNode - fromNode + 1, MPI_INT, rootID, 0, MPI_COMM_WORLD);
            return;
        }
        std::cout << "======" << std::endl;
    }
}
