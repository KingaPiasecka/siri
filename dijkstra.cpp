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

using namespace std;

pair<int, int> calculateNodeRange(int nodesCount, int mpiNodesCount, int mpiNodeId) {
    mpiNodesCount -= 1;
    int fromNode = (nodesCount / mpiNodesCount) * (mpiNodeId - 1);
    int toNode = (nodesCount / mpiNodesCount) * mpiNodeId - 1;
    int otherNodesCount = nodesCount % mpiNodesCount;

    if (mpiNodeId - 1 < otherNodesCount) {
        fromNode += mpiNodeId - 1;
        toNode += mpiNodeId;
    } else {
        fromNode += otherNodesCount;
        toNode += otherNodesCount;
    }

    return pair<int, int>(fromNode, toNode);
}


int getNodeIndex(vector<string> nodes, string nodeName) {
	return find(nodes.begin(), nodes.end(), nodeName) - nodes.begin();
}

bool isNodeVisited(set<int> visited, int node) {
	return visited.find(node) != visited.end();
}

bool isCurrentNodeNeighbour(intVectors weights, int currentNode, int node) { 
	return weights[currentNode][node] != -1; 
}

void dijkstraMain(const Graph *graph, const string& initialNodeName, const string& goalNodeName, const int mpiNodesCount) {
    const intVectors weights = graph->getWeights();
    const vector<string> nodes = graph->getNodes();
    const int nodesCount = nodes.size();

	vector<int> vectorOfDistancesToEachNode(nodesCount, MAX_INT);
	vector<int> vectorOfPreviousVisitedNodes(nodesCount, MAX_INT);

    set<int> setOfVisitedNodes;


    auto initialNode = static_cast<int>(getNodeIndex(nodes, initialNodeName));
    auto currentNode = initialNode;
    auto goalNode = getNodeIndex(nodes, goalNodeName);

    int workerNodes = mpiNodesCount - 1;

    cout << "Sending initial data to workers..." << endl;
    int data[3] = {nodesCount, initialNode, goalNode};
    MPI_Bcast(&data, 3, MPI_INT, rootID, MPI_COMM_WORLD);

    for(auto i=0u; i<nodesCount; ++i)
        MPI_Bcast((int*)&weights[i][0], nodesCount, MPI_INT, rootID, MPI_COMM_WORLD);

    vectorOfDistancesToEachNode[initialNode] = 0;

    while (true) {
        cout << "Sending currNode=" << currentNode << " distance=" << vectorOfDistancesToEachNode[currentNode] << endl;
        int data[2] = {currentNode, vectorOfDistancesToEachNode[currentNode]};
        MPI_Bcast(&data, 2, MPI_INT, rootID, MPI_COMM_WORLD);

        for(int mpiNodeId = 1; mpiNodeId < mpiNodesCount; ++mpiNodeId) {
            const pair<int, int> nodeRanges = calculateNodeRange(nodesCount, mpiNodesCount, mpiNodeId);
			const int fromNode = nodeRanges.first;
			const int toNode = nodeRanges.second;

            MPI_Recv(&vectorOfDistancesToEachNode[fromNode], toNode - fromNode + 1, MPI_INT, mpiNodeId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (currentNode == goalNode) {
            for(int mpiNodeId = 1; mpiNodeId < mpiNodesCount; ++mpiNodeId) {
                const pair<int, int> nodeRanges = calculateNodeRange(nodesCount, mpiNodesCount, mpiNodeId);
                const int fromNode = nodeRanges.first;
                const int toNode = nodeRanges.second; 
                MPI_Recv(&vectorOfPreviousVisitedNodes[fromNode], toNode - fromNode + 1, MPI_INT, mpiNodeId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }


            vector<int> stack;
            while(currentNode != initialNode) {
                stack.push_back(currentNode);

                int prev = vectorOfPreviousVisitedNodes[currentNode];
                currentNode = prev;
            }


            cout << "Total cost: " << vectorOfDistancesToEachNode[goalNode] << endl;

			cout << "Shortest path to each node from initial: ";
			for (int d : vectorOfDistancesToEachNode) {
				cout << d << " ";
			}
			cout << endl;


			cout << "Result path: " << nodes[currentNode];
			reverse(stack.begin(), stack.end());
			for (int node : stack) {
				cout << " -> " << nodes[node];
			}

            cout << endl;

            break;
        }

        setOfVisitedNodes.insert(currentNode);

        auto minCost = MAX_INT; 
        auto nextNode = -1;

        for(int node = 0; node < nodesCount; ++node) {
            int totalCost = vectorOfDistancesToEachNode[node];

            if (!isNodeVisited(setOfVisitedNodes, node) && totalCost < minCost) {
                minCost = totalCost;
                nextNode = node;
            }
        }
        
        currentNode = nextNode;
    }
}


void dijkstraNode(int mpiNodeId, int mpiNodesCount) {
    int data[3];
    MPI_Bcast(&data, 3, MPI_INT, 0, MPI_COMM_WORLD);

    int nodesCount = data[0];
    int initialNode = data[1];
    int goalNode = data[2];

	intVectors weights(nodesCount);
    vector<int> distances(nodesCount, MAX_INT);
    vector<int> prevNodes(nodesCount, MAX_INT);
    set<int> visited;

    for(int i = 0; i < nodesCount; ++i) {
        weights[i].resize(nodesCount);
        MPI_Bcast((int*)&weights[i][0], nodesCount, MPI_INT, rootID, MPI_COMM_WORLD);
    }

    const pair<int, int> nodeRangePair = calculateNodeRange(nodesCount, mpiNodesCount, mpiNodeId);
    const int fromNode = nodeRangePair.first;
    const int toNode = nodeRangePair.second;
    cout << "mpiID=" << mpiNodeId << " from=" << fromNode << " to=" << toNode << endl;

    // real work
    while (1) {
        MPI_Bcast(&data, 2, MPI_INT, rootID, MPI_COMM_WORLD);

        int currentNode = data[0];
        distances[currentNode] = data[1];
        cout << "mpiId=" << mpiNodeId << " bcast recv currNode=" << currentNode << " dist=" << distances[currentNode] << endl;

        // goal node not found
        if (currentNode == -1)
            return;

        for (int node = fromNode; node <= toNode; ++node) {
            if (isNodeVisited(visited, node)) {
                continue;
            }

            if (isCurrentNodeNeighbour(weights, currentNode, node)) {
                auto nodeDistance = weights[currentNode][node];
                auto totalCostToNode = distances[currentNode] + nodeDistance;

				cout << "Node " << node << " is neighbour of " << currentNode << " (distance: " << nodeDistance << ", totalCostToNode: " << totalCostToNode << ")";
                if (totalCostToNode < distances[node]) {
                    distances[node] = totalCostToNode;
                    prevNodes[node] = currentNode;
					cout << "New total cost is less than the old, replacing";
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
    }
}
