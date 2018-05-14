#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <mpi.h>

#include "Graph.h"

const int MAX_INT = 2147483647;
const int rootId = 0;

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
    const intVectors graphWeights = graph->getWeights();
    const vector<string> vectorOfNodes = graph->getNodes();
    const int nodesCount = vectorOfNodes.size();
	const int initialNode = getNodeIndex(vectorOfNodes, initialNodeName);
	const int goalNode = getNodeIndex(vectorOfNodes, goalNodeName);

	vector<int> vectorOfDistancesToEachNode(nodesCount, MAX_INT);
	vector<int> vectorOfPreviousVisitedNodes(nodesCount, MAX_INT);

    set<int> setOfVisitedNodes;


    int currentNode = initialNode;



    cout << "Sending initial data to workers..." << endl;
    int buffer[3] = {nodesCount, initialNode, goalNode};
    MPI_Bcast(&buffer, 3, MPI_INT, rootId, MPI_COMM_WORLD);

    for(int i = 0; i < nodesCount; ++i)
        MPI_Bcast((int*)&graphWeights[i][0], nodesCount, MPI_INT, rootId, MPI_COMM_WORLD);

    vectorOfDistancesToEachNode[initialNode] = 0;

    while(true) {
        cout << "Sending currNode=" << currentNode << " distance=" << vectorOfDistancesToEachNode[currentNode] << endl;
        int buffer[2] = {currentNode, vectorOfDistancesToEachNode[currentNode]};
        MPI_Bcast(&buffer, 2, MPI_INT, rootId, MPI_COMM_WORLD);

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

            vector<int> stackOfNodes;
            while(currentNode != initialNode) {
                stackOfNodes.push_back(currentNode);

                int prev = vectorOfPreviousVisitedNodes[currentNode];
                currentNode = prev;
            }

            cout << "Total cost: " << vectorOfDistancesToEachNode[goalNode] << endl;

			cout << "Shortest path to each node from initial: ";
			for (int d : vectorOfDistancesToEachNode) {
				cout << d << " ";
			}
			cout << endl;

			cout << "Result path: " << vectorOfNodes[currentNode];
			reverse(stackOfNodes.begin(), stackOfNodes.end());
			for (int node : stackOfNodes) {
				cout << " -> " << vectorOfNodes[node];
			}
            cout << endl;

            break;
        }

        setOfVisitedNodes.insert(currentNode);

        int minCost = MAX_INT; 
        int nextNode = -1;

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
    int buffer[3];
    MPI_Bcast(&buffer, 3, MPI_INT, rootId, MPI_COMM_WORLD);

    int nodesCount = buffer[0];
    int initialNode = buffer[1];
    int goalNode = buffer[2];

	intVectors graphWeights(nodesCount);
    vector<int> vectorOfDistancesToEachNode(nodesCount, MAX_INT);
    vector<int> vectorOfPreviousVisitedNodes(nodesCount, MAX_INT);
    set<int> setOfVisitedNodes;

    for(int i = 0; i < nodesCount; ++i) {
        graphWeights[i].resize(nodesCount);

        MPI_Bcast((int*)&graphWeights[i][0], nodesCount, MPI_INT, rootId, MPI_COMM_WORLD);
    }

    const pair<int, int> nodeRangePair = calculateNodeRange(nodesCount, mpiNodesCount, mpiNodeId);
    const int fromNode = nodeRangePair.first;
    const int toNode = nodeRangePair.second;
    cout << "mpiID=" << mpiNodeId << " from=" << fromNode << " to=" << toNode << endl;

    while(true) {
        MPI_Bcast(&buffer, 2, MPI_INT, rootId, MPI_COMM_WORLD);

        int currentNode = buffer[0];
        vectorOfDistancesToEachNode[currentNode] = buffer[1];
        cout << "mpiId=" << mpiNodeId << " bcast recv currNode=" << currentNode << " dist=" << vectorOfDistancesToEachNode[currentNode] << endl;

        if (currentNode == -1)
            return;

        for (int node = fromNode; node <= toNode; ++node) {
            if (isNodeVisited(setOfVisitedNodes, node)) {
                continue;
            }

            if (isCurrentNodeNeighbour(graphWeights, currentNode, node)) {
                int nodeDistance = graphWeights[currentNode][node];
                int totalCostToNode = vectorOfDistancesToEachNode[currentNode] + nodeDistance;

				cout << "Node " << node << " is neighbour of " << currentNode << " (distance: " << nodeDistance << ", totalCostToNode: " << totalCostToNode << ")";
                if (totalCostToNode < vectorOfDistancesToEachNode[node]) {
                    vectorOfDistancesToEachNode[node] = totalCostToNode;
                    vectorOfPreviousVisitedNodes[node] = currentNode;
					cout << "New total cost is less than the old, replacing";
                }
            }
        }

        setOfVisitedNodes.insert(currentNode);
        MPI_Send(&vectorOfDistancesToEachNode[fromNode], toNode - fromNode + 1, MPI_INT, rootId, 0, MPI_COMM_WORLD);

        // goal node found
        if (currentNode == goalNode) {
            MPI_Send(&vectorOfPreviousVisitedNodes[fromNode], toNode - fromNode + 1, MPI_INT, rootId, 0, MPI_COMM_WORLD);
            return;
        }
    }
}