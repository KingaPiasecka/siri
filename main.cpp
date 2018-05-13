#include <iostream>
#include <mpi.h>

#include "Graph.h"
#include "dijkstra.cpp"

int main(int argc, char* argv[]) {
    if (argc == 1) {
        std::cout << "Usage: " << argv[0] << " <testcase file>" << std::endl;
        return -1;
    }

	int mpiNodesCount;
	int mpiNodeId;
    const int rootID = 0;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiNodesCount);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiNodeId);
    std::cout << "Hello from " << mpiNodeId << std::endl;

    if (mpiNodeId == rootID) {
        Graph *graph = Graph::mapGraphFromFile(argv[1]);
        
        auto n = graph->getNodes();
        auto initialNodeName = *n.begin();
        auto goalNodeName = *(n.end()-1);

        std::cout << "Dijkstra search algorithm" << std::endl;
        std::cout << "Starting at node: " << initialNodeName << std::endl;
        std::cout << "Ending at node: " << goalNodeName << std::endl;
        std::cout << "Consecutive nodes (A, B, ...) weights: " << std::endl;
        graph->printAdjacencyMatrix();

        std::cout << "Searching..." << std::endl;
        dijkstraMain(graph, initialNodeName, goalNodeName, mpiNodesCount);
    }
    else
        dijkstraNode(mpiNodeId, mpiNodesCount);

    MPI_Finalize();

    return 0;
}
