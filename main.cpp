#include <iostream>
#include <mpi.h>

#include "Graph.h"
#include "dijkstra.cpp"


int main(int argc, char* argv[]) {
	int mpiNodesCount;
	int mpiNodeId;
    const int rootID = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiNodesCount);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiNodeId);

    if (mpiNodeId == rootID) {
        Graph *graph = Graph::mapGraphFromFile(argv[1]);
       
		dijkstraMain(graph, *graph->getNodes().begin(), *(graph->getNodes().end() - 1), mpiNodesCount);
    }
	else {
		dijkstraNode(mpiNodeId, mpiNodesCount);
	}


    MPI_Finalize();

    return 0;
}
