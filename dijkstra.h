#pragma once

#include <string>

#include "Graph.h"

void dijkstra(const Graph& m, const std::string& initialNodeName, const std::string& goalNodeName, const int mpiNodesCount);

void dijkstraWorker(int mpiNodeId, int mpiWorkerNodes);
