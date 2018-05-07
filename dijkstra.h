#pragma once

#include <string>

#include "map.h"

const int INF = 2147483647;

void dijkstra(const Map& m, const std::string& initialNodeName, const std::string& goalNodeName, const int mpiNodesCount);

void dijkstraWorker(int mpiNodeId, int mpiWorkerNodes);
