#pragma once
#include <vector>
#include <string.h>

using namespace std;

typedef vector<vector<int>> intVectors;


struct Graph {
private:
	Graph();
    Graph(int numberOfVertices);
	
    vector<string> nodes;
    intVectors weights;

public:
    static Graph* mapGraphFromFile(const string& str);
	static bool validateUnsignedInt(const int number);
	
	
	const vector<string> getNodes() const { return nodes; }
    const intVectors getWeights() const { return weights; }
	void setNodes(int numberOfVertices);
	int getNumberOfEdges() const { return weights.size(); }
	
    void printAdjacencyMatrix() const;
	
	static const int noEdge = -1;
	static const char delimiter=' ';
	static const char noWeight = '0';
};
