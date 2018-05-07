#pragma once
#include <vector>

using namespace std;

const int NO_EDGE = -1;

typedef vector<vector<int>> intVectors;


struct Graph {
private:
    Graph(int verticesCount);

    vector<string> nodes;
    intVectors weights;

public:
    static Graph fromFile(const string& str, const char delimiter=',');
    static Graph fromFile(ifstream&& istream, const char delimiter=',');

    int getSize() const { return weights.size(); }
    const decltype(weights) getWeights() const { return weights; }

    const decltype(nodes) getNodes() const { return nodes; }

    void printWeights() const;
};
