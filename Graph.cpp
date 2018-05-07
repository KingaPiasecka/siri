#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>

#include "Graph.h"

Graph::Graph(int verticesCount) : weights(verticesCount) {
    for(int i=0; i<verticesCount; ++i) {
        weights[i].resize(verticesCount);
    }
}

Graph Graph::fromFile(const string& str, const char delimiter) {
    return fromFile(ifstream(str.c_str()), delimiter);
}

Graph Graph::fromFile(ifstream&& istream, const char delimiter) {
    string header;

    if (getline(istream, header, ':') && header == "vertices" && getline(istream, header)) {
        int verticesCount = stoi(header);
        
        Graph m(verticesCount);        

        for(auto i=0; i<verticesCount; ++i) {
            string line;
            getline(istream, line); 

            stringstream linestream(line);
            for(auto j=0; j<verticesCount; ++j) {
                string numstr;
                getline(linestream, numstr, delimiter);
                
                // the `i` node has no edge to `j`
                if (numstr.find("-") != string::npos)
                    m.weights[i][j] = NO_EDGE;
                    
                else
                    m.weights[i][j] = stoi(numstr);
            }
        }

        for(auto i=0; i<verticesCount; ++i) {
            char nodeName[2] = {0};
            nodeName[0] = (char) ( (int)('A') + i );
            cout << "Adding node: " << nodeName << endl;
            m.nodes.push_back(string(nodeName));
        }

        return m;
    }

    throw runtime_error("Wrong Graph file format");
}

void Graph::printWeights() const {
    for(auto i=0u; i<getSize(); ++i) {
        for(auto j=0u; j<getSize(); ++j)
            if (weights[i][j] != -1)
                cout << setw(4) << weights[i][j] << " ";
            else
                cout << setw(4) << "-" << " ";

        cout << endl;
    }
}
