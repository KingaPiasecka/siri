#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>

#include "Graph.h"

Graph::Graph(){}

Graph::Graph(int numberOfVertices) : weights(numberOfVertices) {
    for(int i=0; i<numberOfVertices; ++i) {
        weights[i].resize(numberOfVertices);
    }
}

Graph* Graph::mapGraphFromFile(const string& fileName) {
	string header;
	string line;
	
	std::ifstream ifs (fileName.c_str());
	
	getline(ifs, header);
	
	try {
		int numberOfVertices = stoi(header);
		
		if(Graph::validateUnsignedInt(numberOfVertices)) {
			
			Graph *graph = new Graph(numberOfVertices);
			
			for(int i=0; i<numberOfVertices; ++i) {
				getline(ifs, line); 
				stringstream linestream(line);
				
				for(int j=0; j<numberOfVertices; ++j) {
					string numstr;
					getline(linestream, numstr, Graph::delimiter);
					if (numstr.find(Graph::noWeight) != string::npos) {
						graph->weights[i][j] = Graph::noEdge;
					} else {
						graph->weights[i][j] = stoi(numstr);
					}
				}
				
			}
			
			graph->setNodes(numberOfVertices);
			
			return graph;

		}
		
	} catch(const invalid_argument& ia) {
		throw runtime_error("Wrong input file format - invalid argument for the number of vertices");
	} catch(const out_of_range& oor) {
		throw runtime_error("Wrong input file format - number of vertices out of range");
	} catch ( ... ) {
		// everything else
		throw runtime_error("Something went wrong");
	}

}

void Graph::setNodes(int numberOfVertices) {
	int A = (int)'A';
	for(int i=0; i<numberOfVertices; ++i) {
		int node = A + i;
		this->nodes.push_back(string(1,(char)node));
	}
}

bool Graph::validateUnsignedInt(const int number) {
	return number > 0;
}