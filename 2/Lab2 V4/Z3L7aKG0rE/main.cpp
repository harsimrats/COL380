#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>

#include "graph_partition.h"

using namespace std;

int main(int argc, char *argv[]) {
	graph G;
	int N, E;
	string line;

	int k = atoi(argv[1]);

	getline(cin, line);
	istringstream is(line);
	is >> N >> E;
	
	int node;
	for(int i = 0; i < N; i++) {
		getline(cin, line);
		istringstream is(line);
		vector<int> neighbours;
		while(is >> node) {
			neighbours.push_back(node - 1);
		}
		G.push_back(neighbours);
	}

	// cout << "Graph read\n";

	vector<int> labels = graph_divide_K(G, k, 1);

	// vector<int> count(k, 0);

	for(auto entry : labels) {
		cout << entry << " ";
		// count[entry]++;
	}

	cout << endl;

	// for(auto ct: count) {
	// 	cout << ct << endl;
	// }

	// int edge_cut = 0;
	// for(unsigned int i = 0; i < G.size(); i++) {
	// 	for(auto val: G[i]) {
	// 		if(labels[i] != labels[val])
	// 			edge_cut += 1;
	// 	}
	// }

	// cout << "Edge cut: " << edge_cut << endl;
	
	return 0;
}


