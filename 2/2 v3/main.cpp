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
	for(int i = 1; i <= N; i++) {
		getline(cin, line);
		istringstream is(line);
		unordered_map<int, int> neighbours;
		while(is >> node) {
			neighbours.insert(make_pair(node, 1));
		}
		G.push_back(make_tuple(i, 1, neighbours));
	}

	// cout << "Graph read\n";

	vector<int> labels = graph_partition(G, k);

	vector<int> count(k, 0);

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
	// 	int node = get<0>(G[i]);
	// 	for(auto val: get<2>(G[i])) {
	// 		if(labels[node - 1] != labels[val.first - 1])
	// 			edge_cut += val.second;
	// 	}
	// }

	// cout << "Edge cut: " << edge_cut << endl;
	
	return 0;
}


