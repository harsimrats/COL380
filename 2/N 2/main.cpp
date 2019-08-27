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
		unordered_map<int, int> *neighbours = new unordered_map<int, int>();
		while(is >> node) {
			(*neighbours).insert(make_pair(node, 1));
		}
		G.push_back(make_tuple(i, 1, neighbours));
	}

	vector<int> labels = graph_partition(&G, k);

	for(auto entry : labels) {
		cout << entry << " ";
	}
	cout << endl;
	
	return 0;
}


