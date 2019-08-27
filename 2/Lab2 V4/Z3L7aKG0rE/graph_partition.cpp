#include <iostream>
#include <algorithm>
#include <random>
#include <set>
#include <cmath>
#include <climits>
#include <omp.h>

#include "graph_partition.h"

int NTHREADS = 64;

using namespace std;

int get_edge_cut(graph &G, vector<bool> &visited) {
	int edge_cut = 0;
	for(unsigned int i = 0; i < G.size(); i++) {
		for(auto val: G[i]) {
			if(visited[i] != visited[val])
				edge_cut += 1;
		}
	}

	return edge_cut;
}

vector<bool> partitioning(graph &G) {
	int N = G.size();
	vector<bool> partition;
	int min_edge_cut = INT_MAX;

	random_device random_device;
	mt19937 engine{random_device()};
	uniform_int_distribution<int> dist(0, N - 1);

	omp_set_num_threads(NTHREADS);

	#pragma omp parallel for
	for(int th = 0; th < NTHREADS; th++) {
		vector<bool> visited(N, false);
		vector<int> gain_vector(N, INT_MAX);

		int edge_cut = 0, partition_weight = 0;
		
 
		int s = dist(engine);
		set<pair<int, int> > Q;

		Q.insert(make_pair(0, s));

		while(partition_weight < N / 2) {
			int u;
			if(Q.size() == 0) {
				for(int i = 0; i < N; i++) {
					if(!visited[i]) {
						u = i;
						break;
					}
				}
			} else {
				set<pair<int, int> >::iterator it = Q.begin();
				u = (*it).second;
				Q.erase(it);
			}

			visited[u] = true;
			partition_weight += 1;
			
			vector<int> neighbours = G[u];
			for(auto val: neighbours) {
				if(!visited[val]) {
					set<pair<int, int> >::iterator it = Q.find(make_pair(gain_vector[val], val));
					if(it == Q.end()) {
						int gain = 0;
						for(auto next: G[val]) {
							if(visited[next])
								gain -= 1;
							else
								gain += 1;
						}

						Q.insert(make_pair(gain, val));
						gain_vector[val] = gain;
					} else {
						gain_vector[val] -= 2;
						Q.erase(it);
						Q.insert(make_pair(gain_vector[val], val));
					}	
				}
			}
		}

		edge_cut = get_edge_cut(G, visited);

		#pragma omp critical
		{
			if(edge_cut < min_edge_cut) {
				min_edge_cut = edge_cut;
				partition = visited;
			}
		}
	}

	return partition;
}

pair<graph, graph> make_graph(graph &G, vector<bool> &divide, vector<int> &P1, vector<int> &P2) {
	graph G1;
	graph G2;

	vector<int> inv_p(divide.size());

	int count[2] = {0, 0};
	for(unsigned i = 0; i < divide.size(); i++) {
		if(divide[i]) {
			P1.push_back(i);
			inv_p[i] = count[0]++;
		}
		else {
			P2.push_back(i);
			inv_p[i] = count[1]++;
		}
	}

	for(unsigned i = 0; i < G.size(); i++) {
		vector<int> neighbours;

		for(auto val: G[i]) {
			if(divide[val] == divide[i]) {
				neighbours.push_back(inv_p[val]);
			}
		}

		if(divide[i])
			G1.push_back(neighbours);
		else
			G2.push_back(neighbours);
	}

	return make_pair(G1, G2);
}

vector<int> graph_divide_K(graph &G, int k, int partition_num) {
	if(k == 2) {
		vector<bool> divide = partitioning(G);
		vector<int> partition;

		for(auto val: divide) {
			if(val)
				partition.push_back(partition_num);
			else
				partition.push_back(partition_num - 1);
		}

		return partition;
	}

	vector<bool> divide = partitioning(G);

	vector<int> P1, P2;

	pair<graph, graph> do_graph = make_graph(G, divide, P1, P2);

	vector<int> partition1, partition2;

	partition1 = graph_divide_K(do_graph.first, k / 2, 2 * partition_num - 1);
	partition2 = graph_divide_K(do_graph.second, k / 2, 2 * partition_num + 1);	


	vector<int> partition(partition1.size() + partition2.size());

	for(unsigned int i = 0; i < partition1.size(); i++)
		partition[P1[i]] = partition1[i];

	for(unsigned int i = 0; i < partition2.size(); i++)
		partition[P2[i]] = partition2[i];

	return partition;
}