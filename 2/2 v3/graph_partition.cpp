#include <iostream>
#include <algorithm>
#include <random>
#include <set>
#include <cmath>
#include <climits>
#include <omp.h>

#include "graph_partition.h"

int NTHREADS = 16;

using namespace std;

int get_edge_cut(graph &G, vector<bool> &visited) {
	int edge_cut = 0;
	for(unsigned int i = 0; i < G.size(); i++) {
		int node = get<0>(G[i]);
		for(auto val: get<2>(G[i])) {
			if(visited[node] != visited[val.first])
				edge_cut += val.second;
		}
	}

	return edge_cut;
}

graph transform(graph &G, unordered_map<int, int> *map, unordered_map<int, int> *inv_map) {
	int index = 0;
	graph H;

	for(unsigned int i = 0; i < G.size(); i++) {
		(*map).insert(make_pair(get<0>(G[i]), index));
		(*inv_map).insert(make_pair(index, get<0>(G[i])));
		index++;
	}

	for(unsigned int i = 0; i < G.size(); i++) {
		int node = (*map)[get<0>(G[i])];
		int node_weight = get<1>(G[i]);
		unordered_map<int, int> neighbours;
		for(auto val: get<2>(G[i])) {
			neighbours.insert(make_pair((*map)[val.first], val.second));
		}

		H.push_back(make_tuple(node, node_weight, neighbours));
	}

	return H;
}

split rev_transform(vector<bool>& part, unordered_map<int, int> *inv_map) {
	split P_coarse;

	for(unsigned int i = 0; i < part.size(); i++) {
		P_coarse.insert(make_pair((*inv_map)[i], part[i]));
	}

	return P_coarse;
}

// void refining(graph &H, vector<bool> &P) {
// 	cout << "refining\n";
// 	int N = H.size();
// 	vector<int> gain_vector(N, 0);
// 	set<pair<int, int> > Q1, Q2;

// 	for(int i = 0; i < N; i++) {
// 		int node = get<0>(H[i]), gain = 0;
// 		for(auto neigh : get<2>(H[i])) {
// 			if(P[node] == P[neigh.first]) {
// 				gain += neigh.second;
// 			} else {
// 				gain -= neigh.second;
// 			}
// 		}

// 		gain_vector[i] = gain;
// 		if(P[node])
// 			Q1.insert(make_pair(gain, i));
// 		else
// 			Q2.insert(make_pair(gain, i));
// 	}

// 	int count[2];
// 	count[0] = Q1.size();
// 	count[1] = Q2.size();
// 	int iters = 1;
// 	while(abs(count[0] - count[1]) < 0.04 * N) {	
// 		cout << "Iters: " << iters++ << endl;
// 		set<pair<int, int> >::iterator it1 = Q1.begin(), it2 = Q2.begin();
// 		int u = -1;
// 		if((*it1).first > (*it2).first) {
// 			if((*it2).first < 0) {
// 				u = (*it2).second;
// 				Q2.erase(it2);
// 				count[1]--;
// 				count[0]++;
// 			}
// 		} else {
// 			if((*it1).first < 0) {
// 				u = (*it1).second;
// 				Q1.erase(it1);
// 				count[0]--;
// 				count[1]++;
// 			}
// 		}

// 		if(u != -1) {
// 			for(auto neigh: get<2>(H[u])) {
// 				int gain = gain_vector[neigh.first];
// 				pair<int, int> p = make_pair(gain, neigh.first);
// 				set<pair<int, int> >::iterator ft1 = Q1.find(p), ft2 = Q2.find(p);
// 				if(ft1 != Q1.end()) {
// 					if(P[u] == P[neigh.first])
// 						gain_vector[neigh.first] -= 2 * neigh.second;
// 					else
// 						gain_vector[neigh.first] += 2 * neigh.second; 
// 					Q1.erase(ft1);
// 					Q1.insert(make_pair(gain_vector[neigh.first], neigh.second));
// 				} else if(ft2 != Q2.end()) {
// 					if(P[u] == P[neigh.first])
// 						gain_vector[neigh.first] -= 2 * neigh.second;
// 					else
// 						gain_vector[neigh.first] += 2 * neigh.second; 
// 					Q2.erase(ft2);
// 					Q1.insert(make_pair(gain_vector[neigh.first], neigh.second));
// 				}
// 			}
// 			P[u] = !P[u];
// 		} else {
// 			break;
// 		}
// 	}
// }

vector<bool> partitioning(graph &G) {
	int N = G.size();
	vector<bool> partition;
	int min_edge_cut = INT_MAX, min_partition_weight;

	random_device random_device;
	mt19937 engine{random_device()};
	uniform_int_distribution<int> dist(0, N - 1);

	omp_set_num_threads(NTHREADS);

	#pragma omp parallel for
	for(int th = 0; th < NTHREADS; th++) {
		vector<bool> visited(N, false);
		vector<int> gain_vector(N, INT_MAX);

		int total_weight = 0, edge_cut = 0, partition_weight = 0;
		for(int i = 0; i < N; i++) {
			total_weight += get<1>(G[i]);
		}
 
		int s = dist(engine);
		set<pair<int, int> > Q;

		Q.insert(make_pair(0, s));

		while(partition_weight < total_weight / 2) {
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
			partition_weight += get<1>(G[u]);
			
			unordered_map<int, int> neighbours = get<2>(G[u]);
			for(auto val: neighbours) {
				if(!visited[val.first]) {
					set<pair<int, int> >::iterator it = Q.find(make_pair(gain_vector[val.first], val.first));
					if(it == Q.end()) {
						int gain = 0;
						for(auto next: get<2>(G[val.first])) {
							if(visited[next.first])
								gain -= next.second;
							else
								gain += next.second;
						}

						Q.insert(make_pair(gain, val.first));
						gain_vector[val.first] = gain;
					} else {
						gain_vector[val.first] -= 2 * val.second;
						Q.erase(it);
						Q.insert(make_pair(gain_vector[val.first], val.first));
					}	
				}
			}
		}

		//refining(G, visited);

		edge_cut = get_edge_cut(G, visited);

		#pragma omp critical
		{
			if(edge_cut < min_edge_cut) {
				min_edge_cut = edge_cut;
				partition = visited;
				min_partition_weight = partition_weight;
			}
		}
	}

	return partition;
}

split graph_divide(graph &G) {
	unordered_map<int, int> *map, *inv_map;
	map = new unordered_map<int, int>;
	inv_map = new unordered_map<int, int>;


	graph H = transform(G, map, inv_map);

	vector<bool> divide = partitioning(H);

	split P = rev_transform(divide, inv_map);

	return P;
}

pair<graph, graph> make_graph(graph &G, split P) {
	graph G1;
	graph G2;

	for(unsigned i = 0; i < G.size(); i++) {
		int node = get<0>(G[i]);
		int node_weight = get<1>(G[i]);
		unordered_map<int, int> neighbours;

		for(auto val: get<2>(G[i])) {
			if(!(P[val.first] ^ P[node])) {
				neighbours.insert(val);
			}
		}

		if(P[node])
			G1.push_back(make_tuple(node, node_weight, neighbours));
		else
			G2.push_back(make_tuple(node, node_weight, neighbours));
	}

	return make_pair(G1, G2);
}

unordered_map<int, int> graph_divide_K(graph &G, int k, int partition_num) {
	if(k == 2) {
		split P = graph_divide(G);
		unordered_map<int, int> partition;

		for(auto val: P) {
			if(val.second)
				partition.insert(make_pair(val.first, partition_num));
			else
				partition.insert(make_pair(val.first, partition_num - 1));
		}

		return partition;
	}

	split P = graph_divide(G);

	pair<graph, graph> do_graph = make_graph(G, P);

	unordered_map<int, int> partition1 = graph_divide_K(do_graph.first, k / 2, 2 * partition_num - 1);
	unordered_map<int, int> partition2 = graph_divide_K(do_graph.second, k / 2, 2 * partition_num + 1);

	unordered_map<int, int> partition;

	for(auto val: partition1) {
		partition.insert(val);
	}

	for(auto val: partition2) {
		partition.insert(val);
	}

	return partition;
}

vector<int> graph_partition(graph &G, int k) {
	unordered_map<int, int> partition = graph_divide_K(G, k, 1);

	vector<int> labels(partition.size());

	for(auto val:partition) {
		labels[val.first - 1] = val.second;
	}

	return labels;
}