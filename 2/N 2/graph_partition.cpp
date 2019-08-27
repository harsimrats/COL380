#include <iostream>
#include <algorithm>
#include <random>
#include <set>
#include <climits>

#include "graph_partition.h"

int NTHREADS = 16;

using namespace std;

bool has_key(int key, unordered_map<int, int> *map) {
	unordered_map<int, int>::iterator itr = (*map).find(key);
	if(itr == (*map).end())
		return false;

	return true;
}

bool has_key(int key, unordered_map<int, pair<int, unordered_map<int, int>* > > *map) {
	unordered_map<int, pair<int, unordered_map<int, int>* > >::iterator itr = (*map).find(key);
	if(itr == (*map).end())
		return false;

	return true;
}

coarse_seq coarsening(graph *G) {
	coarse_seq S;

	parent *root = new parent();
	S.push_back(make_pair(G, root));

	graph* G_parent = G;

	random_device rd;
	auto rng = default_random_engine(rd());

	int iters = 0;
	while(true) {
		graph *G_coarse = new graph();
		parent *p = new parent();

		shuffle(begin(*G_parent), end(*G_parent), rng);
		
		for(unsigned i = 0; i < (*G_parent).size(); i++) {
			int node = get<0>((*G_parent)[i]);
			unordered_map<int, int> *neighbours = get<2>((*G_parent)[i]);
			
			if(!has_key(node, p)) {
				int max_weight = 0, max_node = -1;
				for(auto val: *neighbours) {
					int neigh = val.first;
					int weight = val.second;
					if(!has_key(neigh, p)) {
						if(max_weight < weight) {
							max_weight = weight;
							max_node = neigh;
						}
					}
				}
			
				if(max_node != -1) {
					(*p).insert(make_pair(node, node));
					(*p).insert(make_pair(max_node, node));
				}

			}
		}

		int total_node_weight = 0, max_node_weight = 0;

		// Map to provide O(1) access to matched nodes edges
		unordered_map<int, pair<int, unordered_map<int, int>* > > matched_neighbors;

		for(unsigned i = 0; i < (*G_parent).size(); i++) {
			int node = get<0>((*G_parent)[i]);
			int node_weight = get<1>((*G_parent)[i]);
			unordered_map<int, int> *neighbours = get<2>((*G_parent)[i]);
			unordered_map<int, int> *new_neighbours;

			total_node_weight += node_weight;

			if(has_key(node, p))
				node = (*p)[node];

			if(has_key(node, &matched_neighbors)) {
				new_neighbours = matched_neighbors[node].second;
				node_weight += matched_neighbors[node].first;
			} else {
				new_neighbours = new unordered_map<int, int>();
			}

			if(node_weight > max_node_weight) {
				max_node_weight = node_weight;
			}

			for(auto val: *neighbours) {
				int neigh = val.first;
				int weight = val.second;

				if(has_key(neigh, p))
					neigh = (*p)[neigh];

				if(node != neigh) {
					if(has_key(neigh, new_neighbours))
						(*new_neighbours)[neigh] += weight;
					else
						(*new_neighbours).insert(make_pair(neigh, weight));
				}
			}

			if(has_key(node, p)) {
				if(!has_key(node, &matched_neighbors))
					matched_neighbors.insert(make_pair(node, make_pair(node_weight, new_neighbours)));
				else
					(*G_coarse).push_back(make_tuple(node, node_weight, new_neighbours));
			} else {
				(*G_coarse).push_back(make_tuple(node, node_weight, new_neighbours));
			}
		}

		S.push_back(make_pair(G_coarse, p));
		G_parent = G_coarse;

		iters++;

		if(max_node_weight > 0.05 * total_node_weight)
			break;
	}

	return S;
}

graph transform(graph *G, unordered_map<int, int> *map, unordered_map<int, int> *inv_map) {
	int index = 0;
	graph H;

	for(unsigned int i = 0; i < (*G).size(); i++) {
		(*map).insert(make_pair(get<0>((*G)[i]), index));
		(*inv_map).insert(make_pair(index, get<0>((*G)[i])));
		index++;
	}

	for(unsigned int i = 0; i < (*G).size(); i++) {
		int node = (*map)[get<0>((*G)[i])];
		int node_weight = get<1>((*G)[i]);
		unordered_map<int, int> *neighbours = new unordered_map<int, int>();
		for(auto val: *(get<2>((*G)[i]))) {
			(*neighbours).insert(make_pair((*map)[val.first], val.second));
		}

		H.push_back(make_tuple(node, node_weight, neighbours));
	}

	return H;
}

vector<bool> partitioning(const graph &G, int *min_edge_cut) {
	int N = G.size();
	vector<bool> *partition = NULL;
	*min_edge_cut = INT_MAX;

	random_device random_device;
	mt19937 engine{random_device()};
	uniform_int_distribution<int> dist(0, N - 1);

	for(int th = 0; th < NTHREADS; th++) {
		vector<bool> *visited = new vector<bool>(N, false);
		vector<int> node_weight(N);
		vector<int> gain;

		int total_weight = 0;
		for(int i = 0; i < N; i++) {
			int g = 0;
			for(auto val: *(get<2>(G[i]))) {
				g -= val.second;
			}
			gain.push_back(g);
		
			node_weight[i] = get<1>(G[i]);
			total_weight += get<1>(G[i]);
		}

		int s = dist(engine);

		int partition_weight = 0;

		set<pair<int, int> > Q;
		Q.insert(make_pair(gain[s], s));
		while(partition_weight < total_weight / 2) {
			pair<int, int> max = *(Q.rbegin());
			int u = max.second;
			
			(*visited)[u] = true;
			partition_weight += node_weight[u];

			Q.erase(max);

			for(auto next: *(get<2>(G[u]))) {
				int v = next.first, weight = next.second;
				if(!(*visited)[v]) {
					if(Q.find(make_pair(gain[v], v)) != Q.end())
						Q.erase(make_pair(gain[v], v));

					gain[v] += weight;
					Q.insert(make_pair(gain[v], v));
				} else {
					gain[v] += weight;
				}
			}
		}

		int edge_cut = 0;
		for(unsigned int i = 0; i < N; i++) {
			if((*visited)[i])
				edge_cut -= gain[i];
		}

		if(edge_cut < *min_edge_cut) {
			*min_edge_cut = edge_cut;
			partition = visited;
		}
	}

	return *partition;
}

split rev_transform(const vector<bool>& part, unordered_map<int, int> *inv_map) {
	split P_coarse;

	for(unsigned int i = 0; i < part.size(); i++) {
		P_coarse.insert(make_pair((*inv_map)[i], part[i]));
	}

	return P_coarse;
}

split uncoarsening(split& P_coarse, const coarse_seq& S) {
	split *P_coarser = &P_coarse, *P_fine;
	for(unsigned int i = S.size() - 1; i > 0; i--) {
		P_fine = new unordered_map<int, bool>();
		graph *old_graph = S[i - 1].first;
		parent *p = S[i].second;
		for(auto entry: *old_graph) {
			int node = get<0>(entry);
			int key = node;
			if(has_key(node, p)) {
				key = (*p)[node];
			}
			
			(*P_fine).insert(make_pair(node, (*P_coarser)[key]));
		}

		P_coarser = P_fine;
	}
	return *P_fine;
}

void delete_graph(graph *G) {
	for(unsigned i = 0; i < (*G).size(); i++) {
		delete get<2>((*G)[i]);
	}

	delete G;
}

void delete_coarsening_seq(coarse_seq *S) {
	for(int i = 1; i < (*S).size(); i++) {
		delete_graph((*S)[i].first);

		delete (*S)[i].second;
	}

	return;
}

split graph_divide(graph *G) {
	coarse_seq S = coarsening(G);
	
	graph *G_coarse = (S.back()).first;

	unordered_map<int, int> *map, *inv_map;
	map = new unordered_map<int, int>;
	inv_map = new unordered_map<int, int>;

	graph H = transform(G_coarse, map, inv_map);

	int *edge_cut = new int();
	vector<bool> part_H = partitioning(H, edge_cut);

	split P_coarse = rev_transform(part_H, inv_map);

	split P_fine = uncoarsening(P_coarse, S);

	delete_coarsening_seq(&S);

	return P_fine;
}

pair<graph*, graph*> make_graph(graph *G, split P) {
	graph *G1 = new graph();
	graph *G2 = new graph();

	for(unsigned i = 0; i < (*G).size(); i++) {
		int node = get<0>((*G)[i]);
		int node_weight = get<1>((*G)[i]);
		unordered_map<int, int> *neighbours = new unordered_map<int, int>();

		for(auto val: *(get<2>((*G)[i]))) {
			if(!(P[val.first] ^ P[node])) {
				(*neighbours).insert(val);
			}
		}

		if(P[node])
			(*G1).push_back(make_tuple(node, node_weight, neighbours));
		else
			(*G2).push_back(make_tuple(node, node_weight, neighbours));
	}

	return make_pair(G1, G2);
}

unordered_map<int, int>* graph_divide_K(graph *G, int k, int partition_num) {
	if(k == 2) {
		split P = graph_divide(G);
		unordered_map<int, int> *partition = new unordered_map<int, int>();

		for(auto val: P) {
			if(val.second)
				(*partition).insert(make_pair(val.first, partition_num));
			else
				(*partition).insert(make_pair(val.first, partition_num - 1));
		}

		return partition;
	}

	split P = graph_divide(G);

	pair<graph*, graph*> do_graph = make_graph(G, P);

	unordered_map<int, int>* partition1 = graph_divide_K(do_graph.first, k / 2, 2 * partition_num - 1);
	unordered_map<int, int>* partition2 = graph_divide_K(do_graph.second, k / 2, 2 * partition_num + 1);

	unordered_map<int, int>* partition = new unordered_map<int, int>();

	for(auto val: *partition1) {
		(*partition).insert(val);
	}

	for(auto val: *partition2) {
		(*partition).insert(val);
	}

	return partition;
}

vector<int> graph_partition(graph *G, int k) {
	unordered_map<int, int> *partition = graph_divide_K(G, k, 1);

	vector<int> labels((*partition).size());

	for(auto val: *partition) {
		labels[val.first - 1] = val.second;
	}

	return labels;
}