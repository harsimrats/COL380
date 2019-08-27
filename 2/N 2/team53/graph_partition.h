#ifndef GP_H
#define GP_H

#include <vector>
#include <utility>
#include <tuple>
#include <unordered_map>

using namespace std;

typedef vector<tuple<int, int, unordered_map<int, int>* > > graph;
typedef unordered_map<int, int> parent;
typedef vector<pair<graph*, parent *> > coarse_seq;
typedef unordered_map<int, bool> split;

vector<int> graph_partition(graph *, int);

#endif