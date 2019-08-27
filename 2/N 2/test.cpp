#include <iostream>
#include <vector>
#include <algorithm>    // std::random_shuffle
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <unordered_map>
#include <random>
using namespace std;


// bool has_key(int key, unordered_map<int, int> map) {
// 	unordered_map<int, int>::iterator itr = map.find(key);
// 	if(itr == map.end())
// 		return false;

// 	return true;
// }

// bool has_key(int key, unordered_map<int, float> map) {
// 	unordered_map<int, float>::iterator itr = map.find(key);
// 	if(itr == map.end())
// 		return false;

// 	return true;
// }


int main() {

	// srand(unsigned(time(0)));

	// vector<int> A;
	// A.push_back(1);
	// A.push_back(2);
	// A.push_back(3);
	// A.push_back(4);

	// for(auto val: A) {
	// 	cout << val << " ";
	// }

	// vector<int> *B = &A;

	// random_shuffle((*B).begin(), (*B).end());

	// A.

	// for(auto i: A) {
	// 	cout << i << endl;
	// }

	unordered_map<int, int> *p = new unordered_map<int, int>();

	for(int i = 0; i < 10; i++) {
		if(i%2 == 0)
			(*p).insert(make_pair(i, i));
	}

	auto ran = next(begin(*p)), rand(between(0, *p.size()));

	cout<<ran;

	// unordered_map<int, unordered_map<int, int>*> s;
	// s.insert(make_pair(1, p));
	// s.insert(make_pair(1, p));
	// s.insert(make_pair(1, p));
	// s.insert(make_pair(1, p));

	// unordered_map<int, int> *q = s[1];

	// (*q).insert(make_pair(11,11));

	// for(auto val: *(s[1])) {
	// 	cout << val.first << " " << val.second << endl;
	// }

	// std::random_device random_device;
	// std::mt19937 engine{random_device()};
	// std::uniform_int_distribution<int> dist(0, 100 - 1);

	// cout << dist(engine) << " " << dist(engine) << endl;


}