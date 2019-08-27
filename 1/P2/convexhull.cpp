#include "iostream"
#include "vector"
#include "set"
#include "omp.h"

using namespace std;

int counter = 0, nthreads;
vector< pair<int, int> > hull;

int lineSide(pair<int, int> point1, pair<int, int> point2, pair<int, int> point)
{
	int temp = (point.second - point1.second) * (point2.first - point1.first) - (point2.second - point1.second) * (point.first - point1.first);
	if (temp > 0)
		return 1;
	if (temp < 0)
		return -1;
	return 0;
}

int distFromLine(pair<int, int> point1, pair<int, int> point2, pair<int, int> point)
{
	int temp = (point.second - point1.second) * (point2.first - point1.first) - (point2.second - point1.second) * (point.first - point1.first);
	if (temp>=0)
		return temp;
	return -1*temp;
}

void hullRecFunc(vector< vector<int> > image, pair<int, int> point1, pair<int, int> point2, int side)
{
	pair<int, int> ind;
	ind = make_pair(-1, -1);
	int max_dist = 0;

	int max_dist_i[image.size()];
	pair<int, int> ind_i[image.size()];
	#pragma omp parallel for
	for (int i = 0; i < image.size(); i++)
	{
		int var_dist = 0;
		pair<int,int> var_ind = make_pair(-1,-1);
		for (int j = 0; j < image[0].size(); j++)
		{
			if (image[i][j] == 1)
			{
				if (lineSide(point1, point2, make_pair(i,j)) == side)
				{
					int temp = distFromLine(point1, point2, make_pair(i,j));
					if (temp > var_dist)
					{
						var_ind = make_pair(i,j);
						var_dist = temp;
					}
				}
			}
		}
		max_dist_i[i] = var_dist;
		ind_i[i] = var_ind;
	}

	for (int i = 0; i < image.size(); i++)
	{
		if (max_dist < max_dist_i[i])
		{
			max_dist = max_dist_i[i];
			ind = ind_i[i];
		}
	}

	if (ind.first == -1 && ind.second == -1)
	{
		return;
	}
	#pragma omp critical
	{
		hull.push_back(ind);
	}

	if (counter < nthreads)
	{
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				#pragma omp atomic
					counter++;
				hullRecFunc(image, ind, point1, -lineSide(ind, point1, point2));
				#pragma omp atomic
					counter--;
			}
			#pragma omp section
			{
				#pragma omp atomic
					counter++;
				hullRecFunc(image, ind, point2, -lineSide(ind, point2, point1));
				#pragma omp atomic
					counter--;
			}
		}
	}
	else
	{
		hullRecFunc(image, ind, point1, -lineSide(ind, point1, point2));
		hullRecFunc(image, ind, point2, -lineSide(ind, point2, point1));
	}
}

vector< pair<int, int> > calcConvexHull(vector< vector<int> > image, int num_threads)
{
	omp_set_num_threads(num_threads);
	omp_set_nested(1);
	nthreads = num_threads;

	if (image.size() < 3)
	{
		cout << "Convex hull not possible\n";
		return hull;
	}

	pair<int, int> min_x, max_x;
	min_x = make_pair(image.size(), 0);
	max_x = make_pair(0, 0);
	for (int i = 0; i < image.size(); i++)
		for (int j = 0; j < image[0].size(); j++)
			if (image[i][j] == 1)
			{
				if (i < min_x.first)
				{
					min_x = make_pair(i,j);
				}
				if (i > max_x.first)
				{
					max_x = make_pair(i,j);
				}
			}

	hull.push_back(max_x);
	hull.push_back(min_x);
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			#pragma omp atomic
				counter++;
			hullRecFunc(image, min_x, max_x, 1);
			#pragma omp atomic
				counter--;
		}
		#pragma omp section
		{
			#pragma omp atomic
				counter++;
			hullRecFunc(image, min_x, max_x, -1);
			#pragma omp atomic
				counter--;
		}
	}

	set< pair<int, int> > s;
	for(int i = 0; i < hull.size(); i++ )
		s.insert( make_pair(hull[i].first, hull[i].second));
	hull.assign( s.begin(), s.end() );

	return hull;
}