#include <fstream>
#include <iostream>
#include <string.h>
#include <limits>
#include <vector>
#include <algorithm>
#include <mpi.h>

using namespace std;
#define ipair pair < float, char* >

float inf = numeric_limits<float>::max();
vector < vector < ipair > > matrix;
int size, max_len;

void printm()
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < matrix[i].size(); j++)
		{
			cout << matrix[i][j].first << "," << matrix[i][j].second << "  " ;
		}
		cout << endl;
	}
}

void read_file(string basefile)
{
	matrix.reserve(size);

	for (int i = 1; i <= size; i++)
	{
		FILE *f;
		// string filename = basefile + boost::lexical_cast<string>(i);
		string filename = basefile + to_string(static_cast < long long > (i));
		
		char* f_arr;
		f_arr = (char*) malloc(sizeof(filename));
		strcpy(f_arr, filename.c_str());
		f = fopen(f_arr, "rb");

		int n;
		fread(&n, 4, 1, f);

		while(1)
		{
			float num;
			char* str;
			str = (char*) malloc(n);

			if(fread(&num, 4, 1, f) == 0){
				break;
			}
			if(feof(f))	break;
			fread(str, n, 1, f);
			str[n] = '\0';

			matrix[i-1].push_back(make_pair(num, str));
		}
		// matrix[i-1].pop_back();
		fclose(f);
	}
}

int get_max_col_size()
{
	int max_len = 0;
	for (int i = 0; i < size; i++)
	{
		if(max_len < matrix[i].size())
			max_len = matrix[i].size();
	}
	return max_len;
}

void balance_matrix()
{
	int max_len = get_max_col_size();

	char ch = 'a';
	for (int i = 0; i < size; i++)
	{
		int diff = max_len - matrix[i].size();

		if (diff > 0)
		{
			for (int j = 0; j < diff; j++)
			{
				matrix[i].push_back(make_pair(inf, &ch));
			}
		}
	}
}

bool check_sorted()
{
	for (int k = 0; k < max_len; k++)
	{
		vector < ipair > v;
		for (int j = 0; j < size; j++)
		{
			if (matrix[j].size() <= max_len)
			{
				v.push_back(matrix[j][k]);
			}
		}
		if(!is_sorted(v.begin(), v.end()))
			return false;
	}

	return true;
}

void sort_matrix()
{
	while(!check_sorted())
	{
		for (int k = 0; k < max_len ; k++)
		{
			vector < ipair > v;
			for (int j = 0; j < size; j++)
			{
				v.push_back(matrix[j][k]);	
			}
			sort(v.begin(), v.end());
			for (int j = 0; j < v.size(); j++)
			{
				matrix[j][k] = v[j];
			}
		}
		
		for (int i = 0; i < size; i++)
		{
			sort(matrix[i].begin(), matrix[i].end());
		}
	}
}

void gen_out_file(string basefile)
{
	FILE *f;
	string filename = basefile + to_string(static_cast < long long > (0));
	
	char* f_arr;
	f_arr = (char*) malloc(sizeof(filename));
	strcpy(f_arr, filename.c_str());
	f = fopen(f_arr, "wb");

	for (int j = 0; j < max_len; j++)
	{
		for (int i = 0; i < size; i++)
		{
			if (matrix[i][j].first != inf)
			{
				fwrite(&(matrix[i][j].first), 4, 1, f);
				fwrite((matrix[i][j].second), strlen(matrix[i][j].second), 1, f);
			}
		}	
	}

	fclose(f);
}

int main (int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int procid;
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);

	if(procid == 0)
	{
		size = atoi(argv[1]);
		string basefile = argv[2];

		read_file(basefile);

		balance_matrix();

		sort_matrix();

		gen_out_file(basefile);
	}

	MPI_Finalize();

	return 0;
}