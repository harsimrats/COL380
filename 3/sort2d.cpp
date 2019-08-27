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

			if(fread(&num, 4, 1, f) == 0)
				break;

			if(feof(f))	
				break;

			fread(str, n, 1, f);
			str[n] = '\0';

			matrix[i-1].push_back(make_pair(num, str));
		}
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

int main(int argc, char const *argv[])
{
	int rank, size, chunk_size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0)
	{
		size = atoi(argv[1]);
		string basefile = argv[2];

		read_file(basefile);
		balance_matrix();

		chunk_size = max_len / size;
		extra_chunk = max_len % size;

		while(!check_sorted())
		{
			for (int k = 1; k < size-1 ; k++)
			{
				vector < float > v;
				for (int j = 0; j < chunk_size; j++)
				{
					v.push_back(matrix[j+rank*chunk_size][k].first);	
				}

				MPI_Send(&v[0], v.size(), MPI_INT, k+1, 1, MPI_COMM_WORLD);
				
				MPI_Recv(&A[j][0], A[j].size(), pairr, i, 1, MPI_COMM_WORLD, &status);
				
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

	MPI_Finalize();
	return 0;
}