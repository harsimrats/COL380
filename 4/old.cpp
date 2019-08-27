#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string.h>
#include <limits>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <omp.h>

using namespace std;

int max_row, max_col;
int proc_id, num_proc;

struct element
{
	int i, r, c;
	float f;
};

void read_file(vector < vector < struct element > > *matrix, string input_file)
{
	FILE *f;
	f = fopen(input_file.c_str(), "rb");
	int last = -1;
	while(!feof(f))
	{
		struct element e;
		fread(&e.r, 4, 1, f);
		fread(&e.c, 4, 1, f);
		fread(&e.i, 4, 1, f);
		fread(&e.f, 4, 1, f);

		int i = e.r / (max_row / num_proc);
		(*matrix)[i].push_back(e);
		last = i;
	}
	if (last != -1)
		(*matrix)[last].pop_back();
	
	fclose(f);
}

bool sizecom_row1(struct element x, struct element y)
{
	return (x.r < y.r ||  (x.r == y.r && x.f < y.f) || (x.r == y.r && x.f == y.f && x.c < y.c));
}

bool sizecom_row2(struct element x, struct element y)
{
	return (x.r < y.r)	|| (x.r == y.r && x.c < y.c);
}

bool sizecom_col1(struct element x, struct element y)
{
	return (x.c < y.c ||  (x.c == y.c && x.i < y.i) || (x.c == y.c && x.i == y.i && x.r < y.r));
}

bool sizecom_col2(struct element x, struct element y)
{
	return (x.c < y.c)	|| (x.c == y.c && x.r < y.r);
}

void sort_matrix(vector < struct element > *matrix, int rank, int r)
{
	if (r == 1)
	{
		vector< struct element > matrix2(*matrix);
		sort((*matrix).begin(), (*matrix).end(), sizecom_row1);
		sort(matrix2.begin(), matrix2.end(), sizecom_row2);

		#pragma omp parallel for
		for(int i = 0; i < matrix2.size(); i++)
		{
			(*matrix)[i].c = matrix2[i].c;
		}

	}
	
	else
	{
		vector< struct element > matrix2(*matrix);
		sort((*matrix).begin(), (*matrix).end(), sizecom_col1);
		sort(matrix2.begin(), matrix2.end(), sizecom_col2);

		#pragma omp parallel for
		for(int i = 0; i < matrix2.size(); i++)
		{
			(*matrix)[i].r = matrix2[i].r;
		}
	}
}

void gen_out_file(vector < vector < struct element > > *matrix, string out_file)
{
	FILE *f = fopen(out_file.c_str(), "wb");
	for (int i = 0; i < (*matrix).size(); i++)
	{
		for (int j = 0; j < (*matrix)[i].size(); j++)
		{
			fwrite(&((*matrix)[i][j].r), 4, 1, f);
			fwrite(&((*matrix)[i][j].c), 4, 1, f);
			fwrite(&((*matrix)[i][j].i), 4, 1, f);
			fwrite(&((*matrix)[i][j].f), 4, 1, f);
		}	
	}

	fclose(f);
}


int main (int argc, char *argv[])
{
	if (argc != 5)
		return 0;

	string input_file = argv[1];
	string out_file = argv[2];
	max_row = stoi(argv[3]);
	max_col = stoi(argv[4]);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

	MPI_Datatype new_type, old_type[2] = {MPI_INT, MPI_FLOAT};
	int blockcounts[2] = {3, 1};
	MPI_Aint offsets[2] = {0, 12};

	MPI_Type_create_struct(2, blockcounts, offsets, old_type, &new_type);
	MPI_Type_commit(&new_type);

	if (proc_id == 0)
	{
		vector < vector < struct element > > matrix;
		matrix.resize(num_proc);

		read_file(&matrix, input_file);

		for (int sort_time = 0; sort_time < 4; sort_time++)
		{
			// SORT ROWS
			for(int i = 1; i < num_proc; i++)
			{
				int data = matrix[i].size();
				MPI_Send(&data, 1, MPI_INT, i, i, MPI_COMM_WORLD);
			}

			for(int i = 1; i < num_proc; i++)
			{
				MPI_Send(&matrix[i][0], matrix[i].size(), new_type, i, i, MPI_COMM_WORLD);
			}

			sort_matrix(&matrix[proc_id], proc_id, 1);
			for(int i = 1; i < num_proc; i++)
			{
				MPI_Recv(&matrix[i][0], matrix[i].size(), new_type, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			
			// SORT COLUMNS
			vector < vector < struct element > > col_matrix;
			col_matrix.resize(num_proc);
			for(int i = 0; i < matrix.size(); i++)
			{
				for(int j = 0; j < matrix[i].size(); j++)
				{
					int index = matrix[i][j].c / (max_col / num_proc);
					if(index >= num_proc)
						index = num_proc - 1;
					col_matrix[index].push_back(matrix[i][j]);
				}
			}
			for(int i = 1; i < num_proc; i++)
			{
				int data = col_matrix[i].size();
				MPI_Send(&data, 1, MPI_INT, i, i, MPI_COMM_WORLD);
			}

			for(int i = 1; i < num_proc; i++)
			{
				MPI_Send(&col_matrix[i][0], col_matrix[i].size(), new_type, i, i, MPI_COMM_WORLD);
			}

			sort_matrix(&col_matrix[proc_id], proc_id, 0);

			for(int i = 1; i < num_proc; i++)
			{
				MPI_Recv(&col_matrix[i][0], col_matrix[i].size(), new_type, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			vector< vector < struct element > > mm2;
			mm2.resize(num_proc);
			for(int i = 0; i < col_matrix.size(); i++)
			{	
				for(int j = 0; j < col_matrix[i].size(); j++)
				{
					mm2[col_matrix[i][j].r*num_proc/max_row].push_back(col_matrix[i][j]);
				}
			}

			matrix = mm2;
		}
		
		for (int i = 0; i < matrix.size(); i++)
		{
			sort(matrix[i].begin(), matrix[i].end(), sizecom_row2);
		}
		gen_out_file(&matrix, out_file);
	}
	
	else
	{
		for(int sort_time = 0; sort_time < 4; sort_time++)
		{
			MPI_Status stat;
			
			int size;
			MPI_Recv(&size, 1, MPI_INT, 0, proc_id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			vector < struct element > myrow;
			myrow.resize(size);

			MPI_Recv(&myrow[0], size, new_type, 0, proc_id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			sort_matrix(&myrow, proc_id, 1);

			MPI_Send(&myrow[0], size, new_type, 0, proc_id, MPI_COMM_WORLD);

			// COLUMN
			MPI_Recv(&size, 1, MPI_INT, 0, proc_id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			vector<struct element> mycol(size);
			MPI_Recv(&mycol[0], size, new_type, 0, proc_id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			sort_matrix(&mycol, proc_id, 0);

			MPI_Send(&mycol[0], size, new_type, 0, proc_id, MPI_COMM_WORLD);
		}
	}

	MPI_Finalize();

	return 0;
}