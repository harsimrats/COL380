#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <omp.h>
#include <algorithm>
#include <mpi.h>

using namespace std;

#define COUNT 4

int max_rows, max_cols;
int world_size, world_rank;
int row_start, row_end;
int col_start, col_end;

typedef struct node {
	int row, col;
	int key_int;
	float key_flt;
} node;

struct row_comp {
	bool operator() (node a, node b) { 
		if(a.key_flt < b.key_flt)
			return true;
		if(a.key_flt == b.key_flt)
			return a.col < b.col;

		return false;
	}
} row_comp_object;

struct col_comp {
	bool operator() (node a, node b) { 
		if(a.key_int < b.key_int)
			return true;
		if(a.key_int == b.key_int)
			return a.row < b.row;

		return false;
	}
} col_comp_object;

struct ind_comp {
	bool operator() (node a, node b) { 
		if(a.row < b.row)
			return true;
		if(a.row == b.row)
			return a.col < b.col;

		return false;
	}
} ind_comp_object;

void read_data(char *filename, vector<vector<node> > &data);
void print(vector<vector<node> > &data);
void sort_rows(vector<node> &rows);
void make_column_data(vector<vector<node> > &row_data, vector<vector<node> > &col_data);
void sort_cols(vector<node> &cols);
void make_row_data(vector<vector<node> > &col_data, vector<vector<node> > &row_data);
void sort_indices(vector<vector<node> >&row_data);
void write_data(char *filename, vector<vector<node> > &data);

int main(int argc, char* argv[]) {
	// Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	if(argc != 5) {
		cout << "sort2d <inputfile> <outfile> <maxrows> <maxcolumns>\n";
		exit(1);
	}
	
	char* infile = argv[1];
	char* outfile = argv[2];

	max_rows = atoi(argv[3]);
	max_cols = atoi(argv[4]);

	row_start = (max_rows / world_size) * world_rank;
	row_end = (max_rows / world_size) * (world_rank + 1);

	col_start = (max_cols / world_size) * world_rank;
	col_end = (max_cols / world_size) * (world_rank + 1);

	MPI_Datatype node_mpi, oldtypes[2];
	int blockcounts[2];
	MPI_Aint offsets[2];

	offsets[0] = 0; oldtypes[0] = MPI_INT; blockcounts[0] = 3;
	offsets[1] = 3 * 4; oldtypes[1] = MPI_FLOAT; blockcounts[1] = 1;

	MPI_Type_create_struct(2, blockcounts, offsets, oldtypes, &node_mpi);
	MPI_Type_commit(&node_mpi);


	if(world_rank == world_size - 1) {
		row_end = max_rows;
		col_end = max_cols;
	}

	if(world_rank == 0) {
		vector<vector<node> > row_data(world_size);
		read_data(infile, row_data);
		
		for(int count = 0; count < COUNT; count++) {
			MPI_Request reqs[world_size - 1];
			MPI_Status stats[world_size - 1];
			for(int i = 1; i < world_size; i++) {
				MPI_Isend(&row_data[i][0], row_data[i].size(), node_mpi, i, i, MPI_COMM_WORLD, &reqs[i - 1]);
			}

			sort_rows(row_data[0]);
			MPI_Waitall(world_size - 1, reqs, stats);

			for(int i = 1; i < world_size; i++) {
				MPI_Irecv(&row_data[i][0], row_data[i].size(), node_mpi, i, i, MPI_COMM_WORLD, &reqs[i - 1]);
			}
			MPI_Waitall(world_size - 1, reqs, stats);

			// print(row_data);

			vector<vector<node> > col_data(world_size);
			make_column_data(row_data, col_data);

			for(int i = 1; i < world_size; i++) {
				MPI_Isend(&col_data[i][0], col_data[i].size(), node_mpi, i, i, MPI_COMM_WORLD, &reqs[i - 1]);
			}

			sort_cols(col_data[0]);
			MPI_Waitall(world_size - 1, reqs, stats);

			for(int i = 1; i < world_size; i++) {
				MPI_Irecv(&col_data[i][0], col_data[i].size(), node_mpi, i, i, MPI_COMM_WORLD, &reqs[i - 1]);
			}
			MPI_Waitall(world_size - 1, reqs, stats);

			make_row_data(col_data, row_data);

			// print(row_data);
		}

		sort_indices(row_data);
		// print(row_data);
		write_data(outfile, row_data);
		
	} else {
		for(int count = 0; count < COUNT; count++) {
			MPI_Status stat;
			MPI_Probe(0, world_rank, MPI_COMM_WORLD, &stat);
			int size;
			MPI_Get_count(&stat, node_mpi, &size);

			vector<node> rows(size);
			MPI_Recv(&rows[0], size, node_mpi, 0, world_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			sort_rows(rows);
			
			MPI_Send(&rows[0], size, node_mpi, 0, world_rank, MPI_COMM_WORLD);
			
			// for(int i = 0; i < size; i++) {
			// 	cout << rows[i].row << " " << rows[i].col << " " << rows[i].key_int << " " << rows[i].key_flt << endl; 
			// }

			// cout << "********************\n";

			MPI_Probe(0, world_rank, MPI_COMM_WORLD, &stat);
			MPI_Get_count(&stat, node_mpi, &size);

			vector<node> cols(size);
			MPI_Recv(&cols[0], size, node_mpi, 0, world_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			sort_cols(cols);

			MPI_Send(&cols[0], size, node_mpi, 0, world_rank, MPI_COMM_WORLD);

			// for(int i = 0; i < size; i++) {
			// 	cout << cols[i].row << " " << cols[i].col << " " << cols[i].key_int << " " << cols[i].key_flt << endl; 
			// }

			// cout << "********************\n";
		}

	}

	// Finalize the MPI environment.
	MPI_Finalize();
	return 0;
}

void read_data(char *filename, vector<vector<node> > &data) {
	FILE* file = fopen(filename, "rb");
	
	int row, col, key_int;
	float key_flt;
	while(fread(&row, 4, 1, file) == 1) {
		fread(&col, 4, 1, file);
		fread(&key_int, 4, 1, file);
		fread(&key_flt, 4, 1, file);

		node val;
		val.row = row; val.col = col;
		val.key_int = key_int; val.key_flt = key_flt;

		// cout << row << " " << col << " " << key_int << " " << key_flt << endl;
		int index = val.row / (max_rows / world_size);
		if(index >= world_size)
			index = world_size - 1;

		data[index].push_back(val);
	}

	fclose(file);
}

void print(vector<vector<node> > &data) {
	for(int i = 0; i < data.size(); i++) {
		for(int j = 0; j < data[i].size(); j++) {
			cout << data[i][j].row << " " << data[i][j].col << " " << data[i][j].key_int << " " << data[i][j].key_flt << endl;
		}
		cout << "***************************\n";
	}
}

void sort_rows(vector<node> &rows) {
	vector<vector<node> > nodes(row_end - row_start);
	vector<vector<int> > cols(row_end - row_start);

	for(int i = 0; i < rows.size(); i++) {
		nodes[rows[i].row - row_start].push_back(rows[i]);
		cols[rows[i].row - row_start].push_back(rows[i].col);
	}

	omp_set_num_threads(omp_get_num_procs());
	#pragma omp parallel for
	for(int i = 0; i < row_end - row_start; i++) {
		sort(nodes[i].begin(), nodes[i].end(), row_comp_object);
		sort(cols[i].begin(), cols[i].end());

		for(int j = 0; j < nodes[i].size(); j++) {
			nodes[i][j].col = cols[i][j];
		}
	}

	int index = 0;
	for(int i = 0; i < nodes.size(); i++) {
		for(int j = 0; j < nodes[i].size(); j++) {
			rows[index++] = nodes[i][j];
		}
	}

	return;
}

void sort_cols(vector<node> &cols) {
	vector<vector<node> > nodes(col_end - col_start);
	vector<vector<int> > rows(col_end - col_start);

	for(int i = 0; i < cols.size(); i++) {
		nodes[cols[i].col - col_start].push_back(cols[i]);
		rows[cols[i].col - col_start].push_back(cols[i].row);
	}

	omp_set_num_threads(omp_get_num_procs());
	#pragma omp parallel for
	for(int i = 0; i < col_end - col_start; i++) {
		sort(nodes[i].begin(), nodes[i].end(), col_comp_object);
		sort(rows[i].begin(), rows[i].end());

		for(int j = 0; j < nodes[i].size(); j++) {
			nodes[i][j].row = rows[i][j];
		}
	}

	int index = 0;
	for(int i = 0; i < nodes.size(); i++) {
		for(int j = 0; j < nodes[i].size(); j++) {
			cols[index++] = nodes[i][j];
		}
	}

	return;
}

void make_column_data(vector<vector<node> > &row_data, vector<vector<node> > &col_data) {
	for(int i = 0; i < row_data.size(); i++) {
		for(int j = 0; j < row_data[i].size(); j++) {
			int index = row_data[i][j].col / (max_cols / world_size);
			if(index >= world_size)
			index = world_size - 1;

			col_data[index].push_back(row_data[i][j]);
		}
	}
	return;
}


//row_data already filled, needs to over-written
void make_row_data(vector<vector<node> > &col_data, vector<vector<node> > &row_data) {
	vector<int> indices(row_data.size(), 0);
	for(int i = 0; i < col_data.size(); i++) {
		for(int j = 0; j < col_data[i].size(); j++) {
			int index = col_data[i][j].row / (max_rows / world_size);
			if(index >= world_size)
			index = world_size - 1;

			row_data[index][indices[index]++] = col_data[i][j];
		}
	}
	return;
}

void sort_indices(vector<vector<node> >&row_data) {
	#pragma omp parallel for
	for(int i = 0; i < row_data.size(); i++) {
		sort(row_data[i].begin(), row_data[i].end(), ind_comp_object);
	}

	return;
}

void write_data(char *filename, vector<vector<node> > &data) {
	FILE* file = fopen(filename, "wb");
	
	for(int i = 0; i < data.size(); i++) {
		for(int j = 0; j < data[i].size(); j++) {
			fwrite(&(data[i][j].row), 4, 1, file);
			fwrite(&(data[i][j].col), 4, 1, file);
			fwrite(&(data[i][j].key_int), 4, 1, file);
			fwrite(&(data[i][j].key_flt), 4, 1, file);

		}
	}
	fclose(file);
}