#include <stdio.h>
#include <string>
#include <cstring>
#include <stdlib.h>
#include "vector"
#include <algorithm>
#include <iostream>
#include <limits>
#include <math.h>
using namespace std;

typedef struct {
	int r;
	int c;
	int i;
	float f;
}node;

int main(int argc, char *argv[])
{
	FILE *fp;
	char *infile,*outfile;
	infile = argv[1];
	outfile = argv[2];
	int numprocs = 1;
	vector<vector<node>> matrix;
	matrix.resize(numprocs);
	int max_size = 0;
	fp = fopen(infile,"rb");
	while(1){
		node dummy;
		int r=0,c=0,i=0;
		float f=0;
		if(fread(&r,sizeof(int),1,fp)==0){
			// cout << "fuck at 1\n";
			break;
		}
		if(fread(&c,sizeof(int),1,fp)==0){
			cout << "fuck at 2\n";
			// break;
		}
		if(fread(&i,sizeof(int),1,fp)==0){
			// cout << "fuck at 3\n";
			break;
		}
		// cout << i << endl << flush;
		dummy.r = r;
		dummy.c = c;
		dummy.i = i;
		if(fread(&f,sizeof(float),1,fp) == 0){
			// cout << "fuck at 4\n";
			break;
		}
		// cout << f << endl << flush;
		dummy.f = f;
		matrix[0].push_back(dummy);
	}
	fclose(fp);
	fp = fopen(outfile,"w");
	for(int i=0;i<numprocs;i++){
		for(int j=0;j<matrix[i].size();j++){
			fprintf(fp, "%d %d %d %f\n",matrix[i][j].r,matrix[i][j].c,matrix[i][j].i,matrix[i][j].f );;
		}
	}	
	fclose(fp);
	return 0;
}