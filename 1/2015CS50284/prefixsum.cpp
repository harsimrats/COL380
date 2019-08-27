#include "iostream"
#include "vector"
#include "omp.h"
#include "cmath"

using namespace std;

vector<int> calcPrefixSum ( vector<int> input, int num_threads)
{
	int SIZE = input.size();
	int i, j, STEPS, limit, temp, t;
	float upperSteps = log2(SIZE);
	STEPS = upperSteps;
	if (STEPS < upperSteps)
		upperSteps = STEPS + 1;

	temp = 1;
	for (j = 1; j <= STEPS; j++)
	{
		temp = temp*2;
		limit = SIZE/temp;
		#pragma omp parallel for num_threads(num_threads) private(i)
		for(i = 0; i < limit; i++)
		{
			input[i*temp-1+temp] += input[i*temp-1+(temp/2)];
		}
	}

	limit = 0;
	t = pow(2,upperSteps);
	for (j = 1; j < upperSteps; j++)
	{
		limit = limit*2 + 1;
		t = t/2;
		#pragma omp parallel for num_threads(num_threads) private(i)
		for(i = 0; i < limit; i++)
		{
			if((i+1)*t + t/2 - 1 < SIZE)
				input[(i+1)*t + t/2 - 1] += input[(i+1)*t - 1];
		}
	}
	return input;
}