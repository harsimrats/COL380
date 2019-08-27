#include <stdio.h>
#include <stdlib.h>

int prob()
{
	int temp = rand();
	if (temp>0)
		return temp%100;
	else
		return (-1*temp)%100;
}

int main()
{
	int i, j;
	FILE *fp;
	fp = fopen("inp.txt", "w");
	fprintf(fp, "p1\n#cds\n4 3000 3000\n");
	for(i=0;i<4000;i++)
	{
		for (j = 0; j < 4000; j++)
		{
			if (prob()<30)
			{
				fprintf(fp, "1 ");
			}
			else
			{
				fprintf(fp, "0 ");
			}
		}
		fprintf(fp,"\n");
	}	
	return 0;
}