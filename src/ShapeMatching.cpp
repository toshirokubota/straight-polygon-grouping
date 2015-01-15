
#include <ShapeMatching.h>
#include <szmexutilitytemplate.h>
#include <mex.h>

/*
Use dynamic programming to compute the matching cost.
Assume that path2 is not shorter than path1.
*/
float pathCost(const vector<CParticleF>& path1, const vector<CParticleF>& path2)
{
	vector<float> row(path2.size(), 0);
	vector<vector<float>> D(path1.size(), row);
	vector<vector<float>> costs(path1.size(), row);
	for (int i = 0; i<path1.size(); ++i)
	{
		for (int j = 0; j<path2.size(); ++j)
		{
			D[i][j] = Distance(path1[i], path2[j]);
		}
	}
	costs[0][0] = D[0][0];
	for (int i = 1; i<path1.size(); ++i)
	{
		costs[i][0] = costs[i - 1][0] + D[i][0];
	}
	for (int i = 1; i<path2.size(); ++i)
	{
		costs[0][i] = costs[0][i - 1] + D[0][i];
	}
	for (int i = 1; i<path1.size(); ++i)
	{
		for (int j = 1; j<path2.size(); ++j)
		{
			float val1 = costs[i - 1][j - 1] + D[i][j];
			float val2 = costs[i][j - 1] + D[i][j];
			float val3 = costs[i - 1][j] + D[i][j];
			costs[i][j] = Min(val1, Min(val2, val3));
		}
	}
	/*printf("Distance:\n");
	for (int i = 0; i < path1.size(); ++i)
	{
		for (int j = 0; j < path2.size(); ++j)
		{
			printf("%f ", D[i][j]);
		}
		printf("\n");
	}
	printf("Cost:\n");
	for (int i = 0; i < path1.size(); ++i)
	{
		for (int j = 0; j < path2.size(); ++j)
		{
			printf("%f ", costs[i][j]);
		}
		printf("\n");
	}*/
	return costs[path1.size() - 1][path2.size() - 1];
}

