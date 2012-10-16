


double bssa_num_pairs(int sizeI, int sizeJ, bool sameBin)
{
	if(sameBin) return sizeI*.5*(sizeI-1);
	else return sizeI*(double)sizeJ;
}

//NOTE: length of binLens is numBins,
//      length of kmax is numBins*numBins,
//      length of sums array should be numBins+1
void bssa_refresh_sums(int numBins, const int* binLens,
		const double* kmax, double* sums)
{
	double netSum = 0.0;
	for(int i = 0; i < numBins; i++) {
		double rowSum = 0.0;
		int lenI = binLens[i];
		for(int j = 0; j <= i; j++) {
			double kernel = kmax[i*binLens+j]
			rowSum += kernel*bssa_num_pairs(lenI, binLens[j], i == j);
		}
		sums[i] = rowSum;
		netSum += rowSum;
	}
	sums[numBins] = netSum;
}

//call after adding particle, other bin counts should be unchanged
void bssa_on_part_add(int bin, int numBins, const int* binLens, double* sums)
{
	double* kmaxRow = kmax + bin*binLens;
	double netDiff = 0.0;
	for(int i = 0; i < bin; i++) {
		netDiff += binLens[i]*kmaxRow[i];
	}
	int selfPairsDiff = binLens[bin] - 1;
	netDiff += selfPairsDiff*kmaxRow[bin];
	sums[bin] += netDiff;
	for(int j = bin+1; j < numBins; j++) {
		double diff = binLens[j]*kmaxRow[j];
		sums[j] += diff;
		netDiff += diff;
	}
	sums[numBins] += netDiff;
}

//call after removing particle, other bin counts should be unchanged
void bssa_on_part_add(int bin, int numBins, const int* binLens, double* sums)
{

}

//TODO: possibly always keep track of row sums so that it does not
//      need to be re-computed every timestep...
void bssa_coags_for_timestep(int numBins, int* binLens, const double* kmax,
		double duration)
{
	
}