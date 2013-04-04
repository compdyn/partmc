#include <math.h>

//return kernel value between particles at binI-indexI and binJ-indexJ
double bssa_eval_kernel(int binI, int indexI, int binJ, int indexJ);

//coagulate particles at binI-indexI and binJ-indexJ,
//return bin number of coagulated particle
int bssa_coag(int binI, int indexI, int binJ, int indexJ);

//uniform between 0 and 1
double bssa_rand_double();

//uniform from {0, 1, ..., range-1}
int bssa_rand_int(int range);

double bssa_rand_exp(double mean)
{
	return -log(bssa_rand_double())/mean;
}

double bssa_num_pairs(int sizeI, int sizeJ, bool sameBin)
{
	if(sameBin) return sizeI*.5*(sizeI-1);
	else return sizeI*(double)sizeJ;
}

//NOTE: length of binLens is numBins,
//      length of kmax is numBins*numBins,
//      length of sums array should be numBins+1
void bssa_refresh_sums(int numBins, const double* kmax,
		const int* binLens, double* sums)
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

void bssa_on_part_add(int bin, int numBins, const double* kmax,
		int* binLens, double* sums)
{
	double* kmaxRow = kmax + bin*numBins;
	double netDiff = 0.0;
	for(int j = 0; j < bin; j++) {
		netDiff += binLens[j]*kmaxRow[j];
	}
	binLens[bin]++;
	int selfPairsDiff = binLens[bin] - 1;
	netDiff += selfPairsDiff*kmaxRow[bin];
	sums[bin] += netDiff;
	for(int i = bin+1; i < numBins; i++) {
		double diff = binLens[i]*kmaxRow[i];
		sums[i] += diff;
		netDiff += diff;
	}
	sums[numBins] += netDiff;
}

void bssa_on_part_remove(int bin, int numBins, const double* kmax,
		int* binLens, double* sums)
{
	double* kmaxRow = kmax + bin*numBins;
	double netDiff = 0.0;
	for(int j = 0; j < bin; j++) {
		netDiff -= binLens[j]*kmaxRow[j];
	}
	binLens[bin]--;
	int selfPairsDiff = -binLens[bin];
	netDiff += selfPairsDiff*kmaxRow[bin];
	sums[bin] += netDiff;
	for(int i = bin+1; i < numBins; i++) {
		double diff = -binLens[i]*kmaxRow[i];
		sums[i] += diff;
		netDiff += diff;
	}
	sums[numBins] += netDiff;
}

void bssa_try_coag(int numBins, const double* kmax, int* binLens,
		double* sums)
{
	double randVal = bssa_rand_double()*sums[numBins];
	int i, j;
	for(i = numBins - 1; i > 0; i--) {
		if(randVal < sums[i]) break;
		randVal -= sums[i];
	}
	double* kmaxRow = kmax + i*numBins;
	int binLenI = binLens[i];
	for(j = 0; j < i; j++) {
		randVal -= kmaxRow[j]*bssa_num_pairs(binLenI, binLens[j], 0);
		if(randVal < 0.0) break;
	}
	
	int indexI, indexJ;
	if(i == j) {
		if(binLenI <= 1) return;
		indexI = bssa_rand_int(binLenI);
		indexJ = bssa_rand_int(binLenI-1);
		if(indexJ >= indexI) indexJ++;
	}
	else {
		int binLenJ = binLens[j];
		if(binLenI <= 0 || binLenJ <= 0) return;
		indexI = bssa_rand_int(binLenI);
		indexJ = bssa_rand_int(binLenJ);
	}
	
	double kval = bssa_eval_kernel(i, indexI, j, indexJ);
	if(bssa_rand_double()*kmaxRow[j] > kval) return;
	
	int newBin = bssa_coag(i, indexI, j, indexJ);
	//note: j <= i
	bssa_on_part_remove(j, numBins, kmax, binLens, sums);
	if(i != newBin) {
		bssa_on_part_remove(i, numBins, kmax, binLens, sums);
		bssa_on_part_add(newBin, numBins, kmax, binLens, sums);
	}
}

//TODO: possibly always keep track of row sums so that it does not
//      need to be re-computed every timestep...
//TODO: be concerned with possibility of re-binning within a timestep...
void bssa_coags_for_timestep(int numBins, int* binLens, const double* kmax,
		double duration)
{
	double sums[numBins+1];
	bssa_refresh_sums(numBins, kmax, binLens, sums);
	while(1) {
		double rate = sums[numBins];
		if(rate <= 0.0) return;
		duration -= bssa_rand_exp(rate);
		if(duration <= 0.0) break;
		bssa_try_coag(numBins, kmax, binLens, sums);
	}
}
