#include <math.h>

//return kernel value between particles at binI-indexI and binJ-indexJ
double bssa_eval_kernel(int binI, int indexI, int binJ, int indexJ);

//coagulate two particles, possibly using weighted scheme
//-binI: bin of first particle
//-indexI: index of first particle within binI
//-binJ: bin of second particle
//-indexJ: index of second particle within binJ
//-usedI: true iff first particle was removed
//-usedJ: true iff second particle was removed
//-newBin: bin of coagulated particle, or -1 if new particle was not made,
// or -2 if bin structure changed
// (if non-negative, newBin >= max(binI, binJ) is expected)
void bssa_coag(int binI, int indexI, int binJ, int indexJ
		int& usedI, int& usedJ, int& newBin);

//uniform between 0 and 1
double bssa_rand_double();

//uniform from {0, 1, ..., range-1}
int bssa_rand_int(int range);

double bssa_rand_exp(double mean)
{
	return -log(bssa_rand_double())/mean;
}

double bssa_num_pairs(int sizeI, int sizeJ)
{
	return sizeI*(double)sizeJ;
}

double bssa_num_self_pairs(int size)
{
	return size*.5*(size-1);
}

void bssa_refresh_sums(int numBins, const double* kmax,
		const int* binPops, double* sums)
{
	double netSum = 0.0;
	for(int i = 0; i < numBins; i++) {
		double* kmaxRow = kmax + i*numBins;
		double rowSum = 0.0;
		int lenI = binPops[i];
		for(int j = 0; j < i; j++) {
			rowSum += kmaxRow[j]*bssa_num_pairs(lenI, binPops[j]);
		}
		rowSum += kmaxRow[i]*bssa_num_self_pairs(lenI);
		sums[i] = rowSum;
		netSum += rowSum;
	}
	sums[numBins] = netSum;
}

void bssa_on_part_add(int bin, int numBins, const double* kmax,
		int* binPops, double* sums)
{
	double* kmaxRow = kmax + bin*numBins;
	double netDiff = 0.0;
	for(int j = 0; j < bin; j++) {
		netDiff += binPops[j]*kmaxRow[j];
	}
	binPops[bin]++;
	int selfPairsDiff = binPops[bin] - 1;
	netDiff += selfPairsDiff*kmaxRow[bin];
	sums[bin] += netDiff;
	for(int i = bin+1; i < numBins; i++) {
		double diff = binPops[i]*kmaxRow[i];
		sums[i] += diff;
		netDiff += diff;
	}
	sums[numBins] += netDiff;
}

void bssa_on_part_remove(int bin, int numBins, const double* kmax,
		int* binPops, double* sums)
{
	double* kmaxRow = kmax + bin*numBins;
	double netDiff = 0.0;
	for(int j = 0; j < bin; j++) {
		netDiff -= binPops[j]*kmaxRow[j];
	}
	binPops[bin]--;
	int selfPairsDiff = -binPops[bin];
	netDiff += selfPairsDiff*kmaxRow[bin];
	sums[bin] += netDiff;
	for(int i = bin+1; i < numBins; i++) {
		double diff = -binPops[i]*kmaxRow[i];
		sums[i] += diff;
		netDiff += diff;
	}
	sums[numBins] += netDiff;
}

void bssa_select_bin_pair(int numBins, const double* kmax,
		const int* binPops, const double* sums, int& binI, int& binJ)
{
	double randVal = bssa_rand_double()*sums[numBins];
	for(binI = numBins - 1; binI > 0; binI--) {
		if(randVal < sums[binI]) break;
		randVal -= sums[binI];
	}
	double* kmaxRow = kmax + binI*numBins;
	int binPopI = binPops[binI];
	for(binJ = 0; binJ < binI; binJ++) {
		randVal -= kmaxRow[binJ]*bssa_num_pairs(binPopI, binPops[binJ]);
		if(randVal < 0.0) break;
	}
}

//return 1 if can coagulate, 0 otherwise
int bssa_select_index_pair(int popI, int popJ, int isSameBin,
		int& indexI, int& indexJ)
{
	if(isSameBin) {
		if(popI <= 1) return 0;
		indexI = bssa_rand_int(popI);
		indexJ = bssa_rand_int(popI-1);
		if(indexJ == indexI) indexJ = popI-1;
	}
	else {
		if(popI <= 0 || popJ <= 0) return 0;
		indexI = bssa_rand_int(popI);
		indexJ = bssa_rand_int(popJ);
	}
	return 1;
}

//return 0 iff bin structure remains the same
int bssa_try_coag(int numBins, const double* kmax, int* binPops, double* sums)
{
	int binI, binJ, indexI, indexJ;
	bssa_select_bin_pair(numBins, kmax, binPops, sums, binI, binJ);
	canCoag = bssa_select_index_pair(binPops[binI], binPops[binJ],
			binI == binJ, indexI, indexJ);
	if(!canCoag) return 0;
	
	double* kmaxRow = kmax + binI*numBins;
	double kval = bssa_eval_kernel(i, indexI, j, indexJ);
	if(bssa_rand_double()*kmaxRow[j] > kval) return 0;
	
	int usedI, usedJ, newBin;
	bssa_coag(i, indexI, j, indexJ, usedI, usedJ, newBin);
	if(newBin <= -2) return 1;
	
	//note: j <= i
	if(usedJ) {
		if(!usedI && j == newBin) return 0;
		bssa_on_part_remove(j, numBins, kmax, binPops, sums);
	}
	if(usedI) {
		if(i == newBin) return 0;
		bssa_on_part_remove(i, numBins, kmax, binPops, sums);
	}
	if(newBin >= 0) bssa_on_part_add(newBin, numBins, kmax, binPops, sums);
	return 0;
}

//TODO: possibly always keep track of row sums so that it does not
//      need to be re-computed every timestep...

//Performs coagulations for a given amount of time.
//May stop partway through if bssa_coag returns -1, meaning
//the bin structure has changed.
//-numBins: number of bins
//-kmax: flattened 2-d symmetric matrix of kernel bounds for pairs of bins,
// contains numBins*numBins elements
//-binPops: array of length numBins denoting population of each bin,
// should NOT be updated during calls to bssa_coag
//-duration: amount of time to be simulated, which will be zero at the
// end if bin structure did not change, and otherwise is the time
// remaining to be simulated
void bssa_coags_for_timestep(int numBins, const double* kmax, int* binPops,
		double& duration)
{
	if(numBins <= 0) {
		duration = 0.0;
		return;
	}
	double sums[numBins+1];
	bssa_refresh_sums(numBins, kmax, binPops, sums);
	while(1) {
		double rate = sums[numBins];
		if(rate <= 0.0) break;
		duration -= bssa_rand_exp(rate);
		if(duration <= 0.0) break;
		int binsChanged = bssa_try_coag(numBins, kmax, binPops, sums);
		if(binsChanged) return;
	}
	duration = 0.0;
}
