/*
Count the number of spikes in a 100 ms time window
for each value of I0. Convert that count to a frequency
and output the frequency.

If sigma (the noise factor) is 0, then trials should be 1.
Otherwise you're just repeating yourself.
*/

#include "MorrisLecar.h"
#include "RungeKutta.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#define trials 100
#define sigma 0.4

// interval between the tested I0 values
#define I0_step 0.05

#define dt 0.1
#define t_max 200.0
#define t_ignore 100.0

float gasdev(long *idum);

void runTrials(MorrisLecar& neuron,
			   std::vector<double>& psth,
			   std::vector<long>& nspikes,
			   double I0,
			   double I1, 
			   double f,
			   long& seed) {
	double I,
		   Inoise,
		   t;
	
	int trial;
	
	double V_prev,
		   V_threshold = 20.0;
		   	
	// iterate over trials
	for (trial = 0; trial < trials; trial++) {
		// set the starting state
		neuron.state[0] = -21.35; //4.64501;
		neuron.state[1] = 0.2; //0.43051;
		// run each simulation for t_max ms
		for (t=0; t < t_max; t += dt){
			// calculate the input current for the current time
			double noise = (double) gasdev(&seed);
			Inoise = sigma * sqrt(dt) * noise;
			I = I0 + I1 * sin(f * t) + Inoise;
		
			RungeKutta(&neuron, I, 0.0);
			//fprintf(vlog, "%f\t%f\t%f\t%f\t%f\n", I0, t, Inoise, neuron.state[0], neuron.state[1]);
			
			// increment spike count if we're over the threshold
			// but the previous value was under the threshold
			// (we don't want to count the same spike multiple times)
			if (neuron.state[0] > V_threshold &&
				V_prev < V_threshold &&
				t > t_ignore) { // ignore the spikes before t_ignore
				psth[(int) (t/t_max)]++;
				nspikes[trial]++;
			}
			
			V_prev = neuron.state[0];				
		}
	}
}

int main() {	
	double freq,
		I1,
		f,
		avg_spikes;
	long seed = -585;
	MorrisLecar neuron;
	
	long nbins = (int) (t_max / dt);
	std::vector<double> psth(nbins, 0.0);
	std::vector<long> nspikes(trials, 0);
	
	// open our output file
	FILE *log = fopen("ffir_vs_I0.log", "w");
	if (log == NULL) {
		printf("Cannot open the file.\n");
		exit(1);
	}
	
	// print column headings
	fprintf(log, "I0\tfreq\n");
	
	// iterate over f
	for(freq = 0.0; freq <= 120.0; freq += 0.1) {
		// iterate over I1
		for (I1 = 0.0; I1 <= 1.0; I1 += 0.02) {
			// convert freq into something usable for our calculations
			f = 2.0 * M_PI * freq / 1000.0;
			
			// reset spike and psth counters
			psth.assign(nbins, 0.0);
			nspikes.assign(trials, 0);
			
			// run the trials
			runTrials(neuron, psth, nspikes, I0, I1, f, seed);
			
			// find the total variance
			double total_variance = 0.0,
				   avg_psth = 0.0,
				   aux;
			for (int i = 0; i < nbins; i++) {
				aux = ((double) psth[i])/((double) trials);
				total_variance += aux * aux;
				avg_psth += aux;
			}
			avg_psth /= (double) nbins;
			total_variance = (total_variance/(double) nbins) - avg_psth * avg_psth;
			
			// find the average variance?
			double avg_variance = 0.0;
			for (int i = 0; i < trials; i++) {
				aux = ((double) nspikes[i])/((double) nbins);
				aux = aux * (1.0 - aux) / ((double) trials);
				avg_variance += aux;
			}
		}
	}
	
	// close the file
	fclose(log);
	
	return 0;
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

float gasdev(long *idum)
{
	//float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}