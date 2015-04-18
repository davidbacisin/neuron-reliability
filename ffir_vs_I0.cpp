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
#include <windows.h>

#define trials 200
#define sigma 0.4
// #define I0 25.4

// interval between the tested I0 values
#define I0_step 0.01

#define dt 0.1
#define t_pre 100.0
#define t_measure 1000.0
#define t_post 50.0

#define THREAD_COUNT 4

FILE *logfile;
long seed = -585;
float gasdev(long *idum);
double runTrials(MorrisLecar *, double);

struct TrialParameters {
	MorrisLecar neuron;
	double I0;
};

DWORD WINAPI runThread(void *threadData) {
	TrialParameters *tp = (TrialParameters *) threadData;
	
	double avg_spikes = runTrials(&(tp->neuron), tp->I0);
	double freq = (double) 1000.0 * avg_spikes / t_measure;
	fprintf(logfile, "%lf\t%lf\n", tp->I0, freq);
	printf("%lf\t%lf\n", tp->I0, freq);
	
	return 0;
}

double runTrials(MorrisLecar *neuron, double I0) {
	double Inoise;
	double I;
	double t;
	long nspikes;
	double avg_spikes;
	
	int trial;
	
	double V_prev = 0.0,
		   V_threshold = 20.0;
		   	
	avg_spikes = 0;
	// iterate over trials
	for (trial = 0; trial < trials; trial++) {
		nspikes = 0;
		// run each simulation for t_measure ms
		for (t=0; t < t_pre + t_measure + t_post; t += dt){
			// calculate input current for the current time
			float noise = gasdev(&seed);
			Inoise = sigma * sqrt(dt) * noise;
			I = I0;
		
			RungeKutta(neuron, I, Inoise);
			//fprintf(vlogfile, "%f\t%f\t%f\t%f\t%f\n", I0, t, Inoise, neuron.state[0], neuron.state[1]);
			
			// increment spike count if we're over the threshold
			// but the previous value was under the threshold
			// (we don't want to count the same spike multiple times)
			if (neuron->state[0] > V_threshold &&
				V_prev < V_threshold &&
				t > t_pre && t < t_pre + t_measure) { // ignore the spikes before t_ignore
				nspikes++;
			}
			
			V_prev = neuron->state[0];				
		}
		// average the result onto the running average
		avg_spikes = (avg_spikes * trial + nspikes)/(trial+1);
	}
	return avg_spikes;
}

int main() {
	double I0,
		avg_spikes;
	MorrisLecar neuron;
	
	// open our output file
	logfile = fopen("ffir_vs_I0_sub_200t1s.log", "w");
	if (logfile == NULL) {
		printf("Cannot open the file.\n");
		exit(1);
	}
	
	// print column headings
	fprintf(logfile, "I0\tfreq\n");
	// fprintf(logfile, "sigma\tfreq\n");
	
	// set the starting state
	//neuron.state[0] = -21.35; // 4.64501
	//neuron.state[1] = 0.2; // 0.43051
	
	// prepare for running multiple trials simultaneously
	TrialParameters trialData[THREAD_COUNT];
	for (int th=0; th < THREAD_COUNT; th++ ) {
		trialData[th].neuron.state[0] = -21.35;
		trialData[th].neuron.state[1] = 0.2;
	}
		
	HANDLE threads[THREAD_COUNT];
	
	// iterate over I0 from low to high
	for(I0 = 25.0; I0 <= 26.0; ) {
		for (int th=0; th < THREAD_COUNT; th++, I0 += I0_step) {
			trialData[th].I0 = I0;
			threads[th] = CreateThread(NULL, 0, runThread, &(trialData[th]), 0, NULL);
		}
		WaitForMultipleObjects(THREAD_COUNT, threads, TRUE, INFINITE);
		for (int th=0; th < THREAD_COUNT; th++) {
			CloseHandle(threads[th]);
		}
		/*
		avg_spikes = runTrials(&neuron, I0);
		
		// we now have a spike count. Output the frequency in Hz
		double freq = (double) 1000.0 * avg_spikes / t_measure;
		fprintf(logfile, "%lf\t%lf\n", I0, freq);
		printf("%lf\t%lf\n", I0, freq);
		*/
	}
		
	// and back down
	/*for(I0 = 26.0; I0 >= 25.0; I0 -= I0_step) {
		avg_spikes = runTrials(neuron, I0);
		
		// we now have a spike count. Output the frequency in Hz
		double freq = (double) 1000.0 * avg_spikes / t_measure;
		fprintf(logfile, "%lf\t%lf\n", I0, freq);
		printf("%lf\t%lf\n", I0, freq);
	}*/
	
	/*
	for(sigma = 0.0; sigma <= 1.2; sigma += 0.01) {		
		avg_spikes = runTrials(neuron, sigma);
		
		// we now have a spike count. Output the frequency in Hz
		double freq = (double) 1000.0 * avg_spikes / t_measure;
		fprintf(logfile, "%lf\t%lf\n", sigma, freq);
		printf("%lf\t%lf\n", sigma, freq);
	}
	*/
	
	// close the file
	fclose(logfile);
	
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