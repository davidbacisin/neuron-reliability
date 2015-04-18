/*

*/
#ifdef HH
#include "HodgkinHuxley.h"
#define NEURON_MODEL HodgkinHuxley
#else
#include "MorrisLecar.h"
#define NEURON_MODEL MorrisLecar
#endif

#include "RungeKutta.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#ifdef WIN32
#include <windows.h>
#endif

#define trials 5
#define sigma 0.3

// intervals
#define I0_step 0.01
#define I1_step 0.1
#define f_step 1

#ifndef dt // should be defined in RungeKutta.h
#define dt 0.1
#endif

#ifdef HH
#define t_pre 25.0
#define t_measure 50.0
#define t_post 25.0
#else
#define t_pre 50.0
#define t_measure 100.0
#define t_post 50.0
#endif

#ifdef WIN32
#define THREAD_COUNT 1
#else
#define THREAD_COUNT 1
#endif

FILE *logfile,
	 *vlogfile;
float gasdev(long *idum);

struct TrialParameters {
	NEURON_MODEL neuron;
	long seed;
	float I0,
		  I1,
		  f;
};

double getReliability(TrialParameters *tp) {
	double t,
		   I,
		   Inoise,
		   reliability = 0.0;
	long trial,
		 events,
		 subpass,
		 binSubpass = 1.0/dt,
		 bin,
		 binCount = t_measure;
	
	double V_prev = 0.0,
		   V_threshold = 20.0;

	int unstable_trials = 0;
		   
	std::vector<long> psth(binCount);
	std::vector<long> nspikes(trials);
	
	// iterate over trials
	for (trial = 0; trial < trials; trial++) {
		nspikes[trial] = 0;
		events = 0;
		// run each simulation for t_measure ms
		for (t=0; t < t_pre + t_measure + t_post; t += dt, subpass++){
			// calculate input current for the current time
			Inoise = sigma * sqrt(dt) * gasdev(&(tp->seed));
			// divide by 1000 to convert t from ms to sec
			I = tp->I0 + tp->I1 * sin(2.0 / 1000.0 * M_PI * tp->f * t);
		
			RungeKutta(&(tp->neuron), I, Inoise);

			if (isnan(tp->neuron.state[0])){ 
				unstable_trials++;
				break;
			}
			
			//fprintf(vlogfile, "%lf\t%lf\n", t, tp->neuron.state[0]);
			//printf("%lf\t%lf\n", t, tp->neuron.state[0]);
			
			// increment spike count if we're over the threshold
			// but the previous value was under the threshold
			// (we don't want to count the same spike multiple times)
			if (tp->neuron.state[0] > V_threshold &&
				V_prev < V_threshold &&
				t > t_pre &&
				t < t_pre + t_measure) { // count only the spikes in the middle
				events++;
				// printf("Spike!\n");
			}
			
			if (subpass == binSubpass) {
				if(events > 0) {
					bin = (long) (t - t_pre) / (binSubpass * dt);
					if (bin > 0 && bin < binCount) {
						psth[bin]++;
						nspikes[trial]++;
					}
				}
				subpass = 0;
				events = 0;
			}
			
			V_prev = tp->neuron.state[0];				
		}
	}

	if (unstable_trials) {
		printf("%d trials became unstable!\n", unstable_trials);
		return -1;
	}
	
	// now calculate reliability
	double totalVariance = 0.0,
		   psthAvg = 0.0,
		   binAvg;
	for (long i = 1; i < binCount; i++) {
		binAvg = (double) (psth[i])/(double) trials;
		totalVariance += binAvg	* binAvg;
		psthAvg += binAvg;
	}
	psthAvg /= (double) binCount;
	totalVariance /= (double) binCount;
	totalVariance -= psthAvg * psthAvg;
	
	
	double trialVariance = 0.0;
	for(long i = 1; i < trials; i++) {
		binAvg = (double) (nspikes[i])/(double) trials;
		// printf("avg nspikes[bin]: %f\n", binAvg);
		trialVariance += binAvg * (1.0 - binAvg) / ((double) trials);
	}
	
	if (trialVariance > 0.0) {
		reliability = totalVariance/trialVariance;
	}
	else {
		reliability = -1.0; // error, divide by zero
	}
	
	return reliability;
}

#ifdef WIN32
DWORD WINAPI runThread(void *threadData) {
#else
void runThread(void *threadData) {
#endif
	TrialParameters *tp = (TrialParameters *) threadData;
	
	double reliability = getReliability(tp);
	fprintf(logfile, "%lf\t%lf\t%lf\n", tp->I1, tp->f, reliability);
	printf("%lf\t%lf\t%lf\n", tp->I1, tp->f, reliability);

#ifdef WIN32
	return 0;
#endif
}

int main() {
	float I0,
		I1,
		f;
	// prepare for running multiple trials simultaneously
	/*
	MorrisLecar:
		Subthreshold:
			state[0] = -21.35
			state[1] = 0.2
		Suprathreshold:
			state[0] = 4.65
			state[1] = 0.43
	HodgkinHuxley:
		Subthreshold:
			state[0] = -67.2
			state[1] = 0.31
			state[2] = 0.26
			state[3] = 0.01
		Suprathreshold:
			state[0] = 1.0
			state[1] = 0.052 // 0.317
			state[2] = 0.596 // 0.052
			state[3] = 0.317 // 0.596
	*/
	TrialParameters trialData[THREAD_COUNT];
	for (int th=0; th < THREAD_COUNT; th++ ) {
		#ifdef HH
		trialData[th].neuron.state[0] = -67.2;
		trialData[th].neuron.state[1] = 0.31;
		trialData[th].neuron.state[2] = 0.26;
		trialData[th].neuron.state[3] = 0.01;
		#else
		trialData[th].neuron.state[0] = -21.35;
		trialData[th].neuron.state[1] = 0.2;
		#endif
		trialData[th].seed = -585 + th;
	}
	
	// open our output file
	// reliability_hh_bak_I03.6_supra_trials500_sigma0.3_t50-100-50.log	
	logfile = fopen("data/reliability.log", "w");
	if (logfile == NULL) {
		printf("Cannot open the file.\n");
		exit(1);
	}
	
	// to plot voltage vs time
	vlogfile = fopen("data/voltage.log", "w");
	if (vlogfile == NULL) {
		printf("Cannot open the file.\n");
		exit(1);
	}
	
	// print column headings
	fprintf(logfile, "I1\tfreq\treliability\n");
	fprintf(vlogfile, "t\tV\n");

#ifdef WIN32	
	HANDLE threads[THREAD_COUNT];
#endif
	
	#ifdef HH
	I0 = 3.6;
	#else
	I0 = 25.4; // 25.4 for ML, 3.6 for HH
	#endif
	// iterate over I1 from zero to one
	for(I1 = 0.0; I1 <= 1.0; I1 += I1_step) {
		// iterate over frequency
		for(f = 0.0; f <= 120.0; ) {
			for (int th=0; th < THREAD_COUNT; th++, f += f_step) {
				trialData[th].I0 = I0;
				trialData[th].I1 = I1;
				trialData[th].f = f;
				#ifdef WIN32
				threads[th] = CreateThread(NULL, 0, runThread, &(trialData[th]), 0, NULL);
				#else
				runThread(&(trialData[th]));
				#endif
			}
			#ifdef WIN32
			WaitForMultipleObjects(THREAD_COUNT, threads, TRUE, INFINITE);
			for (int th=0; th < THREAD_COUNT; th++) {
				CloseHandle(threads[th]);
			}
			#endif
		}
	}
	
	// close the files
	fclose(logfile);
	fclose(vlogfile);
	
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
