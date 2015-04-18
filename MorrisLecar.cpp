/*
MorrisLecar.cpp
Models a Morris Lecar neuron. The parameters here result in Type II  behavior.
*/

#include "MorrisLecar.h"
#ifdef __AVX__
#include <emmintrin.h>
#include <immintrin.h>
#endif

MorrisLecar::MorrisLecar() {
	eq_count = 2;
	
	diff_eq = (double (**)(double *, double)) malloc(eq_count * sizeof(void *));
	diff_eq[0] = &MorrisLecar::dV;
	diff_eq[1] = &MorrisLecar::dw;
	
	state = (double *) malloc(eq_count * sizeof(double));
}

MorrisLecar::~MorrisLecar() {
	free(diff_eq);
	free(state);
}

double MorrisLecar::dV(double *V, double I) {
	const double C = 1.0;
	const double gCa = 1.1;
	const double gK = 2.0;
	const double gL = 0.5;
	const double ECa = 100.0;
	const double EK = -70.0;
	const double EL = -50.0;
#ifdef __AVX__
/*
AVX is an instruction set from Intel which allows simultaneous operation
on 4 doubles. Use it if we have it.
*/
	double Va[] __attribute__ ((aligned (32))) = {V[0], V[0], V[0], 1.0},
		   Ea[] __attribute__ ((aligned (32))) = {EL, ECa, EK, 0.0},
		   Ga[] __attribute__ ((aligned (32))) = {-gL, -gCa * minf(V[0]), -gK * V[1], I};
	
	// load V
	__m256d Vr = _mm256_load_pd(Va);
	// load E
	__m256d Er = _mm256_load_pd(Ea);
	// load G
	__m256d Gr = _mm256_load_pd(Ga);
	// subtract
	Vr = _mm256_sub_pd(Vr, Er);
	// dot product (why does intel not have _mm256_dp_pd ?)
	Vr = _mm256_mul_pd(Vr, Gr);
	__m256d temp = _mm256_hadd_pd(Vr, Vr);
	__m128d lo128 = _mm256_extractf128_pd(temp, 0);
	__m128d hi128 = _mm256_extractf128_pd(temp, 1);
	__m128d dotproduct = _mm_add_pd(lo128, hi128);
	
	double sseVal;
	// store
	_mm_storel_pd(&sseVal, dotproduct);
	sseVal /= C;
	
	return sseVal;
#else	
	return (-gL * (V[0] - EL) - gCa * minf(V[0]) * (V[0] - ECa)
		 - gK * V[1] * (V[0] - EK) + I) / C;
#endif
}

double MorrisLecar::dw(double *V, double I) {
	const double phi = 0.2;

	return  (phi * (winf(V[0]) - V[1])/tauw(V[0]));;
}

double MorrisLecar::minf(double pot) {
	const double V1 = -1.0;
	const double V2 = 15.0;

	return 0.5 * (1.0 + tanh((pot - V1) / V2));
}

double MorrisLecar::winf(double pot) {
	const double V3 = 0.0;
	const double V4 = 30.0;

	return 0.5 * (1.0 + tanh((pot - V3) / V4));
}

double MorrisLecar::tauw(double pot) {
	const double V3 = 0.0;
	const double V4 = 30.0;

	return 1.0 / cosh((pot - V3) / (2.0 * V4));
}
