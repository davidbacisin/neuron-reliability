#ifndef _HODGKIN_HUXLEY_H_
#define _HODGKIN_HUXLEY_H_

#include "NeuronModel.h"
#include <stdlib.h>
#include <math.h>

/*
V[0] = voltage
V[1] = n gate
V[2] = m gate
V[3] = h gate
*/

class HodgkinHuxley: public NeuronModel {
public:	
	HodgkinHuxley();
	~HodgkinHuxley();
		
private:	
	static double dV(double *V, double I);
	static double dn(double *V, double I);
	static double dm(double *V, double I);
	static double dh(double *V, double I);
};

#endif