#ifndef _MORRIS_LECAR_H_
#define _MORRIS_LECAR_H_

#include "NeuronModel.h"
#include <stdlib.h>
#include <math.h>

class MorrisLecar: public NeuronModel {
public:	
	MorrisLecar();
	~MorrisLecar();
		
private:	
	static double dV(double *V, double I);
	static double dw(double *V, double I);
	
	static double minf(double pot);
	static double winf(double pot);
	static double tauw(double pot);
};

#endif
