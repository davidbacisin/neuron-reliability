#ifndef _NEURON_MODEL_H_
#define _NEURON_MODEL_H_

class NeuronModel {
public:
	double (**diff_eq)(double *V, double I);
	double *state;
	int eq_count;
};

#endif