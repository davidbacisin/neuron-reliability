#include "RungeKutta.h"

#define dt 0.1
#define MAX_EQS 4
/*
Runge-Kutta second-order aka the midpoint method
*/

void RungeKutta(NeuronModel *model, double I, double stochastic) {
    int i;

	double K1[MAX_EQS],
	   K2[MAX_EQS],
	   Vaux[MAX_EQS];

	// copy initial values of V into Vaux
    for (i=0; i < model->eq_count; i++) {
        Vaux[i] = model->state[i];
    }
	
	// calculate the initial derivatives
	for (i=0; i < model->eq_count; i++) {
		K1[i] = (model->diff_eq[i])(Vaux, I);
	}

	// update Vaux with the next value using forward Euler
    for (i=0; i < model->eq_count; i++) {
        Vaux[i] = model->state[i] + dt * K1[i];
    }
	// add noise to the potential before calculating the next derivative
    Vaux[0] += stochastic;

	// calculate the next derivatives
	for (i=0; i < model->eq_count; i++) {
		K2[i] = (model->diff_eq[i])(Vaux, I);
	}

	// combine the parts to output the value
    for (i=0; i < model->eq_count; i++) {
        model->state[i] += dt * (K1[i] + K2[i]) / 2.0;
    }
	// add noise to the potential
    model->state[0] += stochastic;
}