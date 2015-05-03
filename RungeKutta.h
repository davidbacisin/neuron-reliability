#ifndef _RUNGEKUTTA_H_
#define _RUNGEKUTTA_H_

#include "NeuronModel.h"
#include <stdlib.h>

#define dt 0.0006
// #define dt 0.1

void RungeKutta(NeuronModel *model, double I, double stochastic);

#endif
