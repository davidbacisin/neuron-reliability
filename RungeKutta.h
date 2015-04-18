#ifndef _RUNGEKUTTA_H_
#define _RUNGEKUTTA_H_

#include "NeuronModel.h"
#include <stdlib.h>

void RungeKutta(NeuronModel *model, double I, double stochastic);

#endif