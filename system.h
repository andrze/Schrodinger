#ifndef SYSTEM_H
#define SYSTEM_H
#include "stepfunction.h"

class System
{
public:
    System();
    System(StepFunction wavefunction);
    System(StepFunction wavefunction, StepFunction potential_function);

    StepFunction wavefunction, potential_function;
    double mass=1;
    StepFunction hamiltonian();
    StepFunction hamiltonian(StepFunction wavefunction);
    std::array<StepFunction,3> energy_plot();
    double delta_t = 0.00002;
    void rk4_step();


private:

};


#endif // SYSTEM_H
