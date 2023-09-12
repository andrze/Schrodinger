#ifndef SYSTEM_H
#define SYSTEM_H
#include "stepfunction.h"

class System
{
public:
    System();
    System(StepFunction wavefunction, double imaginary_time);
    System(StepFunction wavefunction, StepFunction potential_function, double imaginary_time);

    StepFunction wavefunction, potential_function;
    double mass=1;
    StepFunction hamiltonian();
    StepFunction hamiltonian(StepFunction wavefunction);
    std::array<StepFunction,3> energy_plot();
    double re_delta_t = 1e-4;
    std::complex<double> delta_t = re_delta_t;
    void rk4_step();


private:

};


#endif // SYSTEM_H
