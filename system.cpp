#include "system.h"
#include <complex>
#include <iostream>

System::System()
{

}

System::System(StepFunction wavefunction, double imaginary_time){
    this->wavefunction = wavefunction;

    delta_t = std::complex<double>(re_delta_t, -re_delta_t*imaginary_time);

    potential_function = StepFunction(wavefunction.step_size, std::vector<std::complex<double> >(wavefunction.vals.size(), 0.));
}

System::System(StepFunction wavefunction, StepFunction potential_function, double imaginary_time){
    this->wavefunction = wavefunction;
    this->potential_function = potential_function;
    delta_t = std::complex<double>(re_delta_t, -re_delta_t*imaginary_time);
}

StepFunction System::hamiltonian(){
    return hamiltonian(wavefunction);
}

StepFunction System::hamiltonian(StepFunction wavefunction){
    auto kinetic = (-1./(2*mass))*wavefunction.second_derivative();
    auto potential = wavefunction*potential_function;

    //StepFunction interaction(wavefunction.step_size,0);

    return kinetic+potential;
}

std::array<StepFunction,3> System::energy_plot(){
    auto kinetic = ((-1./(2*mass))*wavefunction.second_derivative());

    auto potential = potential_function*wavefunction;

    //StepFunction interaction(wavefunction.step_size,0);

    return {(kinetic+potential).abs(),kinetic.abs(),potential_function};
}


void System::rk4_step(){

    StepFunction first=wavefunction, second, third, fourth, final;
    StepFunction first_eval, second_eval, third_eval;

    std::complex<double> i(0.,1.);

    first_eval = -i * hamiltonian(first);

    second = first + delta_t/2.*first_eval;
    second_eval = -i * hamiltonian(second);

    third = first + delta_t/2.*second_eval;
    third_eval = -i * hamiltonian(third);

    fourth = first + delta_t*third_eval;

    final = first + (delta_t/6.)*
                    (first_eval
                     + 2*second_eval
                     + 2*third_eval
                     - i * hamiltonian(fourth));

    wavefunction = (1./std::sqrt(final.norm()))*final;
}
