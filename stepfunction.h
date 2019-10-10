#ifndef STEPFUNCTION_H
#define STEPFUNCTION_H
#include <vector>
#include <map>
#include <complex>

class StepFunction
{
public:
    StepFunction();
    StepFunction(double step_size, std::vector<std::complex<double> > y);

    double step_size;
    double domain;
    std::vector<std::complex<double> > vals;
    std::vector<double> xs();

    std::complex<double> operator()(double x);
    std::complex<double>& operator[](int i);
    StepFunction derivative();
    StepFunction second_derivative();
    StepFunction conjugate();
    StepFunction abs();
    StepFunction norm_function();
    StepFunction phase();
    StepFunction fourier_transform();
    StepFunction fast_fourier_transform();
    std::vector<double> re();
    std::vector<double> im();

    std::complex<double> integral();
    double norm();

};

StepFunction operator+(StepFunction lhs, StepFunction rhs);
StepFunction operator*(StepFunction lhs, StepFunction rhs);
StepFunction operator-(StepFunction lhs, StepFunction rhs);

StepFunction operator*(double lhs, StepFunction rhs);
StepFunction operator*(std::complex<double> lhs, StepFunction rhs);
StepFunction operator/(std::complex<double> lhs, StepFunction rhs);
StepFunction operator/(StepFunction lhs, StepFunction rhs);

#endif // STEPFUNCTION_H
