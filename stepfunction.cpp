#include "stepfunction.h"
#include <limits>
#include <iostream>

StepFunction::StepFunction()
{

}

StepFunction::StepFunction(double step_size, std::vector<std::complex<double> > y)
{
    this->step_size = step_size;
    vals = y;
    domain = step_size*vals.size();
}

std::vector<double> StepFunction::xs(){
    std::vector<double> xs;
    double x = -step_size*vals.size()/2.;
    for (size_t i=0; i<vals.size(); i++) {
        xs.push_back(x);
        x += step_size;
    }
    return xs;
}

std::complex<double> StepFunction::operator()(double x){
    if(x < 0){
        x += std::ceil(std::abs(x)/domain)*domain;
        if(x<0){
            throw std::runtime_error("");
        }
    }

    size_t pos = size_t(x/step_size) % vals.size();

    double offset = x/step_size-pos;
    if(offset < 0.01){
        return vals[pos];
    }

    if(pos+1>vals.size()){
        throw std::invalid_argument("Function argument out of range");
    }

    std::complex<double> val = vals[pos]*(1-offset) + vals[pos+1]*offset;

    return val;
}

std::complex<double>& StepFunction::operator[](int i){

    size_t pos = size_t(i) & (vals.size()-1);

    return vals[pos];
}


StepFunction StepFunction::derivative(){
    std::vector<std::complex<double> > new_vals;

    for(int i=0; i<int(vals.size()); i++){
        new_vals.push_back(((*this)[i-2] - 8.*(*this)[i-1] + 8.*(*this)[i+1] - (*this)[i+2])/(12*step_size));
    }

    return StepFunction(step_size, new_vals);
}


StepFunction StepFunction::second_derivative(){

    auto second = derivative().derivative();
    std::vector<std::complex<double> > new_vals;

    new_vals.reserve(vals.size());

    double h2 = step_size*step_size;
    std::array<double,7> coeficients{1./90, -3./20, 3./2, -49./18, 3./2, -3./20, 1./90};
    std::array<double,9> coeficients2{-1./560, 8./315, -1./5, 8./5, -205./72, 8./5, -1./5, 8./315, -1./560};

    for(auto& c: coeficients){
        c /= h2;
    }

    for(int i=0; i<int(vals.size()); i++){
        std::complex<double> val = 0.;
        for(int j=0; j<7; j++){
            val += coeficients[j]*(*this)[i-3+j];
        }
        new_vals.push_back(val);
        //new_vals.push_back((-(*this)[i-2] + 16.*(*this)[i-1] -30.*(*this)[i] + 16.*(*this)[i+1] - (*this)[i+2])/(12.*step_size*step_size));
    }
    /*

    auto der = StepFunction(step_size, new_vals);

    second = second - der;

    for (int i=0; i<int(vals.size()); i++) {
        std::cout<<der[i]<<' '<<second[i]<<' '<<i<<std::endl;
    }*/

    return StepFunction(step_size, new_vals);
}

StepFunction StepFunction::conjugate(){

    std::vector<std::complex<double> > new_vals;
    for(auto&& v: vals){
        new_vals.push_back(std::conj(v));
    }

    return StepFunction(step_size, new_vals);
}

StepFunction StepFunction::abs(){
    std::vector<std::complex<double> > new_vals;
    for(auto&& v: vals){
        new_vals.push_back(std::abs(v));
    }

    return StepFunction(step_size, new_vals);
}

StepFunction StepFunction::norm_function(){
    return (*this) * this->conjugate();
}

std::complex<double> StepFunction::integral(){
    std::complex<double> sum=0;

    for(auto&& v: vals){
        sum += v;
    }

    return sum*step_size;
}

double StepFunction::norm(){
    StepFunction norm_function = this->norm_function();
    std::complex<double> norm = norm_function.integral();
    if(std::imag(norm) > std::numeric_limits<double>::epsilon()){
        throw std::runtime_error("Function norm has complex value");
    }

    return std::real(norm);
}

StepFunction StepFunction::phase(){
    std::complex<double> phase(1,0);
    std::vector<std::complex<double> > new_vals;

    for (size_t i=0; i<vals.size(); i++) {
        auto v = vals[i];
        if(std::abs(v) > 1e-16){
            phase = v/std::abs(v);
        }
        new_vals.push_back(v);
    }
    return StepFunction(step_size, new_vals);
}

StepFunction StepFunction::fourier_transform(){
    std::vector<std::complex<double> > new_vals;

    std::complex<double> i(0.,1.);
    double two_pi = 2*M_PI;
    for(size_t k=0; k<vals.size(); k++){
        double freq = two_pi*k/vals.size();
        std::complex<double> xk;
        for (size_t n=0; n<vals.size();n++){
            xk += vals[n]*(std::cos(freq*n) - i *std::sin(freq*n));
        }
        new_vals.push_back(xk);
    }
    return StepFunction(step_size, new_vals);
}

void ditfft2(std::complex<double>* function, std::complex<double>* target,
 size_t N, size_t s){
    if(N == 1){
        *target = *function;
        return;
    }

    ditfft2(function, target, N/2, 2*s);
    ditfft2(function+s, target+N/2, N/2, 2*s);


    double base_freq = 2.*M_PI/double(N);

    for(size_t k=0; k<N/2; k++){
        std::complex<double> t = *(target+k);
        std::complex<double> component = *(target+k+N/2) * std::polar(1.,-base_freq*double(k));

        *(target+k) = t + component;
        *(target+k+N/2) = t - component;
    }
}

StepFunction StepFunction::fast_fourier_transform(){
    std::vector<std::complex<double> > target(vals.size());
    ditfft2(&vals[0], &target[0], vals.size(), 1);

    return step_size/sqrt(2*M_PI)*StepFunction(M_PI/5, target);
}

std::vector<double> StepFunction::re(){
    std::vector<double > new_vals;
    for(auto&& v: vals){
        new_vals.push_back(std::real(v));
    }

    return new_vals;
}

std::vector<double> StepFunction::im(){
    std::vector<double > new_vals;
    for(auto&& v: vals){
        new_vals.push_back(std::imag(v));
    }

    return new_vals;
}

StepFunction operator+(StepFunction lhs, StepFunction rhs){
    if(std::abs(lhs.step_size - rhs.step_size) > std::numeric_limits<double>::epsilon()){
        throw std::invalid_argument("Step sizes of added functions should be equal");
    }
    if(lhs.vals.size() != rhs.vals.size()){
        throw std::invalid_argument("Values sizes of added functions should be equal");
    }
    std::vector<std::complex<double> > new_vals;

    for(size_t i=0; i<lhs.vals.size(); i++){
        new_vals.push_back(lhs.vals[i]+rhs.vals[i]);
    }

    return StepFunction(lhs.step_size, new_vals);
}

StepFunction operator*(StepFunction lhs, StepFunction rhs){
    if(std::abs(lhs.step_size - rhs.step_size) > std::numeric_limits<double>::epsilon()){
        throw std::invalid_argument("Step sizes of multiplied functions should be equal");
    }
    if(lhs.vals.size() != rhs.vals.size()){
        throw std::invalid_argument("Values sizes of multiplied functions should be equal");
    }
    std::vector<std::complex<double> > new_vals;

    for(size_t i=0; i<lhs.vals.size(); i++){
        new_vals.push_back(lhs.vals[i]*rhs.vals[i]);
    }

    return StepFunction(lhs.step_size, new_vals);
}

StepFunction operator-(StepFunction lhs, StepFunction rhs){

    return lhs + (-1.)*rhs;
}

StepFunction operator*(std::complex<double> lhs, StepFunction rhs){
    std::vector<std::complex<double> > new_vals;

    for(size_t i=0; i<rhs.vals.size(); i++){
        new_vals.push_back(lhs*rhs.vals[i]);
    }

    return StepFunction(rhs.step_size, new_vals);
}

StepFunction operator/(std::complex<double> lhs, StepFunction rhs){
    std::vector<std::complex<double> > new_vals;

    for(size_t i=0; i<rhs.vals.size(); i++){
        new_vals.push_back(lhs/rhs.vals[i]);
    }

    return StepFunction(rhs.step_size, new_vals);
}

StepFunction operator/(StepFunction lhs, StepFunction rhs){
    return lhs * (1./rhs);
}

StepFunction operator*(double lhs, StepFunction rhs){
    std::complex<double> lhs2 = lhs;
    return lhs2*rhs;
}
