#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include <limits>
#include "simpson.h"

double integral(std::function<double(double)> func){
    const double h0=0.01;
    double q=0.;
    double result=0.;
    int num_steps=51;
    double cutoff=300;
    double f0 = std::abs(func(q));

    for(int i=1; i<=30; i++){
        double f1 = std::abs(func(q+i*0.1));
        if(f1 > f0){
            f0 = f1;
        }
    }
    //std::cerr << f0 <<'\n';

    for(int i=0; i<400; i++){
        double coef = std::pow(std::abs(func(q)/f0),0.2);
        double h;
        if(coef >= h0*std::numeric_limits<double>::epsilon()){
            h = h0/coef;
        } else {
            h = h0;
        }
        //std::cout<< coef<< ' '<<h<<std::endl;
        double partial_result=func(q);

        if(q+h*num_steps>cutoff){
            num_steps = 1+2*int((cutoff-q)/(2*h));
        }

        for(int j=0; j<num_steps; j++){
            q += h;
            double f = func(q);
            if(j % 2 == 1){
                partial_result += 2*f;
            } else {
                partial_result += 4*f;
            }
        }
        q+=h;
        partial_result += func(q);
        result += partial_result*h/3;

        if(std::abs(func(q))*h/3 <= std::abs(result)*std::numeric_limits<double>::epsilon()){
            break;
        }
        if(q>=cutoff){
            break;
        }
    }

    return result;
}

