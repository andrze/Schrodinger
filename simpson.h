#ifndef SIMPSON_H
#define SIMPSON_H
#include <functional>

double integral(std::function<double(double)> func);

#endif // SIMPSON_H
