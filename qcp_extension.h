#ifndef QCP_EXTENSION_H
#define QCP_EXTENSION_H
#include "qcustomplot.h"

QCPRange get_x_range(QCustomPlot* plot);

void rescale_axes(QCustomPlot* plot, double x_lower, double x_upper, bool min_at_zero=false);

void add_graph(QCustomPlot* plot, std::vector<double> xvals, std::vector<double> vals,
              bool logplot=false, int parts=1, bool replot=true);


#endif // QCP_EXTENSION_H
