#include "qcp_extension.h"
#include "qcustomplot.h"
#include <cmath>
#include <limits>


static std::vector<QPen> color{QPen(Qt::blue), QPen(Qt::red), QPen(Qt::darkGreen),
                               QPen(Qt::darkYellow), QPen(Qt::darkMagenta), QPen(Qt::cyan)};

QCPRange get_x_range(QCustomPlot* plot){
    double mini = +std::numeric_limits<double>::infinity();
    double maxi = -std::numeric_limits<double>::infinity();

    for(int i=0; i<plot->plottableCount(); i++){
        bool found_range;
        auto p = plot->plottable(i);
        QCPRange range = p->getKeyRange(found_range);

        if(found_range){
            mini = std::min(mini, range.lower);
            maxi = std::max(maxi, range.upper);
        }
    }

    return QCPRange(mini, maxi);
}


void rescale_axes(QCustomPlot* plot, double x_lower, double x_upper, bool min_at_zero){
    double mini = +std::numeric_limits<double>::infinity();
    double maxi = -std::numeric_limits<double>::infinity();

    auto current_range = plot->yAxis->range();
    bool change_scales=false;

    bool found_any = false;

    auto key_range = QCPRange(x_lower, x_upper);
    for(int i=0; i<plot->plottableCount(); i++){
        bool found_range;
        auto p = plot->plottable(i);
        QCPRange range;
        range = p->getValueRange(found_range, QCP::sdBoth, key_range);
        if(found_range){
            mini = std::min(mini, range.lower);

            double limit = std::min(400., range.upper);
            maxi = std::max(maxi, limit);
        }
        found_any = found_range || found_any;
    }


    if(min_at_zero){
        mini=0.;
    } else {
        maxi = std::max(maxi, abs(mini));
        mini = -maxi;
    }

    double upper_margin=1.5, lower_margin=1.5;
    change_scales = (maxi>current_range.upper) || (3*maxi<current_range.upper);
    if(maxi>current_range.upper){
        upper_margin = 1.5;
        lower_margin=1.5;
        if(min_at_zero){
            lower_margin = 0.1;
        }
    }

    if(3*maxi<current_range.upper){
        upper_margin = 0.5;
        lower_margin = 0.5;
        if(min_at_zero){
            lower_margin = 0.1;
        }
    }


    plot->xAxis->setRange(key_range);
    if(found_any && change_scales){
        double diff = maxi-mini;
        plot->yAxis->setRange(mini-diff*lower_margin, maxi+diff*upper_margin);
    }

    plot->replot();
}


void add_graph(QCustomPlot* plot, std::vector<double> xvals, std::vector<double> vals,
              bool logplot, int parts, bool replot){
    QVector<double> x = QVector<double>::fromStdVector(xvals);
    QVector<double> y = QVector<double>::fromStdVector(vals);

    if(logplot){
        plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    } else {
        plot->yAxis->setScaleType(QCPAxis::stLinear);
    }

    int graph_num = plot->plottableCount();
    int limit = 6*parts;
    if(graph_num >= 2*limit){
        for(int i=1; i<=limit; i++){
            plot->removePlottable(graph_num-i);
        }
        graph_num -= limit;
    }

    if (graph_num >= limit){
        graph_num -= limit;
    }

    QCPCurve* curve = new QCPCurve(plot->xAxis, plot->yAxis);

    auto pen = color[size_t(graph_num/parts)];
    pen.setStyle(Qt::SolidLine);
    curve->setPen(pen);
    if(parts>1 && graph_num%2 == 0){
        auto pen = curve->pen();
        pen.setStyle(Qt::DotLine);
        curve->setPen(pen);
    }
    curve->setData(x,y);

    if(logplot){
        QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
        plot->yAxis->setTicker(logTicker);
    }

    if(replot){
        plot->replot();
    }
}
