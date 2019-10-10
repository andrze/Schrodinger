#include <QApplication>
#include "realvector.h"
#include "mainwindow.h"
#include <fenv.h>
#include "simpson.h"
#include <iostream>
#include <cmath>

int main(int argc, char *argv[])
{
    feraiseexcept(FE_INVALID | FE_OVERFLOW);

    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
