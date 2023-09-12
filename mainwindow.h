#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "qcustomplot.h"
#include "system.h"
#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow();

    void plot();

private slots:
    void on_startButton_clicked();

    void on_stopButton_clicked();
    void sim_move();

    void on_harmonic_potential_clicked();
    void on_step_potential_clicked();
    void on_free_particle_potential_clicked();
    void on_mexican_hat_potential_clicked();
    void on_gauss_clicked();
    void on_gauss_sine_clicked();
    void on_gauss_planewave_clicked();
    void on_gauss_pair_clicked();

    void on_potential_a_valueChanged(double);
    void on_potential_b_valueChanged(double);
    void on_potential_x0_valueChanged(double);
    void on_potential_x1_valueChanged(double);
    void on_wavefunction_s_valueChanged(double);
    void on_wavefunction_k_valueChanged(double);
    void on_wavefunction_x0_valueChanged(double);
    void on_wavefunction_x1_valueChanged(double);

    void on_imaginaryTimeSpinBox_valueChanged(double);

    void on_linear_potential_clicked();

private:
    Ui::MainWindow* ui;

    System system;
    std::unique_ptr<QTimer> timer;
    size_t step = 0;

    StepFunction harmonic_potential, no_potential, mexican_hat, step_potential, linear_potential;
    StepFunction gauss, gauss_sine, gauss_plane_wave, gauss_pair;
    double step_size = 10 / 512.;

    void recalculate_functions();
    void restart_system();
};

#endif // MAINWINDOW_H
