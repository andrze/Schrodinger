#include "mainwindow.h"
#include "qcp_extension.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QString>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <thread>
#include <vector>

std::vector<std::complex<double>> normalize(std::vector<std::complex<double>> vals)
{
    double mini = std::numeric_limits<double>::max();
    for (auto&& v : vals) {
        mini = std::min(std::real(v), mini);
    }
    for (size_t i = 0; i < vals.size(); i++) {
        vals[i] -= mini;
    }
    return vals;
}

StepFunction normalize(StepFunction wf)
{
    return 1 / std::sqrt(wf.norm()) * wf;
}

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->harmonic_potential->clicked();
    ui->gauss->clicked();

    timer = std::make_unique<QTimer>(this);
    connect(timer.get(), SIGNAL(timeout()), this, SLOT(sim_move()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::sim_move()
{
    for (size_t i = 0; i < 5; i++) {
        system.rk4_step();
        step++;
    }
    if (step % 40 == 0) {
        plot();
    }
}

void MainWindow::plot()
{
    auto wave = system.wavefunction;
    auto fourier = wave.fast_fourier_transform();
    auto f_norm = fourier.norm_function();
    auto energy = system.energy_plot();
    std::vector<double> xs = wave.xs(),
                        real_y = wave.re(),
                        im_y = wave.im(),
                        dens = wave.norm_function().re();

    ui->re_plot->clearPlottables();
    add_graph(ui->re_plot, xs, real_y, false, 1, false);
    add_graph(ui->re_plot, xs, im_y, false, 1, false);
    rescale_axes(ui->re_plot, -5, 5);

    ui->dens_plot->clearPlottables();
    add_graph(ui->dens_plot, xs, dens, false, 1, false);
    rescale_axes(ui->dens_plot, -5, 5, true);

    ui->energy_plot->clearPlottables();
    add_graph(ui->energy_plot, xs, energy[0].re(), false, 1, false);
    add_graph(ui->energy_plot, xs, energy[1].re(), false, 1, false);
    add_graph(ui->energy_plot, xs, energy[2].re(), false, 1, false);
    rescale_axes(ui->energy_plot, -5, 5, true);

    if (step % 10 == 0) {
        ui->fourier_plot->clearPlottables();
        std::vector<double> f_vals;
        int f_size = int(fourier.vals.size());
        for (int i = 0; i < f_size; i++) {
            if (i < f_size / 2) {
                f_vals.push_back(std::abs(f_norm[i + f_size / 2]));
            } else {
                f_vals.push_back(std::abs(f_norm[i - f_size / 2]));
            }
        }
        add_graph(ui->fourier_plot, fourier.xs(), f_vals, false, 1, false);
        rescale_axes(ui->fourier_plot, -fourier.domain / 16, fourier.domain / 16, true);
    }

    ui->spatial_density->setValue(system.wavefunction.norm());
    ui->momentum_density->setValue(fourier.norm());
}

void MainWindow::on_startButton_clicked()
{
    restart_system();

    step = 0;
    timer->start(1);
}

void MainWindow::on_stopButton_clicked()
{
    if (timer->isActive()) {
        timer->stop();
        ui->stopButton->setText(QString("Wznów"));
    } else {
        timer->start();
        ui->stopButton->setText(QString("Zatrzymaj"));
    }
}

void MainWindow::recalculate_functions()
{
    double a = ui->potential_a->value();
    double b = ui->potential_b->value();
    double x0 = ui->potential_x0->value();
    double x1 = ui->potential_x1->value();

    double s = ui->wavefunction_s->value();
    double k = ui->wavefunction_k->value();
    double wf_x0 = ui->wavefunction_x0->value();

    std::vector<std::complex<double>> gauss, gauss_sine, gauss_plane_wave, two_gausses;
    std::vector<std::complex<double>> parabolic_potential, tunneling_potential, no_potential, mexican_potential, linear_potential;

    for (double x = -5; x < 5; x += step_size) {
        std::complex<double> i(0., 1.);
        std::complex<double> x2 = x * x,
                            delta_x2 = (x-wf_x0) * (x-wf_x0),
                             ex2 = std::exp(-s * delta_x2);

        gauss.push_back(ex2);
        gauss_sine.push_back(ex2 * std::sin(M_PI * x * k));
        gauss_plane_wave.push_back(ex2 * std::exp(i * M_PI * x * k));
        two_gausses.push_back(std::exp(-s * std::pow(x - wf_x0, 2)) + std::exp(-s * std::pow(x + wf_x0, 2)));

        parabolic_potential.push_back(a * x2);
        mexican_potential.push_back(a * std::pow(x, 4) - b * std::pow(x, 2));
        no_potential.push_back(0.);

        if (x < -4.6 || x > 4.6) {
            linear_potential.push_back(10000.);
            tunneling_potential.push_back(10000.);
        } else {
            linear_potential.push_back((x+4.6)*a);

            if (x > x0 && x < x1) {
                tunneling_potential.push_back(a);
            } else {
                tunneling_potential.push_back(0.);
            }
        }
    }

    mexican_hat = StepFunction(step_size, normalize(mexican_potential));
    harmonic_potential = StepFunction(step_size, normalize(parabolic_potential));
    this->no_potential = StepFunction(step_size, normalize(no_potential));
    step_potential = StepFunction(step_size, normalize(tunneling_potential));
    this->linear_potential = StepFunction(step_size, normalize(linear_potential));

    this->gauss = normalize(StepFunction(step_size, gauss));
    this->gauss_sine = normalize(StepFunction(step_size, gauss_sine));
    this->gauss_plane_wave = normalize(StepFunction(step_size, gauss_plane_wave));
    gauss_pair = normalize(StepFunction(step_size, two_gausses));
}

void MainWindow::restart_system()
{
    StepFunction potential, wavefunction;

    recalculate_functions();

    if (ui->free_particle_potential->isChecked()) {
        potential = no_potential;
    } else if (ui->harmonic_potential->isChecked()) {
        potential = harmonic_potential;
    } else if (ui->mexican_hat_potential->isChecked()) {
        potential = mexican_hat;
    } else if (ui->step_potential->isChecked()){
        potential = step_potential;
    } else {
        potential = linear_potential;
    }

    if (ui->gauss->isChecked()) {
        wavefunction = gauss;
    } else if (ui->gauss_planewave->isChecked()) {
        wavefunction = gauss_plane_wave;
    } else if (ui->gauss_sine->isChecked()) {
        wavefunction = gauss_sine;
    } else {
        wavefunction = gauss_pair;
    }

    system = System(wavefunction, potential, ui->imaginaryTimeSpinBox->value());
    plot();
}

void MainWindow::on_harmonic_potential_clicked()
{
    ui->potential_a->setEnabled(true);
    ui->potential_a->setValue(80.);
    ui->potential_a->setSingleStep(5.);
    ui->potential_b->setEnabled(false);
    ui->potential_x0->setEnabled(false);
    ui->potential_x1->setEnabled(false);
    restart_system();
}

void MainWindow::on_step_potential_clicked()
{
    ui->potential_a->setEnabled(true);
    ui->potential_a->setValue(200.);
    ui->potential_a->setSingleStep(50.);
    ui->potential_b->setEnabled(false);
    ui->potential_x0->setEnabled(true);
    ui->potential_x0->setValue(3.);
    ui->potential_x1->setEnabled(true);
    ui->potential_x1->setValue(3.3);
    restart_system();
}

void MainWindow::on_free_particle_potential_clicked()
{
    ui->potential_a->setEnabled(false);
    ui->potential_b->setEnabled(false);
    ui->potential_x0->setEnabled(false);
    ui->potential_x1->setEnabled(false);
    restart_system();
}

void MainWindow::on_mexican_hat_potential_clicked()
{
    ui->potential_a->setEnabled(true);
    ui->potential_a->setSingleStep(0.1);
    ui->potential_a->setValue(1.5);
    ui->potential_b->setEnabled(true);
    ui->potential_b->setSingleStep(1.);
    ui->potential_b->setValue(25.);
    ui->potential_x0->setEnabled(false);
    ui->potential_x1->setEnabled(false);
    restart_system();
}

void MainWindow::on_gauss_clicked()
{
    ui->wavefunction_s->setEnabled(true);
    ui->wavefunction_s->setValue(8.);
    ui->wavefunction_k->setEnabled(false);
    ui->wavefunction_x0->setEnabled(true);
    ui->wavefunction_x0->setValue(0.);
    ui->wavefunction_x1->setEnabled(false);
    restart_system();
}

void MainWindow::on_gauss_sine_clicked()
{
    ui->wavefunction_s->setEnabled(true);
    ui->wavefunction_s->setValue(4.);
    ui->wavefunction_k->setEnabled(true);
    ui->wavefunction_k->setSingleStep(1.);
    ui->wavefunction_k->setValue(12.);
    ui->wavefunction_x0->setEnabled(false);
    ui->wavefunction_x1->setEnabled(false);
    restart_system();
}

void MainWindow::on_gauss_planewave_clicked()
{
    ui->wavefunction_s->setEnabled(true);
    ui->wavefunction_s->setValue(4.);
    ui->wavefunction_k->setEnabled(true);
    ui->wavefunction_k->setSingleStep(1.);
    ui->wavefunction_k->setValue(12.);
    ui->wavefunction_x0->setEnabled(false);
    ui->wavefunction_x1->setEnabled(false);
    restart_system();
}

void MainWindow::on_gauss_pair_clicked()
{
    ui->wavefunction_s->setEnabled(true);
    ui->wavefunction_s->setValue(4.);
    ui->wavefunction_k->setEnabled(false);
    ui->wavefunction_x0->setEnabled(true);
    ui->wavefunction_x0->setValue(2.);
    ui->wavefunction_x1->setEnabled(false);
    // ui->wavefunction_x1->setValue(-2.);
    restart_system();
}

void MainWindow::on_potential_a_valueChanged(double)
{
    restart_system();
}

void MainWindow::on_potential_b_valueChanged(double)
{
    restart_system();
}

void MainWindow::on_potential_x0_valueChanged(double)
{
    restart_system();
}

void MainWindow::on_potential_x1_valueChanged(double)
{
    restart_system();
}

void MainWindow::on_wavefunction_s_valueChanged(double)
{
    restart_system();
}

void MainWindow::on_wavefunction_k_valueChanged(double)
{
    restart_system();
}

void MainWindow::on_wavefunction_x0_valueChanged(double)
{
    restart_system();
}

void MainWindow::on_wavefunction_x1_valueChanged(double)
{
    restart_system();
}

void MainWindow::on_imaginaryTimeSpinBox_valueChanged(double)
{
    system.delta_t = std::complex<double>(system.re_delta_t, -system.re_delta_t*ui->imaginaryTimeSpinBox->value());
}

void MainWindow::on_linear_potential_clicked()
{
    ui->potential_a->setEnabled(true);
    ui->potential_a->setSingleStep(0.1);
    ui->potential_a->setValue(1.5);
    ui->potential_b->setEnabled(false);
    ui->potential_x0->setEnabled(false);
    ui->potential_x1->setEnabled(false);
    restart_system();
}

