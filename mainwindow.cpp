#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <iostream>
#include <QThread>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    QObject::connect(this,SIGNAL(setname()), this, SLOT(set_noise_filename()));
    QObject::connect(this,SIGNAL(setname()), this, SLOT(set_pulse_filename()));
    QObject::connect(ui->Generate,SIGNAL(clicked()), this, SLOT(compute()));
    proc = new Processing(1024,40);
    setWindowTitle("Calibration");
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::open_noise()
{
    noise_filename = QFileDialog::getOpenFileName(this, "Ouvrir un fichier", QString());
    emit setname();
}

void MainWindow::open_pulse()
{
    pulse_filename = QFileDialog::getOpenFileName(this, "Ouvrir un fichier", QString());
    emit setname();
}

void MainWindow::set_noise_filename()
{
    ui->noise_file->clear();
    ui->noise_file->append(noise_filename);
}

void MainWindow::set_pulse_filename()
{
    ui->pulse_file->clear();
    ui->pulse_file->append(pulse_filename);
}

void MainWindow::compute()
{
    /*QThread *thread = new QThread();
    proc.moveToThread(QApplication::instance()->thread());
    connect(thread, SIGNAL(finished()), &proc, SLOT(deleteLater()));
    connect(thread, SIGNAL(started()), &proc, SLOT(calibrate()));
    thread->start();*/
    proc->calibrate(pulse_filename,noise_filename);
}
