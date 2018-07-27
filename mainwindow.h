#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QString>
#include "processing.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void open_noise();
    void open_pulse();
    void set_noise_filename();
    void set_pulse_filename();
    void compute();

signals:
    void setname();

private:
    Ui::MainWindow *ui;
    QString noise_filename;
    QString pulse_filename;
    Processing *proc;
};

#endif // MAINWINDOW_H
