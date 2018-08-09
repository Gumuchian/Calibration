#include "mainwindow.h"
#include <QApplication>
#include "processing.h"
#include <fstream>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    //MainWindow w;
    Processing proc(2048,40);
    //QStringList args = QCoreApplication::arguments();
    proc.calibrate("C:/Users/Paul/Desktop/test_Laurent/20180713_174251_0000_IQ_ALL_DUMP_Channel_0_Energy_Calibration@7keV.dat","C:/Users/Paul/Desktop/test_Laurent/20180713_174246_0000_IQ_ALL_DUMP_Channel_0_Noise_Calibration.dat");
    //w.show();
    a.quit();
    std::cout << "done" << std::endl;
    //return a.exec();
    return 0;
}
