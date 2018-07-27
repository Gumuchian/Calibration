#include "mainwindow.h"
#include <QApplication>
#include "processing.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    //MainWindow w;
    Processing proc(2048,40);
    QStringList args;
    args=QCoreApplication::arguments();
    QString k=QString::fromStdString("C:/Users/Paul/Desktop/test_Laurent/20180713_174246_0000_IQ_ALL_DUMP_Channel_0_Noise_Calibration.dat");
    QString kk=QString::fromStdString("C:/Users/Paul/Desktop/test_Laurent/20180713_174251_0000_IQ_ALL_DUMP_Channel_0_Energy_Calibration@7keV.dat");
    QString q=args[1];
    QString qq=args[2];
    std::cout << q.toStdString() << std::endl;
    std::cout << qq.toStdString() << std::endl;
    proc.calibrate(kk,k);
    //w.show();

    return a.exec();
}
