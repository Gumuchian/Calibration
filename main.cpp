#include "mainwindow.h"
#include <QApplication>
#include "processing.h"
#include <fstream>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    //MainWindow w;
    QStringList args = QCoreApplication::arguments();
    Processing proc(args[3].toInt(),args[4].toInt());
    proc.calibrate(args[1],args[2]);
    //w.show();
    a.quit();
    //return a.exec();
    return 0;
}
