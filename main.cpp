#include "mainwindow.h"
#include <QApplication>
#include "processing.h"
#include <fstream>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    //MainWindow w;
    Processing proc(2048,40);
    QStringList args;
    args=QCoreApplication::arguments();
    proc.calibrate(args[1],args[2]);
    //w.show();
    a.quit();
    //return a.exec();
    return 0;
}
