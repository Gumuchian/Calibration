#ifndef PROCESSING_H
#define PROCESSING_H
#include "event_processor.h"

using namespace boost::numeric::ublas;

class Processing
{
public:
    Processing();
    Processing(int Nbre, double threshold);
    void calibrate(QString pulse_path, QString noise_path);
private:
    int N;
    Event_Processor EP;
    int thres;
};

#endif // PROCESSING_H
