#include "event_processor.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
#include <iostream>
#include <fstream>
#define PI 3.14159265358979


Event_Processor::Event_Processor()
{

}

Event_Processor::Event_Processor(int Npattern):Trigger_coeff(8,0),Buffer(8,0),corr_coeff(3,0),Record(Npattern+2,0),OutputFilter(3,0),ImpulseResponse(Npattern,0),Z(3,3),pulse_fft(Npattern),noise_fft(Npattern),pulse_phase(Npattern),IR(Npattern)
{
    counter = 0;
    count = 0;
    wait = false;
    factor_count = 0;
    factor = 0;
    energy = 0;
    RecordSize = Npattern;
    recording = false;
    ReadyToCompute = false;
    corr_coeff(0)=0;
    corr_coeff(1)=0;
    corr_coeff(2)=7000;
    for (int k=0;k<8;k++)
    {
        if (k<4)
        {
            Trigger_coeff(k) = -1;
        }
        else
        {
            Trigger_coeff(k) = 1;
        }
    }
    matrix<double> X(3,3);
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            X(i,j)=pow(i-3/2,2-j);
        }
    }
    InvertMatrix(X,Z);
    for (int k=0;k<Npattern;k++)
    {
        pulse_fft[k]=0;
        noise_fft[k]=0;
    }
}

Event_Processor::~Event_Processor()
{

}

void Event_Processor::trigger_function()
{
    //Trigger_output = inner_prod(Buffer,Trigger_coeff);
    Trigger_output = offset - Buffer(0);
    if (!recording && ReadyToCompute)
    {
        ReadyToCompute = false;
    }
    if (std::abs(Trigger_output) > Threshold && count == 0)
    {
        wait = true;
    }
    if (wait && std::abs(Trigger_output) < Threshold)
    {
        wait = false;
    }
    count ++;
    if (std::abs(Trigger_output) > Threshold && !wait)
    {
        recording = true;
    }
    if (recording)
    {
        Record(counter) = offset - Buffer(2);
        counter++;
        if (counter == RecordSize+2)
        {
            recording = false;
            ReadyToCompute = true;
            counter = 0;
            count = 0;
        }
    }
}

void Event_Processor::computeOptimalFilter()
{
    for(int i=0;i<3;i++)
    {
        OutputFilter(i) = inner_prod(subslice(Record,i,1,(int)(Record.size()-2)),ImpulseResponse);
    }
}

void Event_Processor::computeFit()
{
    Poly_coeff=prod(Z,OutputFilter);
    energy=7000.0*(Poly_coeff(2)-pow(Poly_coeff(1),2)/(2*Poly_coeff(0)))/factor;
    t0=-Poly_coeff(1)/(2*Poly_coeff(0));
    energy/=(corr_coeff(0)*pow(t0,2)+corr_coeff(1)*t0+corr_coeff(2))/7000;
}

void Event_Processor::setInput(double input)
{
    for (int i=7;i>0;i--)
    {
        Buffer(i) = Buffer(i-1);
    }
    Buffer(0) = input;
}

double Event_Processor::getEnergy()
{
    return energy;
}

void Event_Processor::computeEventProcessor()
{
    trigger_function();
    if (ReadyToCompute)
    {
        computeOptimalFilter();
        computeFit();
    }
}

void Event_Processor::setThreshold(double thres)
{
    Threshold = thres;
}

void Event_Processor::setFactor(double f)
{
    factor = f;
}

bool Event_Processor::getRecording()
{
    return ReadyToCompute;
}

void Event_Processor::setImpulseResponse(vector<double> IR)
{
    ImpulseResponse = IR;
}

void Event_Processor::recordImpulseResponse()
{
    trigger_function();
    if (ReadyToCompute)
    {
        CArray module(RecordSize);
        Complex c[RecordSize];
        for (int i=0;i<RecordSize;i++)
        {
            c[i]=Record(i+1);
            module[i]=c[i];
        }
        fft(module);
        if (mode)
        {
            pulse_fft+=abs(module);
            for (int j=0;j<RecordSize;j++)
            {
                pulse_phase[j]=arg(module[j]);
            }
        }
        else
        {
            noise_fft+=abs(module);
        }
    }
}

void Event_Processor::computeImpulseResponse()
{
    IR = pulse_fft/noise_fft;
    const Complex const_i(0,1);
    for (int i=0;i<RecordSize;i++)
    {
        IR[i]*=std::exp(const_i*pulse_phase[i]);
    }
    ifft(IR);
    for (int i=0;i<RecordSize;i++)
    {
        IR[i]=std::real(IR[i]);
    }
}

void Event_Processor::recordFactor()
{
    trigger_function();
    if (ReadyToCompute)
    {
        factor_count+=1;
        computeOptimalFilter();
        factor+=OutputFilter(1);
    }
}

void Event_Processor::computeFactor()
{
    factor/=factor_count;
    factor_count = 0;
}

template<class T> bool Event_Processor::InvertMatrix(const matrix<T>& input, matrix<T>& inverse)
{
    typedef permutation_matrix<std::size_t> pmatrix;
    matrix<T> A(input);
    pmatrix pm(A.size1());
    int res = lu_factorize(A, pm);
    if (res != 0) return false;
    inverse.assign(identity_matrix<T> (A.size1()));
    lu_substitute(A, pm, inverse);
    return true;
}

void Event_Processor::fft(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;

    CArray even = x[std::slice(0,N/2,2)];
    CArray  odd = x[std::slice(1,N/2,2)];

    fft(even);
    fft(odd);

    for (size_t k=0;k<N/2;++k)
    {
        Complex t = std::polar(1.0,-2*PI*k/N)*odd[k];
        x[k]=even[k]+t;
        x[k+N/2]=even[k]-t;
    }
}

void Event_Processor::ifft(CArray& x)
{
    x = x.apply(std::conj);
    fft(x);
    x = x.apply(std::conj);
    x /= x.size();
}

void Event_Processor::setMode(bool mod)
{
    mode = mod;
}

CArray Event_Processor::getIR()
{
    return IR;
}

void Event_Processor::setRecording()
{
    ReadyToCompute = false;
    recording = false;
}

double Event_Processor::getFactor()
{
    return factor;
}

double Event_Processor::gett0()
{
    return t0;
}

void Event_Processor::setOffset(double off)
{
    offset=off;
}

void Event_Processor::setCorr_coeff(vector<double> v)
{
    for (int i=0;i<3;i++)
    {
        corr_coeff(i)=v(i);
    }
}

double Event_Processor::getData(char buffer[])
{
    int I,Q;
    I=(unsigned char)buffer[1]*256 + (unsigned char)buffer[0];
    Q=(unsigned char)buffer[3]*256 + (unsigned char)buffer[2];
    if (I>31767)
    {
        I=I-65536;
    }
    if (Q>31767)
    {
        Q=Q-65536;
    }
    return sqrt(pow(I,2)+pow(Q,2));
}
