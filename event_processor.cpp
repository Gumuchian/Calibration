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

Event_Processor::Event_Processor(int Npattern):Trigger_coeff(8,0),Buffer(200,0),p_coeff(3,0),b_coeff(2,0),Record(Npattern+2,0),OutputFilter(3,0),ImpulseResponse(Npattern,0),Z(3,3),pulse_fft(Npattern),noise_fft(Npattern),pulse_phase(Npattern),IR(Npattern)
{
    index = 0;
    counter = 0;
    count = 0;
    wait = false;
    factor_count = 0;
    factor = 0;
    energy = 0;
    RecordSize = Npattern;
    recording = false;
    ReadyToCompute = false;
    b_coeff(0)=0;
    b_coeff(1)=7000;
    p_coeff(0)=0;
    p_coeff(1)=0;
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
    matrix<double> X(3,3),Xp(3,3);
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            X(i,j)=pow(i-3/2,2-j);
        }
    }
    Xp=prod(X,trans(X));
    InvertMatrix(Xp,Z);
    Z=prod(trans(X),Z);
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
    if (!recording)
    {
        offset=computeMean();
    }
    Trigger_output = - Buffer(index);
    if (!recording && ReadyToCompute)
    {
        ReadyToCompute = false;
    }
    if (std::abs(offset+Trigger_output) > Threshold && count == 0)
    {
        wait = true;
    }
    if (wait && std::abs(offset+Trigger_output) < Threshold)
    {
        wait = false;
    }
    count ++;
    if (std::abs(offset+Trigger_output) > Threshold && !wait)
    {
        recording = true;
    }
    if (recording)
    {
        Record(counter) = - Buffer((index+1)%200);
        counter++;
        if (counter == RecordSize+2)
        {
            recording = false;
            ReadyToCompute = true;
            counter = 0;
            count = 0;
        }
    }
    index ++;
    index = index%200;
}

void Event_Processor::store_noise()
{
    Trigger_output = Buffer(index);
    if (!recording && ReadyToCompute)
    {
        ReadyToCompute = false;
    }
    if (!recording && !ReadyToCompute)
    {
        recording = true;
    }
    if (recording)
    {
        Record(counter) = Buffer((index+1)%200);
        counter++;
        if (counter == RecordSize+2)
        {
            recording = false;
            ReadyToCompute = true;
            counter = 0;
        }
    }
    index ++;
    index = index%200;
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
    energy=7000.0*(Poly_coeff(2)-pow(Poly_coeff(1),2)/(4*Poly_coeff(0)))/factor;
    t0=-Poly_coeff(1)/(2*Poly_coeff(0));
    energy-=p_coeff(0)*pow(t0,2)+p_coeff(1)*t0-(p_coeff(0)*pow(t0_mean,2)+p_coeff(1)*t0_mean);
    energy-=b_coeff(0)*(offset-mean_offset)/10000.0;
}

void Event_Processor::setInput(double input)
{
    Buffer(index) = input;
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
    CArray module(RecordSize);
    Complex c[RecordSize];
    for (int i=0;i<RecordSize;i++)
    {
        c[i]=Record(i+1);
        module[i]=c[i];
    }
    if (mode)
    {
        trigger_function();
        if (ReadyToCompute)
        {
            pulse_fft+=module;
        }
    }
    else
    {
        store_noise();
        if (ReadyToCompute)
        {
            fft(module);
            noise_fft+=pow(abs(module),2);
        }
    }
}

void Event_Processor::computeImpulseResponse()
{
    fft(pulse_fft);
    pulse_fft[0]=0;
    IR = pulse_fft/noise_fft;

    std::fstream file_noise,file_pulse;
    file_noise.open("Noise_spectrum.txt",std::ios::out);
    file_pulse.open("Pulse_spectrum.txt",std::ios::out);
    ifft(IR);
    for (int i=0;i<RecordSize;i++)
    {
        IR[i]=real(IR[i]);
        file_noise << abs(noise_fft[i]) << std::endl;
        file_pulse << abs(pulse_fft[i]) << std::endl;
    }
    file_pulse.close();
    file_noise.close();
}

void Event_Processor::recordFactor()
{
    trigger_function();
    if (ReadyToCompute)
    {
        factor_count++;
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

double Event_Processor::getOffset()
{
    return offset;
}

void Event_Processor::setOffset(double off)
{
    mean_offset=off;
}

void Event_Processor::setPhase(double phase)
{
    t0_mean=phase;
}

void Event_Processor::setP_coeff(vector<double> v)
{
    for (int i=0;i<(int)p_coeff.size();i++)
    {
        p_coeff(i)=v(i);
    }
}

void Event_Processor::setB_coeff(vector<double> v)
{
    for (int i=0;i<(int)b_coeff.size();i++)
    {
        b_coeff(i)=v(i);
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

double Event_Processor::computeMean()
{
    double sum=0;
    for (int i=0;i<170;i++)
    {
        sum+=Buffer((index+1+i)%200);
    }
    return sum/170.0;
}
