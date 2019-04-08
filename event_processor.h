#ifndef EVENT_PROCESSOR_H
#define EVENT_PROCESSOR_H
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <valarray>
#include <QString>

using namespace boost::numeric::ublas;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

class Event_Processor
{
    public:
        Event_Processor();
        Event_Processor(int Npattern);
        virtual ~Event_Processor();
        void trigger_function();
        void store_noise();
        void computeOptimalFilter();
        void computeFit();
        void computeEventProcessor();
        double getEnergy();
        template<class T> bool InvertMatrix(const matrix<T>& input, matrix<T>& inverse);
        void setInput(double input);
        void setFactor(double f);
        void setThreshold(double thres);
        void setImpulseResponse(vector<double> IR);
        bool getRecording();
        void computeImpulseResponse();
        void recordImpulseResponse();
        void ExtractData(QString path, char *buffer, int &ind, int &size);
        double getData(char buffer[]);
        void fft(CArray& x);
        void ifft(CArray& x);
        void setMode(bool mod);
        CArray getIR();
        void setRecording();
        void recordFactor();
        void computeFactor();
        double getFactor();
        double gett0();
        double getOffset();
        void setOffset(double off);
        void setCorr_coeff(vector<double> v);
        double computeMean();

    private:
        int index;
        int counter;
        int count;
        bool wait;
        double Threshold;
        bool ReadyToCompute;
        bool recording;
        vector<double> Trigger_coeff;
        vector<double> Buffer;
        vector<double> corr_coeff;
        vector<double> Record;
        vector<double> OutputFilter;
        vector<double> ImpulseResponse;
        vector<double> Poly_coeff;
        double Trigger_output;
        int RecordSize;
        double factor;
        double energy;
        double offset;
        double t0;
        matrix<double> Z;
        int noise_size;
        int pulse_size;
        int noise_offset;
        int pulse_offset;
        char *noise;
        char *pulse;
        bool mode;
        int factor_count;
        CArray pulse_fft;
        CArray noise_fft;
        CArray pulse_phase;
        CArray IR;
};

#endif // EVENT_PROCESSOR_H
