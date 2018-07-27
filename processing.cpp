#include "processing.h"
#include "event_processor.h"
#include <iostream>
#include <fstream>
#include <QString>
#include <iomanip>

Processing::Processing()
{

}

Processing::Processing(int Nbre, double threshold):N(Nbre),EP(Nbre),thres(threshold)
{

}

void Processing::calibrate(QString pulse_path, QString noise_path)
{
    std::string str;
    std::fstream pulse_file,noise_file,save_RI,save_f,test,save_noise;
    std::vector<double> energy,t0;

    //Set offset
    char IQ[4];
    int j=0,s=0,Npix=41,pas=(Npix+2)*2;
    double s_offset=0;
    bool start=true;
    noise_file.open(noise_path.toStdString(), std::ios::in|std::ios::binary);
    noise_file.seekg(0, std::ios::end);
    int size = noise_file.tellg();
    noise_file.close();

    noise_file.open(noise_path.toStdString(), std::ios::in|std::ios::binary);
    noise_file >> std::noskipws;

    while(j<size/4)
    {
        if((int)(((unsigned char)IQ[1]*256 + (unsigned char)IQ[0]))!=56026 && start)
        {
            noise_file >> IQ[0];
            noise_file >> IQ[1];
        }
        else
        {
            start=false;
            if (j==0)
            {
                noise_file >> IQ[0];
                noise_file >> IQ[1];
                for (int i=0;i<Npix+1;i++)
                {
                    noise_file >> IQ[0];
                    noise_file >> IQ[1];
                    noise_file >> IQ[2];
                    noise_file >> IQ[3];
                }
            }
            else
            {
                for (int i=0;i<4;i++)
                {
                    noise_file >> IQ[i];
                }

                if (j%pas==0)
                {
                    s_offset+=EP.getData(IQ);
                    s++;
                }
            }
            j++;
        }
    }
    noise_file.close();
    EP.setOffset(s_offset/s);


    // Record pulse for IR
    pulse_file.open(pulse_path.toStdString(), std::ios::in|std::ios::binary);
    pulse_file.seekg(0, std::ios::end);
    size = pulse_file.tellg()/4;
    pulse_file.close();

    pulse_file.open(pulse_path.toStdString(), std::ios::in|std::ios::binary);
    pulse_file >> std::noskipws;

    EP.setMode(true);
    EP.setThreshold(thres);
    j=0;
    start=true;
    while(j<size)
    {
        if((int)(((unsigned char)IQ[1]*256 + (unsigned char)IQ[0]))!=56026 && start)
        {
            pulse_file >> IQ[0];
            pulse_file >> IQ[1];
        }
        else
        {
            start=false;
            if (j==0)
            {
                pulse_file >> IQ[0];
                pulse_file >> IQ[1];
                for (int i=0;i<Npix+1;i++)
                {
                    pulse_file >> IQ[0];
                    pulse_file >> IQ[1];
                    pulse_file >> IQ[2];
                    pulse_file >> IQ[3];
                }
            }
            else
            {
                for (int i=0;i<4;i++)
                {
                    pulse_file >> IQ[i];
                }

                if (j%pas==0)
                {
                    EP.setInput(EP.getData(IQ));
                    EP.recordImpulseResponse();
                }
            }
            j++;
        }
    }
    EP.setRecording();
    pulse_file.close();


    // Record noise for IR
    noise_file.open(noise_path.toStdString(), std::ios::in|std::ios::binary);
    noise_file.seekg(0, std::ios::end);
    size = noise_file.tellg()/4;
    noise_file.close();

    noise_file.open(noise_path.toStdString(), std::ios::in|std::ios::binary);
    noise_file >> std::noskipws;

    EP.setMode(false);
    EP.setThreshold(0);
    j=0;
    while(j<size)
    {
        if((int)(((unsigned char)IQ[1]*256 + (unsigned char)IQ[0]))!=56026 && start)
        {
            noise_file >> IQ[0];
            noise_file >> IQ[1];
        }
        else
        {
            start=false;
            if (j==0)
            {
                noise_file >> IQ[0];
                noise_file >> IQ[1];
                for (int i=0;i<Npix+1;i++)
                {
                    noise_file >> IQ[0];
                    noise_file >> IQ[1];
                    noise_file >> IQ[2];
                    noise_file >> IQ[3];
                }
            }
            else
            {
                for (int i=0;i<4;i++)
                {
                    noise_file >> IQ[i];
                }

                if (j%pas==0)
                {
                    EP.setInput(EP.getData(IQ));
                    EP.recordImpulseResponse();
                }
            }
            j++;
        }
    }
    EP.setRecording();
    noise_file.close();
    //Compute IR
    EP.computeImpulseResponse();


    CArray OF(N);
    vector<double> IR(N);
    OF = EP.getIR();
    save_RI.open("Pattern.txt",std::ios::out);
    for (int i=0;i<N;i++)
    {
        save_RI << real(OF[i]) << std::endl;
        IR(i)=real(OF[i]);
    }
    save_RI.close();
    EP.setImpulseResponse(IR);
    EP.setRecording();


    //Record factor
    pulse_file.open(pulse_path.toStdString(), std::ios::in|std::ios::binary);
    pulse_file.seekg(0, std::ios::end);
    size = pulse_file.tellg()/4;
    pulse_file.close();

    pulse_file.open(pulse_path.toStdString(), std::ios::in|std::ios::binary);
    pulse_file >> std::noskipws;
    EP.setThreshold(thres);
    j=0;
    start=true;
    while(j<size)
    {
        if((int)(((unsigned char)IQ[1]*256 + (unsigned char)IQ[0]))!=56026 && start)
        {
            pulse_file >> IQ[0];
            pulse_file >> IQ[1];
        }
        else
        {
            start=false;
            if (j==0)
            {
                pulse_file >> IQ[0];
                pulse_file >> IQ[1];
                for (int i=0;i<Npix+1;i++)
                {
                    pulse_file >> IQ[0];
                    pulse_file >> IQ[1];
                    pulse_file >> IQ[2];
                    pulse_file >> IQ[3];
                }
            }
            else
            {
                for (int i=0;i<4;i++)
                {
                    pulse_file >> IQ[i];
                }

                if (j%pas==0)
                {
                    EP.setInput(EP.getData(IQ));
                    EP.recordFactor();
                }
            }
            j++;
        }
    }
    EP.setRecording();
    pulse_file.close();


    //Compute factor
    EP.computeFactor();
    save_f.open("Factor.txt",std::ios::out);
    save_f << EP.getFactor();
    save_f.close();

    //E=f(t0)
    pulse_file.open(pulse_path.toStdString(), std::ios::in|std::ios::binary);
    pulse_file >> std::noskipws;
    EP.setThreshold(thres);
    j=0;
    start=true;
    while(j<size)
    {
        if((int)(((unsigned char)IQ[1]*256 + (unsigned char)IQ[0]))!=56026 && start)
        {
            pulse_file >> IQ[0];
            pulse_file >> IQ[1];
        }
        else
        {
            start=false;
            if (j==0)
            {
                pulse_file >> IQ[0];
                pulse_file >> IQ[1];
                for (int i=0;i<Npix+1;i++)
                {
                    pulse_file >> IQ[0];
                    pulse_file >> IQ[1];
                    pulse_file >> IQ[2];
                    pulse_file >> IQ[3];
                }
            }
            else
            {
                for (int i=0;i<4;i++)
                {
                    pulse_file >> IQ[i];
                }

                if (j%pas==0)
                {
                    EP.setInput(EP.getData(IQ));
                    EP.computeEventProcessor();
                    if (EP.getRecording())
                    {
                        energy.push_back(EP.getEnergy());
                        t0.push_back(EP.gett0());
                    }
                }
            }
            j++;
        }
    }
    pulse_file.close();
    test.open("Test.txt",std::ios::out);
    vector<double> E((int)energy.size()),coeff(3);
    matrix<double> T((int)energy.size(),3),Tinv(3,3),Tin(3,3);
    test << std::setprecision(10);
    for (int i=0;i<(int)energy.size();i++)
    {
        E(i)=energy[i];
        test << energy[i] << "\t" << t0[i] << std::endl;
        for (int j=0;j<3;j++)
        {
            T(i,j)=pow(t0[i],2-j);
        }
    }
    test.close();
    Tin=prod(trans(T),T);
    EP.InvertMatrix(Tin,Tinv);
    coeff=prod(prod(Tinv,trans(T)),trans(E));
}
