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
    std::fstream pulse_file,noise_file,save_RI,save_f;
    std::vector<double> energy_p,energy_b,t0,off;


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
    s_offset/=s;
    EP.setOffset(s_offset);



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
    j=0;
    start=true;
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
    save_f << EP.getFactor() << std::endl;
    save_f.close();



    //Compute baseline-energy correlation
    pulse_file.open(pulse_path.toStdString(), std::ios::in|std::ios::binary);
    pulse_file.seekg(0, std::ios::end);
    size = pulse_file.tellg()/4;
    pulse_file.close();

    pulse_file.open(pulse_path.toStdString(), std::ios::in|std::ios::binary);
    pulse_file >> std::noskipws;
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
                        energy_b.push_back(EP.getEnergy());
                        off.push_back(EP.getOffset());
                    }
                }
            }
            j++;
        }
    }
    EP.setRecording();
    pulse_file.close();


    vector<double> AB(2),coeff(2);
    matrix<double> M(2,2);
    double m00=0,m10=0,A=0,B=0,M_det=0,mean_offset=0;
    for (int i=0;i<(int)energy_b.size();i++)
    {
        m00+=pow(off[i]/10000.0,2);
        m10+=off[i]/10000.0;
        A+=energy_b[i]*off[i]/10000.0;
        B+=energy_b[i];
        mean_offset+=off[i];
    }
    mean_offset/=(double)energy_b.size();
    M_det=(m00*((double)energy_b.size())-m10*m10);
    M(0,0)=(1.0/M_det)*(double)energy_b.size();
    M(0,1)=-(1.0/M_det)*m10;
    M(1,0)=-(1.0/M_det)*m10;
    M(1,1)=(1.0/M_det)*m00;
    AB(0)=A;
    AB(1)=B;
    coeff=prod(M,AB);

    save_f.open("Factor.txt",std::ios::out|std::ios::app);
    save_f << mean_offset << std::endl;
    for (int i=0;i<2;i++)
    {
        save_f << coeff(i) << std::endl;
    }
    save_f.close();

    EP.setB_coeff(coeff);
    EP.setOffset(mean_offset);



    //E=f(t0)
    pulse_file.open(pulse_path.toStdString(), std::ios::in|std::ios::binary);
    pulse_file >> std::noskipws;
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
                        energy_p.push_back(EP.getEnergy());
                        t0.push_back(EP.gett0());
                    }
                }
            }
            j++;
        }
    }
    pulse_file.close();

    vector<double> E((int)energy_p.size()),p_coeff(3);
    matrix<double> T((int)energy_p.size(),3),Tinv(3,3),Tin(3,3);
    double mean_phase=0;
    for (int i=0;i<(int)energy_p.size();i++)
    {
        mean_phase+=t0[i];
        E(i)=energy_p[i];
        for (int j=0;j<3;j++)
        {
            T(i,j)=pow(t0[i],2-j);
        }
    }
    mean_phase/=(double)t0.size();

    Tin=prod(trans(T),T);
    EP.InvertMatrix(Tin,Tinv);
    p_coeff=prod(prod(Tinv,trans(T)),trans(E));
    save_f.open("Factor.txt",std::ios::out|std::ios::app);
    for (int i=0;i<3;i++)
    {
        save_f << p_coeff(i) << std::endl;
    }
    save_f << mean_phase << std::endl;
    save_f.close();

    EP.setP_coeff(p_coeff);
    EP.setPhase(mean_phase);

}
