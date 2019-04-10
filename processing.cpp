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
    std::fstream pulse,noise;
    pulse.open("Pulse.txt",std::ios::in);
    noise.open("Noise.txt",std::ios::in);

    std::string str;
    std::fstream pulse_file,noise_file,save_RI,save_f;
    std::vector<double> energy,t0,off;

    //Set offset
    char IQ[4];
    int j=0,s=0,Npix=41,pas=(Npix+2)*2;
    double s_offset=0,var_noise=0;
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
                    double noi;
                    noise >> noi;
                    s_offset+=noi;
                    //s_offset+=EP.getData(IQ);
                    s++;
                }
            }
            j++;
        }
    }

    noise.close();

    noise_file.close();
    s_offset/=s;
    EP.setOffset(s_offset);

    //Get noise threshold
    noise.open("Noise.txt",std::ios::in);
    noise_file.open(noise_path.toStdString(), std::ios::in|std::ios::binary);
    noise_file.seekg(0, std::ios::end);
    size = noise_file.tellg();
    noise_file.close();

    noise_file.open(noise_path.toStdString(), std::ios::in|std::ios::binary);
    noise_file >> std::noskipws;
    j=0;
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
                    double noi;
                    noise >> noi;
                    var_noise+=pow(noi-s_offset,2);
                    //s_offset+=EP.getData(IQ);
                }
            }
            j++;
        }
    }

    var_noise=sqrt(var_noise/s);
    noise.close();
    noise_file.close();

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
                    double in;
                    pulse >> in;
                    EP.setInput(in);
                    //EP.setInput(EP.getData(IQ));
                    EP.recordImpulseResponse();
                }
            }
            j++;
        }
    }
    EP.setRecording();
    pulse_file.close();

    pulse.close();
    noise.open("Noise.txt",std::ios::in);

    // Record noise for IR
    noise_file.open(noise_path.toStdString(), std::ios::in|std::ios::binary);
    noise_file.seekg(0, std::ios::end);
    size = noise_file.tellg()/4;
    noise_file.close();

    noise_file.open(noise_path.toStdString(), std::ios::in|std::ios::binary);
    noise_file >> std::noskipws;

    EP.setMode(false);
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
                    double in;
                    noise >> in;
                    EP.setInput(in);
                    //EP.setInput(EP.getData(IQ));
                    EP.recordImpulseResponse();
                }
            }
            j++;
        }
    }
    EP.setRecording();
    noise_file.close();

    noise.close();
    pulse.close();

    //Compute IR
    EP.computeImpulseResponse();


    CArray OF(N);
    vector<double> IR(N);
    OF = EP.getIR();
    save_RI.open("Pattern.txt",std::ios::out);
    std::fstream fff;
    fff.open("Pat.txt",std::ios::in);
    for (int i=0;i<N;i++)
    {
        save_RI << real(OF[i]) << std::endl;
        IR(i)=real(OF[i]);
        double in;
        fff >> in;
        IR(i)=in;
    }
    save_RI.close();
    EP.setImpulseResponse(IR);
    EP.setRecording();


    //Record factor
    pulse.open("Pulse.txt",std::ios::in);

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
                    double in;
                    pulse >> in;
                    EP.setInput(in);
                    //EP.setInput(EP.getData(IQ));
                    EP.recordFactor();
                }
            }
            j++;
        }
    }
    EP.setRecording();
    pulse_file.close();
    pulse.close();


    //Compute factor
    EP.computeFactor();
    save_f.open("Factor.txt",std::ios::out);
    save_f << EP.getFactor() << std::endl << s_offset << std::endl;
    save_f.close();


    //Compute Offset/Energy correlation
    pulse.open("Pulse.txt",std::ios::in);
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
                    double in;
                    pulse >> in;
                    EP.setInput(in);
                    EP.computeEventProcessor();
                    if (EP.getRecording())
                    {
                        energy.push_back(EP.getEnergy());
                        off.push_back(EP.getOffset());
                    }
                }
            }
            j++;
        }
    }
    EP.setRecording();
    pulse_file.close();
    pulse.close();


    std::fstream file_o;
    file_o.open("offset_energy.txt",std::ios::out);
    for (int i=0;i<(int)energy.size();i++)
    {
        file_o << energy[i] << "\t" << off[i] << std::endl;
    }
    file_o.close();

    vector<double> AB(2),coeff(2);
    matrix<double> M(2,2);
    double m00=0,m10=0,A=0,B=0,M_det=0,mean_offset=0;
    for (int i=0;i<(int)energy.size();i++)
    {
        m00+=pow(off[i]/10000.0,2);
        m10+=off[i]/10000.0;
        A+=energy[i]*off[i]/10000.0;
        B+=energy[i];
        mean_offset+=off[i];
    }
    mean_offset/=(double)energy.size();
    M_det=(m00*((double)energy.size())-m10*m10);
    M(0,0)=(1.0/M_det)*(double)energy.size();
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

    EP.setCorr_coeff(coeff);
    EP.setOffset(mean_offset);
   /*
    //E=f(t0)
    pulse.open("Pulse.txt",std::ios::in);

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
                    double in;
                    pulse >> in;
                    EP.setInput(in);
                    //EP.setInput(EP.getData(IQ));
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

    vector<double> E((int)energy.size()),coeff(3);
    matrix<double> T((int)energy.size(),3),Tinv(3,3),Tin(3,3);
    for (int i=0;i<(int)energy.size();i++)
    {
        E(i)=energy[i];
        for (int j=0;j<3;j++)
        {
            T(i,j)=pow(t0[i],2-j);
        }
    }

    Tin=prod(trans(T),T);
    EP.InvertMatrix(Tin,Tinv);
    coeff=prod(prod(Tinv,trans(T)),trans(E));
    save_f.open("Factor.txt",std::ios::out|std::ios::app);
    for (int i=0;i<3;i++)
    {
        save_f << coeff(i) << std::endl;
    }
    save_f.close();
    */
    std::fstream file,file_s;
    std::vector<double> ener,offf;
    EP.setRecording();
    file.open("Measure.txt",std::ios::in);
    for (int i=0;i<5220914;i++)
    {
        double in;
        file >> in;
        EP.setInput(in);
        EP.computeEventProcessor();
        if (EP.getRecording())
        {
            ener.push_back(EP.getEnergy());
            offf.push_back(EP.getOffset());
        }
    }
    file.close();
    file_s.open("List_energies.txt",std::ios::out);
    for (int i=0;i<(int)ener.size();i++)
    {
        file_s << ener[i] << "\t" << offf[i] << std::endl;
    }
    file_s.close();
}
