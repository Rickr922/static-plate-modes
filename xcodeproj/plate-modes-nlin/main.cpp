//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Modal Plate Reverb
// C++ implementation by: Riccardo Russo, University of Bologna
//
// Date: 11-July-2023
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include "plateModalData.h"
#include "../../nemus-cpp-audio-tools/src/audio.h"
#include <chrono>

#define PI 3.141592653589793

int main(int argc, const char * argv[])
{
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    // Excitation Signal

#define osFac 1
#define durSec 2
#define baseSR 48000
    const unsigned int timeSamples = baseSR * osFac * durSec;
    float excit[timeSamples];

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    float SR = osFac * baseSR;
    float k = 1.0 / SR;
    
    float a = 30;
    float gain = 10000.f;
    
    std::cout << "Number of Modes: " << modesNumber << '\n';
    
    //// Simulation
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    /// Initializing vectors

    float output[timeSamples];

    float q[modesNumber] {0};
    float p[modesNumber] {0};

    for (int i = 0; i < timeSamples; ++i)
    {
        excit[i] = 0.f;
        output[i] = 0.f;
    }
    
    for (int i = 0; i < modesNumber; ++i)
    {
        dampCoeffs[i] = 0.5f * k * dampCoeffs[i] / a;
        eigenFreqs[i] = 0.5f * k * eigenFreqs[i];
    }

    excit[0] = 1.f;
    
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();
    
    for (int n = 0;  n < timeSamples; ++n)
    {
        double exc = excit[n] * gain;
        
        for (int m = 0 ; m < modesNumber; ++m)
        {
            float& omega_m =  eigenFreqs[m];
            float& sigma_m =  dampCoeffs[m];
            float& q_m =  q[m];
            float& p_m =  p[m];
            float gi_m = modesIn[m] * k;
            
            float eta = (a * p_m);
            float d = 1.0 + eta*eta;
            float lambda = 1.0 + 3.0*eta*eta;
            
            float detCoeff = 1.0 / (1.0 + (sigma_m * lambda * a) + (omega_m * omega_m));
                    float b1 = q_m + omega_m * p_m;
                    float b2 = -omega_m * q_m + p_m + sigma_m * eta * (lambda - 2 * d) + gi_m * exc;

                    float qn = detCoeff * (
                        (1 + sigma_m * a * lambda) * b1
                        + omega_m * b2
                        );
                    float pn = detCoeff * (b2 - omega_m * b1);
            
            q_m = qn;
            p_m = pn;
            
            output[n] += q_m * modesOut[m];
        }
        
    }
    
    std::chrono::steady_clock::time_point toc  = std::chrono::steady_clock::now();
    double computeTime = std::chrono::duration<double >(std::chrono::duration_cast<std::chrono::nanoseconds>(toc - tic)).count();
    std::cout << "Compute Time: " << computeTime << '\n';
    
    //Debugging
    std::ofstream outFile;
    outFile.open("output.txt");

    if (!outFile.fail())
    {
        for (int n = 0; n < timeSamples; ++n)
        {
            outFile << output[n] << '\n';
        }

        outFile.close();
    }

    /*float* outPtr = &output[0];
    float* downsampled;
    if (osFac > 1)
    {
        downsampled = resample(outPtr, (double)SR, (double)baseSR, (int)timeSamples);
    }
    else
    {
        downsampled = output;
    }*/

    //audiowrite("plate-saturator-class", downsampled, timeSamples, baseSR);
    audiowrite("nlin_plate", output, timeSamples, baseSR);

    return 0;
}
