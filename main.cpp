//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Modal Plate Reverb
// C++ implementation by: Riccardo Russo, University of Bologna
//
// Date: 11-July-2023
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include "plateModalData.h"
#include "nemus-cpp-audio-tools/src/audio.h"
#include <chrono>

#define PI 3.141592653589793

int main(int argc, const char * argv[])
{
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    // Excitation Signal

#define osFac 1
#define durSec 2
#define baseSR 44100
    const unsigned int timeSamples = baseSR * osFac * durSec;
    float excit[timeSamples];

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    double SR = osFac * baseSR;
    double k = 1.0 / SR;
    
    std::cout << "Number of Modes: " << modesNumber << '\n';
    
    //// Simulation
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    /// Initializing vectors

    float output[timeSamples];

    float x[modesNumber];
    float xPrev[modesNumber];
    float xNext[modesNumber];

    for (int i = 0; i < timeSamples; ++i)
    {
        excit[i] = 0.f;
        output[i] = 0.f;

        if (i < modesNumber)
        {
            x[i] = 0.f;
            xPrev[i] = 0.f;
            xNext[i] = 0.f;
        }
    }

    excit[0] = 1.f;
    
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();
    
    for (int n = 0;  n < timeSamples; ++n)
    {
        double exc = excit[n];
        
        for (int m = 0 ; m < modesNumber; ++m)
        {
            float c1 = (2 - eigenFreqs[m] * eigenFreqs[m] * k * k) / (dampCoeffs[m] * k + 1.f);
            float c2 = (dampCoeffs[m] * k - 1.f) / (dampCoeffs[m] * k + 1.f);
            float c3 = k * k / (dampCoeffs[m] * k + 1.f);

            xNext[m] = c1 * x[m] + c2 * xPrev[m] + c3 * exc * modesIn[m];

            xPrev[m] = x[m];
            x[m] = xNext[m];
            
            output[n] += xNext[m] * modesOut[m];
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
    audiowrite("plate-saturator-class", output, timeSamples, baseSR);

    return 0;
}
