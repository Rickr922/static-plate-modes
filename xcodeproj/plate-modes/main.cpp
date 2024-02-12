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
#define exactOsc true
    
#define osFac 1
#define durSec 2
#define baseSR 44100
    const unsigned int timeSamples = baseSR * osFac * durSec;
    
#if exactOsc
    float excit[timeSamples + 2];
    excit[0] = 0.f;
    excit[1] = 0.f;
#else
    float excit[timeSamples];
#endif

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    double SR = osFac * baseSR;
    double k = 1.0 / SR;
    
    int modesNumber = modesNumberFull;

#if !exactOsc
    // STABILITY CONDITION FOR NON EXACT OSCILLATOR!!!
    for(int i = 0; i < modesNumberFull; ++i)
    {
        if(eigenFreqs[i] > 2/k)
        {
            modesNumber = i - 1;
            break;
        }
            
    }
#endif
    
    std::cout << "Number of Modes: " << modesNumber << '\n';
    
    //// Simulation
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    /// Initializing vectors

    float output[timeSamples];

    float x[modesNumber];
    float xPrev[modesNumber];
    float xNext[modesNumber];
    
    float c1[modesNumber];
    float c2[modesNumber];
    float c3[modesNumber];
    
#if exactOsc
    float c1exc[modesNumber];
    float c2exc[modesNumber];
    float c3exc[modesNumber];
#endif

    for (int i = 0; i < timeSamples; ++i)
    {
#if exactOsc
        excit[i + 2] = 0.f;
#else
        excit[i] = 0.f;
#endif
        output[i] = 0.f;

        if (i < modesNumber)
        {
            x[i] = 0.f;
            xPrev[i] = 0.f;
            xNext[i] = 0.f;
            
#if exactOsc
            c1[i] = 2.f * exp(-dampCoeffs[i] * k)*cos(sqrt((eigenFreqs[i] * eigenFreqs[i]) - (dampCoeffs[i] * dampCoeffs[i])) * k);
            c2[i] = -exp(-2.f * dampCoeffs[i] * k);
            c3[i] = k * k * modesIn[i];
            
            c1exc[i] = (1.f + k * dampCoeffs[i]) / 12.f;
            c2exc[i] = - (5.f / 6.f) - dampCoeffs[i] * k + (2.f / 3.f) * k * k * dampCoeffs[i] * dampCoeffs[i] - k * k * eigenFreqs[i] * eigenFreqs[i] / 12.f;
            c3exc[i] = (1.f - k * dampCoeffs[i]) / 12.f;
#else
            c1[i] = (2.f - eigenFreqs[i] * eigenFreqs[i] * k * k) / (dampCoeffs[i] * k + 1.f);
            c2[i] = (dampCoeffs[i] * k - 1.f) / (dampCoeffs[i] * k + 1.f);
            c3[i] = k * k * modesIn[i] / (dampCoeffs[i] * k + 1.f);
#endif
        }
    }

    excit[7] = 1.f;
    
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();
    
    for (int n = 0;  n < timeSamples; ++n)
    {
#if !exactOsc
        double exc = excit[n];
#endif
        for (int m = 0 ; m < modesNumber; ++m)
        {
#if exactOsc
            double exc = excit[n + 2]*c1exc[m] + excit[n + 1]*c2exc[m] + excit[n]*c3exc[m];
#endif
            xNext[m] = c1[m] * x[m] + c2[m] * xPrev[m] + c3[m] * exc;

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
    audiowrite("plate", output, timeSamples, baseSR);

    return 0;
}
