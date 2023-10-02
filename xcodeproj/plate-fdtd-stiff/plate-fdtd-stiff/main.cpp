//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Modal Plate Reverb
// C++ implementation by: Riccardo Russo, University of Bologna
//
// Date: 11-July-2023
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include "plateFDTDStiffData.h"
#include "../../../nemus-cpp-audio-tools/src/audio.h"
#include <chrono>

#define PI 3.141592653589793

int main(int argc, const char * argv[])
{
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    // Excitation Signal

#define durSec 2
    
    const unsigned int timeSamples = SR * durSec;
    float excit[timeSamples] = {0};

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    static double k = 1.0 / SR;
    static double k2 = k * k;

    //static float lambda = C * k / h;
    //static float lambda2 = lambda * lambda;
    static float sigma0 = 5.f;
    static float dampCoeff = 1 + sigma0 * k;

    static float inPoint[] = {
        0.52f * Lx,
        0.53f * Ly
    };
    static float outPoint[] = {
        0.52f * Lx,
        0.53f * Ly
    };

    
    std::cout << "Nx: " << Nx << '\n';
    std::cout << "Ny: " << Ny << '\n';
    
    //// Simulation
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    /// Initializing vectors

    float output[timeSamples];

    static float u[totN]      = {0};
    static float uPrev[totN]  = {0};
    static float uNext[totN]  = {0};
    static float Jcoeff[Nx-1][Ny-1] = {0};
    static float Jvec[totN] = {0};
    
    //Interpolation operators
    int lo = (int) inPoint[0] / h;
    int mo = (int) inPoint[1] / h;
    float alphax = inPoint[0] / h - lo;
    float alphay = inPoint[1] / h - mo;
    Jcoeff[lo][mo] = (1 - alphax) * (1 - alphay) / (h * h);
    Jcoeff[lo][mo + 1] = (1-alphax)*alphay / (h * h);
    Jcoeff[lo + 1][mo] = alphax*(1-alphay) / (h * h);
    Jcoeff[lo + 1][mo + 1] = alphax*alphay / (h * h);
    
    for (int i = 0; i < Nx-2; ++i) {
         for (int j = 0; j < Ny-2; ++j) {
             Jvec[j + (Ny-1) * i] = Jcoeff[i][j];
         }
    }
    
    //Exctitation signal
    excit[0] = 1.f;

    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();
    
    for (int n = 0;  n < timeSamples; ++n)
    {
        double exc = excit[n];
        // Compute uNext
        for (int i = 0; i < totN; ++i)
        {
            float sum1 = 0;
            float sum2 = 0;

            for (int j = 0; j < totN; ++j)
            {
                sum1 = sum1 + A[i][j] * u[j];
                sum2 = sum2 + S[i][j] * uPrev[j];
            }
            uNext[i] =   sum1 + sum2 + k2 * Jvec[i] * exc/dampCoeff;
        }
        // Compute outputs
        output[n] = uNext[(int)totN/2];
        // Swap buffers
        for (int i = 0; i < totN; ++i)
        {
             uPrev[i] = u[i];
             u[i] = uNext[i];
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
    audiowrite("platestiff-fdtd", output, timeSamples, SR);

    return 0;
}
