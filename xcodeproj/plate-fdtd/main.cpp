//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Modal Plate Reverb
// C++ implementation by: Riccardo Russo, University of Bologna
//
// Date: 11-July-2023
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include "plateFDTDData.h"
#include "../../nemus-cpp-audio-tools/src/audio.h"
#include <chrono>

#define PI 3.141592653589793

int main(int argc, const char * argv[])
{
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    // Excitation Signal

#define durSec 2
    
    const unsigned int timeSamples = SR * durSec;
    float excit[timeSamples];

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    double k = 1.0 / SR;
    double k2 = k * k;
    
    float lambda = C * k / h;
    float lambda2 = lambda * lambda;
    
    float sigma0 = 5.f;
    float dampCoeff = 1 + sigma0 * k;
    
    float inPoint[] = {0.52f * Lx, 0.53f * Ly};
    float outPoint[] = {0.52f * Lx, 0.53f * Ly};
    
    std::cout << "Nx: " << Nx << '\n';
    std::cout << "Ny: " << Ny << '\n';
    
    //// Simulation
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    /// Initializing vectors

    float output[timeSamples];

    float u[Nx][Ny];
    float uPrev[Nx][Ny];
    float uNext[Nx][Ny];
    
    float Jcoeff[Nx][Ny];
    
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            u[i][j] = 0.f;
            uPrev[i][j] = 0.f;
            uNext[i][j] = 0.f;
            Jcoeff[i][j] = 0.f;
        }
    }
    
    //Interpolation operators
    int lo = (int) inPoint[0] / h;
    int mo = (int) inPoint[1] / h;
    float alphax = inPoint[0] / h - lo;
    float alphay = inPoint[1] / h - mo;
    Jcoeff[lo][mo] = (1 - alphax) * (1 - alphay) / (h * h);
    Jcoeff[lo][mo + 1] = (1-alphax)*alphay / (h * h);
    Jcoeff[lo + 1][mo] = alphax*(1-alphay) / (h * h);
    Jcoeff[lo + 1][mo + 1] = alphax*alphay / (h * h);
    
    int loO = (int) outPoint[0] / h;
    int moO = (int) outPoint[1] / h;
    float alphaxO = outPoint[0] / h - loO;
    float alphayO = outPoint[1] / h - moO;

    //Exctitation signal
    for (int i = 0; i < timeSamples; ++i)
    {
        excit[i] = 0.f;
        output[i] = 0.f;
    }
    excit[0] = 1.f;

    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();
    
    for (int n = 0;  n < timeSamples; ++n)
    {
        double exc = excit[n];
        
        for (int i = 1; i < Nx - 1; ++i)
            for (int j = 1; j < Ny - 1; ++j)
            {
                int I = 0;
                if ((i ==5) && (j ==5)) I = 1;
                uNext[i][j] = 2 * (1 - 2 * lambda2) * u[i][j] / dampCoeff + (sigma0 * k - 1) * uPrev[i][j] / dampCoeff + lambda2 * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1]) / dampCoeff + k2 * Jcoeff[i][j] * exc / dampCoeff;
            }
        
        output[n] = (1 - alphaxO) * (1 - alphayO) * uNext[loO][moO] + (1 - alphaxO) * alphayO*uNext[loO][moO + 1] + alphaxO * (1 - alphayO)*uNext[loO + 1][moO] + alphaxO * alphayO * uNext[loO + 1][moO + 1];
        
        for (int i = 1; i < Nx - 1; ++i)
        {
            for (int j = 1; j < Ny - 1; ++j)
            {
                uPrev[i][j] = u[i][j];
                u[i][j] = uNext[i][j];
            }
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
    audiowrite("2dwave-fdt", output, timeSamples, SR);

    return 0;
}
