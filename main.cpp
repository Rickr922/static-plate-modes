//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Modal Plate Reverb
// C++ implementation by: Riccardo Russo, University of Bologna
//
// Date: 11-July-2023
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include <cmath>
#include <complex>
#include "nemus-cpp-audio-tools/src/audio.h"
#include <chrono>

///
enum class NonLinearityFunction
{
    Cubic, Tanh, Sinh, Exp, Linear
};

///
enum class PlateMaterial
{
    Steel1, Steel2, Gold, Silver, Copper, Aluminium,
};

///
struct GridPoint
{
    double x;
    double y;
};

struct ModalState
{
    double q = 0;
    double p = 0;
};

///
struct MaterialCoefficients
{
    ///
    double E;
    ///
    double ni;
    ///
    double rho;
    ///
    double R1;
    ///
    double C1;
};

MaterialCoefficients getMaterialData(PlateMaterial);
double computeOmega(double m1, double m2, double T, double rho, double D, double Lx, double Ly, double Lz);
double computeMode(double xp,double yp, double m1, double m2, double Lx, double Ly);
double computeTHDamp(double R1, double C1, double* freqs, double thick);
double computeRadDamp(double CF, double* freqs, double rhoA, double  cA, double  rhoPlate, double thick, double widthX, double widthY);
double computePorousDamp(double freq, double rhoA, double cA, bool modelTypeCA, double thick, double rhoPlate, double Dplate, double modes, double Dist, double widthPor);
int getUpperModeCount(double dim, double upperAngularFreqLimit, double T, double D, double rho, double Lz);

int main(int argc, const char * argv[])
{
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    // Excitation Signal
    constexpr double amp = 50000.0;
    constexpr double osFac = 2.0;
    constexpr double durSec = 1.0;
    double SR = 48000.0;
    unsigned int timeSamples = std::floor(SR * osFac * durSec);
    double* excit = new double[timeSamples];
    std::fill(excit, excit + timeSamples, 0);
    excit[0] = amp;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    SR *= osFac;
    double k = 1.0 / SR;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    // Custom Parameters
    NonLinearityFunction nLinType = NonLinearityFunction::Cubic;
    
    // sigma0 = 0;
    double cCoeff = 0.1;
    double a = 30;
    
    double Lx = 2.0;    //[m] Hor length
    double Ly = 1.0;    //[m] Ver lentgh
    double Lz = 5.0e-4; //[m] Thickness
    
    MaterialCoefficients materialData = getMaterialData(PlateMaterial::Steel2);
    
    double& E   = materialData.E;
    double& ni  = materialData.ni;
    double& rho = materialData.rho;
    double T = 600.0;                      //[N] Tension
    double D = E * (Lz*Lz*Lz) / (12.0 * (1.0 - (ni*ni)));     //Flexural Rigidity
    
    constexpr double omegaUpperLimit = 10000.0 * 2.0 * nemus::pi; // 2/k;
    constexpr double omegaLowerLimit = (20.0 * 2.0 * nemus::pi);
    double modeLimit = 0.6;
    
    GridPoint inPoint   = {0.52 * Lx, 0.53 * Ly};
    GridPoint outPoint = {0.47 * Lx, 0.62 * Ly};
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    // Computing eigenfrequencies and eigenvectors
    
    // Computing max number of modes
    unsigned int Nx = getUpperModeCount(Lx, omegaUpperLimit, T, D, rho, Lz);
    unsigned int Ny = getUpperModeCount(Ly, omegaUpperLimit, T, D, rho, Lz);
    unsigned int Mx = 1;
    unsigned int My = 1;
    double omegaCurr = 0;
    
    std::cout << "Nx: " << (Nx-Mx) << " Ny: " << (Ny-My) << '\n';
    
    // Computing eigenfrequencies

    int M = 0;
    
    for (int i = Mx; i <= Nx; i++)
    {
        for (int j = My; j <= Ny; j++)
        {
            omegaCurr = computeOmega(i, j, T, rho, D, Lx, Ly, Lz);
            if (    omegaCurr < omegaUpperLimit
                and omegaCurr > omegaLowerLimit)
            {
                M++;
            }
        }
    }
    
    std::cout << "Number of Modes: " << M << '\n';
    
    // Computing Modes for in and out points
    constexpr double rhoAir = 1.225;  // [kg/m^3]
    constexpr double cAir   = 343.0;    // [m/s]
    
    const double dampCoef1 = (materialData.R1 * materialData.C1);
    const double Lz2 = (Lz*Lz);
    const double dampCoef3 = (materialData.C1*materialData.C1) / (Lz2);
    const double kappa = std::sqrt(D / (rho * Lz));
    const double criticalFreq = (cAir*cAir) / (2.0 * nemus::pi * kappa);
    const double icf = 1.0 / criticalFreq;
    constexpr double dist = 0.005;   //[m] distance between plate and medium
    constexpr double H = 0.1;        //[m] Medium thickness
    
    alignas(32) double* omega  = new double[M];
    alignas(32) double* sigma =  new double [M];
    alignas(32) double* g_mi = new double[M];
    alignas(32) double* g_mo = new double[M];
    
    int m = 0;
    for (int i = Mx; i <= Nx; i++)
    {
        for (int j = My; j <= Ny; j++)
        {
            double omegaCurr = computeOmega(i, j, T, rho, D, Lx, Ly, Lz);
            
            if (    omegaCurr < omegaUpperLimit
                and omegaCurr > omegaLowerLimit)
            {
                
                double modeIn  = computeMode(inPoint.x,
                                             inPoint.y,
                                             i,
                                             j,
                                             Lx,
                                             Ly);
                double modeOut = computeMode(outPoint.x,
                                             outPoint.y,
                                             i,
                                             j,
                                             Lx,
                                             Ly);
                
                //Thermoelastic Damping
                double freq = omegaCurr;
                double freq2 = freq*freq;
                double sigmaTH = (freq2 * dampCoef1) / (2.0 * ((freq2 * Lz2) + dampCoef3));
                
                //Radiation Damping
                double psi = std::sqrt((freq * 0.5 / nemus::pi) * icf);
                double g =
                ((1.0 - (psi*psi)) * std::log((1.0 + psi) / (1.0 - psi)) + 2.0 * psi)
                / std::pow((1.0 - psi), 1.5); //std::sqrt((1.0 - psi)*(1.0 - psi)*(1.0 - psi));
                
                double sigmaRad = (1.0 / (4.0 * nemus::pi*nemus::pi))
                * (rhoAir * cAir / (rho * Lz))
                * (2.0 * (Lx + Ly) / (Lx * Ly))
                * (cAir * icf)
                * g;
                
                //Radiation into a Porous Medium
                double sigmaPorous = computePorousDamp(omegaCurr, rhoAir, cAir, true, Lz, rho, D, 1.0, dist, H);
                
                //Full damping coeffs
                //divided by 2 because here 2sigma0 is considered (Bilbao convention)
                double sigma0s = (sigmaRad + sigmaTH + sigmaPorous) * 0.5;
                
                omega[m] = 0.5 * k * omegaCurr;
                sigma[m] = 0.5 * k * sigma0s;
                g_mi [m] = modeIn * k / (rho * Lz);
                g_mo [m] = modeOut / omegaCurr;
                
                m++;
            }
        }
    }
    
    
    //// Simulation
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    /// Initializing vectors
    ModalState* x = new ModalState[M];
    //    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    /// Simulation
    double* output = new double[timeSamples];
    std::fill(output, output+timeSamples, 0.0);
    alignas(32) double* p = new double [M];
    std::fill(p, p+M, 0);
    alignas(32) double* q = new double [M];
    std::fill(q, q+M, 0);
    
    double maxSample = 0.0;
    
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();
    
    for (int n = 0;  n < timeSamples; ++n)
    {
        double exc = excit[n];
        
        for (int m = 0 ; m < M; ++m)
        {
            double& omega_m =  omega[m];
            double& sigma_m =  sigma[m];
            double& q_m =  q[m];
            double& p_m =  p[m];
            double& gi_m =  g_mi[m];
            
            double eta = (a * p_m);
            double d = 1.0 + eta*eta;
            double lambda = 1.0 + 3.0*eta*eta;
            
            double detCoeff = 1.0 / (1.0 + (a * sigma_m * lambda) + (omega_m * omega_m));
            
            double qn = detCoeff * (
                                    + omega_m * (
                                                 + (gi_m   * exc)
                                                 - (omega_m * q_m)
                                                 + (2.0f * p_m * (1.0f + (sigma_m * lambda) - (sigma_m * d)))
                                                 )
                                    + (q_m * (1.0f + sigma_m * lambda))
                                    );
            
            double pn = detCoeff * (
                                    + (gi_m * exc)
                                    - (2.0f * omega_m * q_m)
                                    + p_m * (
                                             + (1.0f)
                                             + (sigma_m * lambda)
                                             - (omega_m * omega_m)
                                             - (2.0f * sigma_m * d)
                                             )
                                    );
            
            q_m = qn;
            p_m = pn;
            
            output[n] += q_m * g_mo[m];
        }
        
        if (std::abs(output[n]) > maxSample) maxSample =  std::abs(output[n]);
    }
    
    std::chrono::steady_clock::time_point toc  = std::chrono::steady_clock::now();
    double computeTime = std::chrono::duration<double >(std::chrono::duration_cast<std::chrono::nanoseconds>(toc - tic)).count();
    std::cout << "Compute Time: " << computeTime << '\n';
    
    maxSample = 1.0 / maxSample;
    
    
    for (int i = 0; i < timeSamples; i++) {
        output[i] *= maxSample;
    }
    
    float* downsampled = resample(output, SR / osFac, SR, timeSamples);
    diff(downsampled, timeSamples);
    audiowrite("plate-saturator", downsampled, timeSamples, SR / osFac);
}

//-------------------------------------------------------------------------
// Functions

double computeOmega(double m1, double m2, double T, double rho, double D, double Lx, double Ly, double Lz)
{
    double C1 = T / (rho * Lz);
    double C2 = D / (rho * Lz);
    double C3 =
    (((m1*m1) * (nemus::pi*nemus::pi)) / (Lx*Lx))
    + (((m2*m2) * nemus::pi*nemus::pi) / (Ly*Ly));
    
    
    return std::sqrt((C1 * C3) + (C2 * (C3*C3)));
}

double computeMode(double xp,double yp, double m1, double m2, double Lx, double Ly)
{
    return std::sqrt(4.0 / (Lx * Ly)) * std::sin(m1 * nemus::pi * xp / Lx) * std::sin(m2 * nemus::pi * yp / Ly);
}

double computeTHDamp(double R1, double C1, double freq, double thick)
{
    double coeff1 =  (freq * freq) * (R1*C1);
    double coeff2 = ((freq * freq) * (thick*thick));
    return coeff1 / (2.0 * (coeff2 + ((C1*C1) / thick * thick)));
    
}
//
double computeRadDamp(double CF, double* freqs, double rhoA, double  cA, double  rhoPlate, double thick, double widthX, double widthY)
{
    //    double psi = std::sqrt((freqs / 2.0 / nemus::pi) / CF);
    //    double g = ((1.0 - (psi*psi)) .* std::log((1.0 + psi) ./ (1.0 - psi)) + 2.0 * psi)
    //    ./ std::pow((1.0 - psi), 1.5); //std::sqrt((1.0 - psi)*(1.0 - psi)*(1.0 - psi));
    //    return (1.0 / 4.0 / nemus::pi2)
    //    * (rhoA * cA / rhoPlate / thick)
    //    * (2.0 * (widthX + widthY) / (widthX * widthY))
    //    * (cA / CF)
    //    * g;
    
    return 0.0;
}

double computePorousDamp(double freq, double rhoA, double cA, bool modelTypeCA, double thick, double rhoPlate, double Dplate, double mode, double Dist, double widthPor)
{
    using namespace std::complex_literals;
    std::complex<double> waveNumAir = freq / cA;
    std::complex<double> rhoMed = 0;
    std::complex<double> cMed = 0;
    if (modelTypeCA)
    {
        // Craik-Allard
        constexpr double flowRes = 20000.0; // [Rayl/m]
        double csi = rhoA * (freq / (2.0 * nemus::pi)) / flowRes;
        
        constexpr std::complex<double> c1{     1.2000, 0.0};
        constexpr std::complex<double> c2{    -0.0364, 0.0};
        constexpr std::complex<double> c3{     0.0000, 0.1144};
        constexpr std::complex<double> c4{101320.0000, 0.0};
        constexpr std::complex<double> c5{     0.0000, 29.64};
        constexpr std::complex<double> c6{     0.0000, 21.17};
        constexpr std::complex<double> c7{     2.8200, 0.0};
        constexpr std::complex<double> c8{     0.0000, 24.9};
        std::complex<double> c9 = std::sqrt(c7 / (csi*csi) + c8 / csi);
        
        rhoMed = c1 + std::sqrt((c2 / (csi*csi)) - (c3 / csi));
        
        cMed = std::sqrt((c4 / rhoMed) * (c5 + c9) / (c6 + c9));
    }
    std::complex<double> waveNumPlate = std::sqrt(freq / std::sqrt(Dplate / (rhoPlate * thick)));
    std::complex<double> waveNumPorous = freq / cMed;
    
    std::complex<double>  waveNumGap = std::sqrt((waveNumAir*waveNumAir) - (waveNumPlate*waveNumPlate));
    std::complex<double>  waveNumMed = std::sqrt((waveNumPorous*waveNumPorous)
                                                 - (waveNumPlate*waveNumPlate));
    
    std::complex<double> constant = (waveNumMed * rhoA) / (waveNumGap * rhoMed);
    std::complex<double> constant2 = (constant + 1.0) / (constant - 1.0);
    std::complex<double> constant3 = (constant - 1.0) / (constant + 1.0);
    std::complex<double> BOverA = (-std::exp(-2.0 * 1i * waveNumGap * Dist)
                                   * (1.0 - std::exp(-2.0 * 1i * waveNumMed * widthPor)))
    / (constant2 - std::exp(-2.0 * 1i * waveNumMed * widthPor) * constant3);
    std::complex<double> F = (1.0 / waveNumGap) * ((1.0 + BOverA) / (1.0 - BOverA));
    std::complex<double> radEfficiency = waveNumAir * std::real(F);
    
    std::complex<double> t60 = 13.82 * (rhoPlate * thick) / (rhoA * cA * radEfficiency);
    
    return 6.0 * std::log(10.0) / std::real(t60);
}

MaterialCoefficients getMaterialData(PlateMaterial material)
{
    switch (material)
    {
        case PlateMaterial::Steel1:
            return {2e11, 0.3, 8.05e3, 9.416e-3, 0.14965e-3};
            break;
        case PlateMaterial::Steel2:
            return {2e11, 0.3, 7.872e3, 9.664e-3, 0.1855e-3};
            break;
        case PlateMaterial::Gold:
            return {1.6e11, 0.42, 7.872e3, 4.727e-3, 1.270e-3};
            break;
        case PlateMaterial::Silver:
            return {1e11, 0.38, 10.49e3, 8.403e-3, 1.679e-3};
            break;
        case PlateMaterial::Copper:
            return {1.4e11,0.35,8.96e3,5.691e-3,1.148e-3};
            break;
        case PlateMaterial::Aluminium:
            return {8e10, 0.34, 2.7e3, 9.975e-3, 0.976e-3};
            break;
        default:
            return {2e11, 0.3, 8.05e3, 9.416e-3, 0.14965e-3};
    }
}

int getUpperModeCount(double dim, double upperAngularFreqLimit, double T, double D, double rho, double Lz)
{
    double b = (T / (rho * Lz))*(nemus::pi*nemus::pi) / (dim*dim);
    double c = (D / (rho * Lz))*((nemus::pi*nemus::pi)/(dim*dim))*((nemus::pi*nemus::pi) / (dim*dim));
    double a = upperAngularFreqLimit*upperAngularFreqLimit;
    return std::floor(std::sqrt((std::sqrt((4.0 * a * c) + b*b) - b) / (2*c)));
}