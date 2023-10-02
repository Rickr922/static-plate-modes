/************************************************************************
 ************************************************************************
    Syfala compilation flow
    Copyright (C) 2022 INSA-LYON, INRIA, GRAME-CNCM
---------------------------------------------------------------------
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 ************************************************************************
 ************************************************************************/

#include <iostream>
#include <algorithm>
#include <cmath>
#include <inttypes.h>
#include <string.h>
#include <syfala/utilities.hpp>
#include "plate-fdtd.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#define PI 3.141592653589793

#define FAUST_INT_CONTROLS 0
#define FAUST_REAL_CONTROLS 0
#define FAUST_INT_ZONE 0
#define FAUST_FLOAT_ZONE 0
#define FAUST_INPUTS 0
#define FAUST_OUTPUTS 2
#define FAUST_ACTIVES 0
#define FAUST_PASSIVES 0

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

static float excit[SYFALA_BLOCK_NSAMPLES] = {0};

static float u[totN]      = {0};
static float uPrev[totN]  = {0};
static float uNext[totN]  = {0};
static float Jcoeff[Nx-1][Ny-1] = {0};
static float Jvec[totN] = {0};

static int lo = (int) inPoint[0] / h;
static int mo = (int) inPoint[1] / h;
static float alphax = inPoint[0] / h - lo;
static float alphay = inPoint[1] / h - mo;

static void initialize_dsp() {
    Jcoeff[lo][mo] = (1 - alphax) * (1 - alphay) / (h * h);
    Jcoeff[lo][mo + 1] = (1-alphax)*alphay / (h * h);
    Jcoeff[lo + 1][mo] = alphax*(1-alphay) / (h * h);
    Jcoeff[lo + 1][mo + 1] = alphax*alphay / (h * h);
    for (int i = 0; i < Nx-2; ++i) {
         for (int j = 0; j < Ny-2; ++j) {
             Jvec[j + (Ny-1) * i] = Jcoeff[i][j];
         }
    }
    excit[0] = 1.f;
}

static int loO = (int) outPoint[0] / h;
static int moO = (int) outPoint[1] / h;
static float alphaxO = outPoint[0] / h - loO;
static float alphayO = outPoint[1] / h - moO;

static bool initialize = true;

// ----------------------------------------------------------------------------
// HLS TOP-LEVEL FUNCTION
// ----------------------------------------------------------------------------

void syfala (
    sy_ap_int audio_out_0[SYFALA_BLOCK_NSAMPLES],
    sy_ap_int audio_out_1[SYFALA_BLOCK_NSAMPLES],
        float arm_control_f[2],
          int arm_control_i[2],
        float arm_control_p[2],
         int* control_block,
          int arm_ok,
        bool* i2s_rst,
       float* mem_zone_f,
         int* mem_zone_i,
         bool bypass,
         bool mute,
         bool debug
){
#pragma HLS INTERFACE ap_fifo port=audio_out_0
#pragma HLS INTERFACE ap_fifo port=audio_out_1
#pragma HLS INTERFACE s_axilite port=arm_control_f
#pragma HLS INTERFACE s_axilite port=arm_control_i
#pragma HLS INTERFACE s_axilite port=arm_control_p
#pragma HLS INTERFACE s_axilite port=arm_ok
#pragma HLS INTERFACE s_axilite port=control_block
#pragma HLS INTERFACE m_axi port=mem_zone_f latency=30 bundle=ram
#pragma HLS INTERFACE m_axi port=mem_zone_i latency=30 bundle=ram

    // Active high reset, this HAVE TO BE DONE FIRST (crash with *some* dsp if not)
    *i2s_rst = !arm_ok;

    /* RAM must be enabled by ARM before any computation */
    if (arm_ok) {
        if (initialize) {
            /* first iteration: constant initialization */
            initialize_dsp();
            initialize = false;
        } else {
            /* All other iterations:
             * - update controllers values from IP ports
             * - compute one block of N samples
             * - write back passive controller values */
            if (bypass || mute) {
                for (int n = 0; n < SYFALA_BLOCK_NSAMPLES; ++n) {
                    audio_out_0[n] = 0;
                    audio_out_1[n] = 0;
                }
            } else {
                static float output[SYFALA_BLOCK_NSAMPLES] = {0};

                for (int n = 0;  n < SYFALA_BLOCK_NSAMPLES; ++n) {
                    double exc = excit[n];
                    // Compute uNext
                    for (int i = 0; i < totN-1; ++i) {
                        for (int j = 0; j < totN-1; ++j) {
                             uNext[i] =   A[i][j] * u[j]
                                        + S[i][j] * uPrev[j]
                                        + k2 * Jvec[i] * exc/dampCoeff
                                        + uNext[i];
                        }
                    }
                    // Compute outputs
                    output[n] = uNext[(int)totN/2];
                    // Swap buffers
                    for (int i = 0; i < totN-1; ++i) {
                         uPrev[i] = u[i];
                         u[i] = uNext[i];
                    }
                    float out = output[n] * SCALE_FACTOR;
                    audio_out_0[n] = sy_ap_int(out);
                    audio_out_1[n] = sy_ap_int(out);
                }
            }
        }
    }
}
