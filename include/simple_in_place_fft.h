/*
 * Copyright (C) 2019 EveryWhereLab. All rights reserved.
 *
 * Simple FFT/IFFT functions and tests are provided.
 * This implementation is based on the algorithms described in chapter 12 of the numerical recipes in C book.
 * There are slight differences between the original implementation and this one in the following aspects:
 * 1. The array index starts from zero in this implementation.
 * 2. A init. function is added to create the initial trigonometric tables for later use.
 *
 * Licensed under the Apache License, Version 2.0 (the License); you may
 * not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __SIMPLE_INPLACE_FFT_H__
#define __SIMPLE_INPLACE_FFT_H__
#include <stdlib.h>
#include <math.h>

typedef enum { 
    FFT_FORWARD,
    FFT_BACKWARD
} FftDirection_t;

typedef struct {
    int fftWindowLength;
    double *initTrigTable[2];
} FftConfig_t;

/*
 * Calling this Init. function for preparing trigonometric tables
 * This function must be called prior to other FFT calculations.
 *
 * Parameters
 * ----------
 *  int maxWindowSize
 *    The maixmum fft length will be used in sebsequent FFT calculations.
 *
 * Retrun
 *    A ponter points to a FftConfig struct for sebsequent FFT calculations.
 */
FftConfig_t *fftInit(int maxWindowSize);

/*
 * Fast Fourier transform function for complex samples
 * This function is a radix-2 in-place implementation.It replaces complex input samples by its FFT result if direction is
 * FF_FORWARD; or replaces complex input samples by fftWindowLength times its IFFT result if direction is FFT_BACKWARD.
 *
 * Parameters
 * ----------
 *  FftConfig *config
 *    A ponter points to a FftConfig struct returned from the fftInit() function.
 *  float data[]
 *    The in/out array containing the complex samples with
 *    real/imaginary parts interleaved [Re(0), Im(0), ..., Re(fftWindowLength-1), Im(fftWindowLength-1)]
 *    The fftWindowLength should be a power of 2 and must not be larger than the maixmum fft length in your configuration 
 *  unsigned long fftWindowLength
 *    The FFT length of the input data.
 *  int direction
 *    This value should be FFT_FORWARD or FFT_BACKWARD.
 *
 * Return
 *    The returned value will be log2(fftWindowLength) in a successful case. Or it will return -1 to indicate an improper
 *    configuration.
 */
int doComplexFFT(FftConfig_t *config, float data[], unsigned long fftWindowLength, int direction);

/*
 * Fast Fourier transform function for real samples
 * This function is a radix-2 in-place implementation. It replaces real input samples by its complex FFT result if direction
 * is FFT_FORWARD. If input is a FFT result of real samples and set direction to FFT_BACKWARD, the complex input samples will
 * be replaced by a recovered result of IFFT(the result still needs to be multiplied by fftWindowLength/2).
 *
 * Parameters
 * ----------
 *  FftConfig *config
 *    A ponter points to a FftConfig struct returned from the fftInit() function.
 *  float data[]
 *    The in/out array containing the samples.
 *    If direction is FFT_FORWARD, input data will be interpreted as real samples [Re(0),..,Re(fftWindowLength-1)]. The real input
 *    samples will be replaced by its coplex FFT result[Re(0),Im(0),..,Re(fftWindowLength/2-1),Im(fftWindowLength/2-1)]. The real
 *    part of Nyquist frequency bin Re(fftWindowLength/2) is stored at Im(0). If direction is FFT_BACKWARD, input data will be
 *    interpreted as complex samples [Re(0),Im(0),..,Re(fftWindowLength/2-1),Im(fftWindowLength/2-1)]. The complex input samples
 *    will be replaced by its IFFT result.
 *  unsigned long fftWindowLength
 *    The FFT length of the input data.
 *  int direction
 *    This value should be FFT_FORWARD or FFT_BACKWARD.
 *
 * Return
 *    The returned value will be log2(fftWindowLength) in a successful case. Or it will return -1 to indicate an improper
 *    configuration.
 */
int doRealFFT(FftConfig_t *config, float data[], unsigned long fftWindowLength, int direction);

/*
 * Destory function for releasing resources
 * This function should be called when you don't need the allocated FftConfig resource.
 *
 * Parameters
 * ----------
 * FftConfig *config
 *  A ponter points to a FftConfig struct to be free.
 */
void fftDestroy(FftConfig_t *config);

#endif
