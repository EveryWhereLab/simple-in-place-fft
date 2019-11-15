#include "simple_in_place_fft.h"

#define SWAP(a,b) { float t; t = a; a = b; b = t; }

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
FftConfig_t *fftInit(int maxWindowSize) {
    double theta;
    unsigned long bsep = 2;
    int fftOrder;
    double tempTrig;
    FftConfig_t *config = (FftConfig_t *)malloc(sizeof(FftConfig_t));
    fftOrder = 0;
    config->fftWindowLength =1;
    while(config->fftWindowLength < (unsigned int)maxWindowSize){
        config->fftWindowLength = config->fftWindowLength << 1;
        fftOrder++;
    }
    // Use double precision initial values to alleviate the accumulated error of the trigonometric recurrences.
    config->initTrigTable[0] = (double *)malloc(sizeof(double) *(fftOrder));
    config->initTrigTable[1] = (double *)malloc(sizeof(double) *(fftOrder));
    for(int i=0;i<fftOrder;i++) {
        theta=(6.28318530717959/bsep);
        tempTrig = -sin(0.5*theta);
        config->initTrigTable[0][i]=-2.0 * tempTrig*tempTrig;
        config->initTrigTable[1][i]=-sin(theta);
        bsep *=2;
    }
    return config;
}

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
int doComplexFFT(FftConfig_t *config, float data[], unsigned long fftWindowLength, int direction){
    unsigned long n, bsep, m, j, bwidth, i;
    double wtemp, wr, wpr, wpi, wi, theta; 
    float tempr, tempi;
    int stage = 0;
    if(fftWindowLength > config->fftWindowLength)
        return -1;
    n = fftWindowLength << 1;
    // In-place bit reverse shuffling of data.
    j = 0;
    for (i = 0; i < n; i+=2) {
        if (i < j) {
            tempr   = data[j];
            data[j] = data[i];
            data[i] = tempr;
            tempi     = data[j+1];
            data[j+1] = data[i+1];
            data[i+1] = tempi;
        }
        //Find j using Gold Radar's algorithm.
        m = n >> 1;
        while (m >= 2 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    bsep = 2;
    while (n > bsep) { // Outer loop: log2(fftWindowLength) stages will be executed.
        bwidth = bsep << 1;
        // Get the initial rotation factor for this stage.
        wpr = config->initTrigTable[0][stage];
        wpi= direction == FFT_FORWARD ? config->initTrigTable[1][stage] : -1.0 * config->initTrigTable[1][stage];
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < bsep; m += 2) {
            for (i = m; i <= n; i += bwidth) {
                // The computation below can be described as a butterfly diagram.
                // Index j points to bottom node and index i points to top node.
                j = i + bsep;
                tempr = wr * data[j-1] - wi * data[j]; // Temporarily Keep the rotated value of bottom node in the butterfly diagram.
                tempi = wr * data[j] + wi * data[j-1]; 
                data[j-1] = data[i-1] - tempr; // Update the value of bottom node.
                data[j] = data[i] - tempi;
                data[i-1] += tempr; // Update the value of top node.
                data[i] += tempi;
            }
            // Obtain the next rotation factor by trigonometric recurrences.
            wtemp = wr;
            wr = wtemp * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        bsep = bwidth;
        stage+=1;
    }
    return stage;
}

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
int doRealFFT(FftConfig_t *config, float data[], unsigned long fftWindowLength, int direction){
    unsigned long i, r_idx, i_idx, mirror_r_idx, mirror_i_idx, np3;
    float c1 = 0.5, c2, h1r, h1i, h2r, h2i;
    double wr, wi, wpr, wpi, wtemp, theta;
    int stage;
    int half_n = (fftWindowLength>>1);
    int quarter_n = (fftWindowLength>>2);
    if (direction == FFT_FORWARD) {
        c2 = -0.5;
        stage = doComplexFFT(config, data,half_n,FFT_FORWARD);
        if(stage==-1)
            return -1;
        wpr = config->initTrigTable[0][stage];
        wpi= config->initTrigTable[1][stage];
    } else {
        c2 = 0.5;
        stage = 0;
        unsigned long bsep;
        bsep = 2;
        while(fftWindowLength > (unsigned int)bsep){
            bsep <<= 1;
            stage++;
        }
        wpr = config->initTrigTable[0][stage];
        wpi= -1.0 *config->initTrigTable[1][stage];
    }
    wr = 1.0 + wpr;
    wi = wpi;
    np3 = fftWindowLength + 1;
    for (i = 1; i< quarter_n; ++i) {
        mirror_i_idx = 1+(mirror_r_idx = np3-(i_idx = 1+(r_idx = i+i)));
        h1r =  c1*(data[r_idx] + data[mirror_r_idx]); 
        h1i =  c1*(data[i_idx] - data[mirror_i_idx]);
        h2r = -c2*(data[i_idx] + data[mirror_i_idx]);
        h2i =  c2*(data[r_idx] - data[mirror_r_idx]);
        data[r_idx] =  h1r + wr * h2r - wi * h2i;
        data[i_idx] =  h1i + wr * h2i + wi * h2r;
        data[mirror_r_idx] =  h1r - wr * h2r + wi * h2i;
        data[mirror_i_idx] = -h1i + wr * h2i + wi * h2r;
        wtemp = wr;
        wr = wtemp * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }
    // Change the sign of the image part at Nyquist bin.
    data[half_n+1] = -1 *data[half_n+1];
    h1r = data[0];
    if (direction == FFT_FORWARD) {
        data[0] = h1r + data[1]; 
        data[1] = h1r - data[1];
    } else {
        data[0] = c1 * (h1r + data[1]);
        data[1] = c1 * (h1r - data[1]);
        return doComplexFFT(config, data, half_n, FFT_BACKWARD)+1; 
    }
    return stage+1;
}

/*
 * Destory function for releasing resources
 * This function should be called when you don't need the allocated FftConfig resource.
 *
 * Parameters
 * ----------
 * FftConfig *config
 *  A ponter points to a FftConfig struct to be free.
 */
void fftDestroy(FftConfig_t *config){
    if(config==NULL)
        return;
    if(config->initTrigTable[0])
        free(config->initTrigTable[0]);
    if(config->initTrigTable[1])
        free(config->initTrigTable[1]);
    free(config);
}

