#include "simple_in_place_fft.h"
#include <time.h>
#include <stdio.h>

#define TEST_FFT_LENGTH 256

int main( int argc, char *argv[] ) 
{
    FftConfig_t *config = (FftConfig_t *) fftInit(TEST_FFT_LENGTH);
    float dataArrayForTest[TEST_FFT_LENGTH];
    float dataArrayPattern[TEST_FFT_LENGTH];
    float a = 10000.0;
    int stage;
    int recoveredErrorCount = 0;
    srand( time(NULL) );
    printf("Generating a test sequence...\r\n");
    for(int i=0;i<TEST_FFT_LENGTH;i++) {
        dataArrayForTest[i] = (float)rand()/(float)(RAND_MAX/a);
        dataArrayPattern[i] = dataArrayForTest[i];
    }
    printf("Applying forward FFT on the test sequence...\r\n");
    stage = doRealFFT(config, dataArrayForTest,TEST_FFT_LENGTH,FFT_FORWARD);
    if(stage==-1) {
        printf("Error occurred in invoking applyRealFFT!\r\n");
        exit(-1);
    }
    else
        printf("Successfully invoking applyRealFFT for %d-stage FFT calculation.\r\n",stage);

    printf("Applying backward FFT on the transform-domain...\r\n");
    stage = doRealFFT(config, dataArrayForTest,TEST_FFT_LENGTH,FFT_BACKWARD);
    if(stage==-1) {
        printf("Error occurred in invoking applyRealFFT!\r\n");
        exit(-1);
    }
    else
        printf("Successfully invoking applyRealFFT for %d-stage IFFT calculation.\r\n",stage);



    printf("Normalizing the signal...\r\n");
    for(int i=0;i<TEST_FFT_LENGTH;i++) {
        dataArrayForTest[i]/=(TEST_FFT_LENGTH>>1);
    }
    printf("Checking if the recovered sequence is the same as the original one...\r\n");
    for(int i=0;i<TEST_FFT_LENGTH;i++) {
        if(abs(dataArrayForTest[i] - dataArrayPattern[i]) > 1e-5) {
            recoveredErrorCount++;
        }
    }
    if(recoveredErrorCount>0)
        printf("Error occurred! %d points are different.\r\n", recoveredErrorCount);
    else
        printf("The recovered sequence is almost the same as the original one.\r\n");
    fftDestroy(config);
    return 0;
}

