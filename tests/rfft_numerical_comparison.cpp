#include <time.h>
#include <Python.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

extern "C" {
#include "simple_in_place_fft.h"
}

using namespace std; 

#define CHECK_NULL_ASSERT(p) \
    if (NULL == (p)) {\
        PyErr_Print();\
        assert(0);\
    }

void init_numpy(){
    import_array();
}

#define TEST_FFT_LENGTH 16

int main( int argc, char *argv[] ) 
{
    FftConfig_t *config = (FftConfig_t *) fftInit(TEST_FFT_LENGTH);
    float realFftInOutArray[TEST_FFT_LENGTH];
    float fftpackInArray[TEST_FFT_LENGTH];
    float fftpackRealOutArray[TEST_FFT_LENGTH];
    float fftpackImagOutArray[TEST_FFT_LENGTH];
    srand( time(NULL) );
    float a = 1000.0;
    int stage;
    cout<<"Generating a test sequence..."<<endl;
    cout<< setiosflags(ios::fixed);
    cout<<"{";
    for(int i=0;i<TEST_FFT_LENGTH;i++) {
        realFftInOutArray[i] = (float)rand()/(float)(RAND_MAX/a);
        fftpackInArray[i] = (double)realFftInOutArray[i];
        cout<<realFftInOutArray[i];
        cout<< (i==TEST_FFT_LENGTH-1? "" : ",");
        if((i+1)%8==0 && (i+1)<TEST_FFT_LENGTH)
            cout<<endl;
    }
    cout<<"}"<<endl;
    cout<<endl<<"Applying forward FFT on the test sequence using the applyRealFFT function..."<<endl;
    stage = doRealFFT(config, realFftInOutArray,TEST_FFT_LENGTH,FFT_FORWARD);
    if(stage==-1) {
        cout<<"Error occurred in invoking applyComplexFFT!"<<endl;
        exit(EXIT_FAILURE);
    }
    else
        cout<<"Successfully invoking applyRealFFT for "<<stage<<"-stage FFT calculation."<<endl;
    cout<<endl<<"Applying forward FFT on the test sequence using the fftpack..."<<endl;
    // Prepare Python/C environment.
    setenv("PYTHONPATH", "./", 1);
    Py_SetProgramName((char *)argv[0]);
    Py_Initialize();
    // Import numpy.
    init_numpy();
    // Import scipy.
    PyObject* pModule = PyImport_ImportModule("scipy");
    CHECK_NULL_ASSERT(pModule);
    PyObject* pDict = PyModule_GetDict(pModule);
    CHECK_NULL_ASSERT(pDict);
    // Find the fft function.
    PyObject* pFunc = PyDict_GetItemString(pDict, "fft");
    CHECK_NULL_ASSERT(pFunc);
    // Prepare the input data for invoking the fft function.
    PyArrayObject *np_arg;
    npy_intp dims[1];
    dims[0] = TEST_FFT_LENGTH;
    np_arg = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNewFromData(1, &dims[0], NPY_FLOAT, 
        reinterpret_cast<void*>(fftpackInArray)));
    PyObject *pArgs;
    pArgs = PyTuple_New(1);
    PyTuple_SetItem(pArgs, 0, reinterpret_cast<PyObject*>(np_arg));
    // Call fftpack for calculating FFT.
    PyObject* PyResult =  PyObject_CallObject(pFunc, pArgs);
    PyArrayObject* pContArray = PyArray_GETCONTIGUOUS((PyArrayObject*)PyResult);
    if (NULL != pContArray) {
        cout<<"Got a result from fftpack."<<endl;
        // Copy the result to conventional arrays for easy comparison.
        int ndim = PyArray_NDIM(pContArray);
        long int *dims = PyArray_DIMS(pContArray);
        uint8_t *pData = (uint8_t *)PyArray_DATA(pContArray);
        for(int i = 0; i < ndim; i++){
            for(int j = 0; j < dims[i]; j++){
                Py_complex * c = (Py_complex *)(pData + i * PyArray_STRIDE(pContArray,1) + j * PyArray_STRIDE(pContArray,0));
                fftpackRealOutArray[j]=(float)c->real;
                fftpackImagOutArray[j]=(float)c->imag;
            }
        }
        cout<<endl<<"Now comparing results..."<<endl;
        cout<<endl<<"real parts(fftpack,applyRealFFT):"<<endl;
        for(int i=0;i<(TEST_FFT_LENGTH>>1)+1;i++) {
            if(i==(TEST_FFT_LENGTH>>1)) {
                // The real part of Nyquist bin will be stored in the image part of DC in applyRealFFT function.
                // The imaginary part of DC and Nyquist bins will be zero in real signal.
                cout<<"("<<fftpackRealOutArray[i]<<","<<realFftInOutArray[1]<<"),abs diff.="<<abs(fftpackRealOutArray[i] -  realFftInOutArray[1])<<endl;
            }
            else
                cout<<"("<<fftpackRealOutArray[i]<<","<<realFftInOutArray[2*i]<<"),abs diff.="<< abs(fftpackRealOutArray[i] -  realFftInOutArray[2*i])<<endl;
        }
        cout<<endl<<"image parts(fftpack,applyRealFFT):"<<endl;
        for(int i=0;i<(TEST_FFT_LENGTH>>1)+1;i++) {
            if(i==0 || i==(TEST_FFT_LENGTH>>1))
               cout<<"("<<fftpackImagOutArray[i]<<","<<0.0<<"),abs diff.="<< abs(fftpackImagOutArray[i] - 0.0)<<endl;
            else
                cout<<"("<<fftpackImagOutArray[i]<<","<<realFftInOutArray[2*i+1]<<"),abs diff.="<< abs(fftpackImagOutArray[i] - realFftInOutArray[2*i+1])<<endl;
        }
    }
    else {
         cout<<"Failed to get results using fftpack. Please check your environment."<<endl;
    }
    fftDestroy(config);
    Py_DECREF(pArgs);
    if (pContArray != NULL)
        Py_DECREF(pContArray);
    if (PyResult != NULL)
        Py_DECREF(PyResult);
    Py_DECREF(pModule);
    Py_Finalize();
    return 0;
}

