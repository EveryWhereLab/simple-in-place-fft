# Simple in-place FFT

This repository contains simple FFT/IFFT functions.

This implementation is based on the algorithms described in chapter 12 of the numerical recipes in C book.
There are slight differences between the original implementation and this one in the following aspects:
1. The array index starts from zero in this implementation.
2. A init. function is added to create the initial trigonometric tables for later use.

## To build the library and try a test, run:

```
mkdir build && cd build
cmake ..
cmake --build .
./fft_ifft_test
```

## To try comparison tests:
python-dev and scipy packages must be installed in the environment, and run:
```
mkdir build && cd build
cmake -DBUILD_COMPARISON_TEST=ON ..
cmake --build .
./cfft_numerical_comparison
./rfft_numerical_comparison
./fft_ifft_test
```

