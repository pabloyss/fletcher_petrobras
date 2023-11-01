#ifndef __CUDA_DEFINES
#define __CUDA_DEFINES

#define restrict __restrict__
#define BSIZE_X 32
#define BSIZE_Y 16
#define NPOP 4
#define TOTAL_X (BSIZE_X+2*NPOP)
#define TOTAL_Y (BSIZE_Y+2*NPOP)

#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <stdio.h>

#define HIP_CALL(call) do{      \
   const hipError_t err=call;         \
   if (err != hipSuccess)       \
   {                             \
     fprintf(stderr, "CUDA ERROR: %s on %s:%d\n", hipGetErrorString(err), __FILE__, __LINE__);\
     exit(1);                    \
   }}while(0)

#endif
