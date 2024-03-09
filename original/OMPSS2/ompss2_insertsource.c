#include "ompss2_insertsource.h"


// InsertSource: compute and insert source value at index iSource of arrays p and q

#pragma oss task inout (p[iSource], q[iSource])
void OPENMP_InsertSource(float dt, int it, int iSource, 
			  float *p, float*q, float src) {
  {
//      src  = Source(dt, it);
    p[iSource]+=src;
     q[iSource]+=src;
  }
}
