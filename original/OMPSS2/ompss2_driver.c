#include "../driver.h"
#include "ompss2_propagate.h"
#include "ompss2_insertsource.h"
#include "../sample.h"

void DRIVER_Initialize(const int sx, const int sy, const int sz, const int bord,
		       float dx, float dy, float dz, float dt,
		       float * restrict vpz, float * restrict vsv, float * restrict epsilon, float * restrict delta,
		       float * restrict phi, float * restrict theta,
		       float * restrict pp, float * restrict pc, float * restrict qp, float * restrict qc)
{
//    printf("Driver Initializa\n");

}


void DRIVER_Finalize()
{
}


void DRIVER_Update_pointers(const int sx, const int sy, const int sz, float *pc)
{
}

#pragma oss task in (pc[bord;sx*sy*sz],qc[bord:(sx*sy*sz)-1]) inout (pp[bord:(sx*sy*sz)-1], qp[bord:(sx*sy*sz)-1])
void DRIVER_Propagate(const int sx, const int sy, const int sz, const int bord,
	       const float dx, const float dy, const float dz, const float dt, const int it, 
	       float * restrict pp, float * restrict pc, float * restrict qp, float * restrict qc)
{

	OPENMP_Propagate (  sx,   sy,   sz,   bord,
                                  dx,   dy,   dz,   dt,   it,
                                  pp,   pc,   qp,   qc);

}

#pragma oss task inout (p[iSource], q[iSource])
void DRIVER_InsertSource(float dt, int it, int iSource, float *p, float*q, float src)
{
    OPENMP_InsertSource(dt,it,iSource,p,q,src);
}

