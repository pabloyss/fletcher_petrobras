#ifndef _OPENACC_PROPAGATE
#define _OPENACC_PROPAGATE

// Propagate: using Fletcher's equations, propagate waves one dt,
//            either forward or backward in time


void OPENACC_PropagateDer1Der1LM(int sx, int sy, int sz, int bord,
	       float dx, float dy, float dz, float dt, int it, 
	       float * restrict pp, float * restrict pc, float * restrict qp, float * restrict qc);

#endif
