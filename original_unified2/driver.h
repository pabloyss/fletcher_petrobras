#ifndef __driver_h__
#define __driver_h__

#ifdef __cplusplus
extern "C"
{
#endif

	void DRIVER_Initialize(const int sx, const int sy, const int sz, const int bord,
						   float dx, float dy, float dz, float dt,
						   float *restrict vpz, float *restrict vsv, float *restrict epsilon, float *restrict delta,
						   float *restrict phi, float *restrict theta,
						   float *restrict pp, float *restrict pc, float *restrict qp, float *restrict qc);

	void DRIVER_Finalize(float *restrict vpz, float *restrict vsv);

	void DRIVER_Propagate(const int sx, const int sy, const int sz, const int bord,
						  const float dx, const float dy, const float dz, const float dt, const int it,
						  float *pp, float *pc, float *qp, float *qc);

	void DRIVER_Update_pointers(const int sx, const int sy, const int sz, float *pc);

	void DRIVER_InsertSource(float dt, int it, int iSource, float *p, float *q, float src);

#ifdef __cplusplus
}
#endif
#endif
