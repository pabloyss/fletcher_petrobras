#include "utils.h"
#include "source.h"
#include "driver.h"
#include "fletcher.h"
#include "walltime.h"
#include "model.h"
#ifdef PAPI
#include "ModPAPI.h"
#endif

#define MODEL_GLOBALVARS
#include "precomp.h"
#undef MODEL_GLOBALVARS


void ReportProblemSizeCSV(const int sx, const int sy, const int sz,
			  const int bord, const int st, 
			  FILE *f){
  fprintf(f,
	  "sx; %d; sy; %d; sz; %d; bord; %d;  st; %d; \n",
	  sx, sy, sz, bord, st);
}

void ReportMetricsCSV(double walltime, double MSamples,
		      long HWM, char *HWMUnit, FILE *f){
  fprintf(f,
	  "walltime; %lf; MSamples; %lf; HWM;  %ld; HWMUnit;  %s;\n",
	  walltime, MSamples, HWM, HWMUnit);
}

void Model(const int st, const int iSource, const float dtOutput, SlicePtr sPtr,
           const int sx, const int sy, const int sz, const int bord,
           const float dx, const float dy, const float dz, const float dt, const int it, 
	   float * restrict pp, float * restrict pc, float * restrict qp, float * restrict qc,
	   float * restrict vpz, float * restrict vsv, float * restrict epsilon, float * restrict delta,
	   float * restrict phi, float * restrict theta)
{

  float tSim=0.0;
  int nOut=1;
  float tOut=nOut*dtOutput;

  const long samplesPropagate=(long)(sx-2*bord)*(long)(sy-2*bord)*(long)(sz-2*bord);
  const long totalSamples=samplesPropagate*(long)st;

#ifdef PAPI
  long long values[NCOUNTERS];
  long long ThisValues[NCOUNTERS];
  for (int i=0; i<NCOUNTERS; i++) {
    values[i]=0LL;
    ThisValues[i]=0LL;
  }

  const int eventset=InitPAPI_CreateCounters();
#endif

#define MODEL_INITIALIZE
#ifdef MODEL_GLOBALVARS
    float *ch1dxx=NULL;  // isotropy simetry deep angle
float *ch1dyy=NULL;  // isotropy simetry deep angle
float *ch1dzz=NULL;  // isotropy simetry deep angle
float *ch1dxy=NULL;  // isotropy simetry deep angle
float *ch1dyz=NULL;  // isotropy simetry deep angle
float *ch1dxz=NULL;  // isotropy simetry deep angle
float *v2px=NULL;  // coeficient of H2(p)
float *v2pz=NULL;  // coeficient of H1(q)
float *v2sz=NULL;  // coeficient of H1(p-q) and H2(p-q)
float *v2pn=NULL;  // coeficient of H2(p)
#endif

#ifdef MODEL_INITIALIZE
// Precalcula campos abaixo

// coeficients of derivatives at H1 operator
    ch1dxx = (float *) malloc(sx*sy*sz*sizeof(float));
    ch1dyy = (float *) malloc(sx*sy*sz*sizeof(float));
    ch1dzz = (float *) malloc(sx*sy*sz*sizeof(float));
    ch1dxy = (float *) malloc(sx*sy*sz*sizeof(float));
    ch1dyz = (float *) malloc(sx*sy*sz*sizeof(float));
    ch1dxz = (float *) malloc(sx*sy*sz*sizeof(float));


// coeficients of H1 and H2 at PDEs
    v2px = (float *) malloc(sx*sy*sz*sizeof(float));
    v2pz = (float *) malloc(sx*sy*sz*sizeof(float));
    v2sz = (float *) malloc(sx*sy*sz*sizeof(float));
    v2pn = (float *) malloc(sx*sy*sz*sizeof(float));

#pragma oss taskloop grainsize(284 * 284) in(theta[i],phi[i]) out(ch1dxx[i],ch1dyy[i],ch1dzz[i],ch1dxy[i],ch1dyz[i],ch1dxz[i])
    for (int i=0; i<sx*sy*sz; i++) {
        float sinTheta=sin(theta[i]);
        float cosTheta=cos(theta[i]);
        float sin2Theta=sin(2.0*theta[i]);
        float sinPhi=sin(phi[i]);
        float cosPhi=cos(phi[i]);
        float sin2Phi=sin(2.0*phi[i]);
        ch1dxx[i]=sinTheta*sinTheta * cosPhi*cosPhi;
        ch1dyy[i]=sinTheta*sinTheta * sinPhi*sinPhi;
        ch1dzz[i]=cosTheta*cosTheta;
        ch1dxy[i]=sinTheta*sinTheta * sin2Phi;
        ch1dyz[i]=sin2Theta         * sinPhi;
        ch1dxz[i]=sin2Theta         * cosPhi;
    }
#ifdef _DUMP
//    {
//        const int iPrint=ind(bord+1,bord+1,bord+1);
//        printf("ch1dxx=%f; ch1dyy=%f; ch1dzz=%f; ch1dxy=%f; ch1dxz=%f; ch1dyz=%f\n",
//               ch1dxx[iPrint], ch1dyy[iPrint], ch1dzz[iPrint], ch1dxy[iPrint], ch1dxz[iPrint], ch1dyz[iPrint]);
//    }
#endif


#pragma oss taskloop grainsize(284 * 284) in(vsv[i],vpz[i]) out(v2sz[i],v2pz[i],v2px[i],v2pn[i])
    for (int i=0; i<sx*sy*sz; i++){
        v2sz[i]=vsv[i]*vsv[i];
        v2pz[i]=vpz[i]*vpz[i];
        v2px[i]=v2pz[i]*(1.0+2.0*epsilon[i]);
        v2pn[i]=v2pz[i]*(1.0+2.0*delta[i]);
    }

#ifdef _DUMP
    {
//
//        const int iPrint=ind(bord+1,bord+1,bord+1);
//        printf("vsv=%e; vpz=%e, v2pz=%e\n",
//               vsv[iPrint], vpz[iPrint], v2pz[iPrint]);
//
//        printf("v2sz=%e; v2pz=%e, v2px=%e, v2pn=%e\n",
//               v2sz[iPrint], v2pz[iPrint], v2px[iPrint], v2pn[iPrint]);
    }

#endif

#endif
#undef MODEL_INITIALIZE

    // DRIVER_Initialize initialize target, allocate data etc
  DRIVER_Initialize(sx,   sy,   sz,   bord,
		      dx,  dy,  dz,  dt,
		      vpz,    vsv,    epsilon,    delta,
		      phi,    theta,
		      pp,    pc,    qp,    qc);

  double walltime=0.0;

//#pragma oss taskwait
// Run N steps

    const double t0=wtime();

  for (int it=1; it<=st; it++) {

              // IN: dt
          float src  = Source(dt, it-1);
//#pragma oss task in (dt,src) out (pc, qc)
//      printf("\n src: %d \n", iSource);
      DRIVER_InsertSource(dt,it-1,iSource,pc,qc,src);

//      DRIVER_InsertSource(dt,it-1,iSource,pc,qc,src);
    // OUT: pc, qc

#ifdef PAPI
    StartCounters(eventset);
#endif

    // Propaga na grid

    // IN: pp,pc,qp,qc
      DRIVER_Propagate(  sx,   sy,   sz,   bord,
		       dx,   dy,   dz,   dt,   it,
		       pp,    pc,    qp,    qc);
    // OUT: pp,qp


    // IN: pp,pc,qp,qc
      SwapArrays(&pp, &pc, &qp, &qc);
      // OUT: pp,pc,qp,qc


#ifdef PAPI
    StopReadCounters(eventset, ThisValues);
    for (int i=0; i<NCOUNTERS; i++) {
      values[i]+=ThisValues[i];
    }
#endif

    tSim=it*dt;
    if (tSim >= tOut) {

        // não faz nada
      DRIVER_Update_pointers(sx,sy,sz,pc);

      // IN: sx,sy,pc,sPtr
#pragma oss task in (pc[0:22906303]) inout(sPtr)
        DumpSliceFile(sx,sy,sz,pc,sPtr);
      // OUT: sPtr

      tOut=(++nOut)*dtOutput;

#ifdef _DUMP
      // nada
        DRIVER_Update_pointers(sx,sy,sz,pc);
      // DumpSliceSummary(sx,sy,sz,sPtr,dt,it,pc,src);
#endif
    }
  }
#pragma oss taskwait

    walltime+=wtime()-t0;

  // get HWM data

#define MEGA 1.0e-6
#define GIGA 1.0e-9
  const char StringHWM[6]="VmHWM";
  char line[256], title[12],HWMUnit[8];
  const long HWM;
  const double MSamples=(MEGA*(double)totalSamples)/walltime;
  
  FILE *fp=fopen("/proc/self/status","r");
  while (fgets(line, 256, fp) != NULL){
    if (strncmp(line, StringHWM, 5) == 0) {
      sscanf(line+6,"%ld %s", &HWM, HWMUnit);
      break;
    }
  }
  fclose(fp);

  // Dump Execution Metrics
  
  printf ("Execution time (s) is %lf\n", walltime);
  printf ("MSamples/s %.0lf\n", MSamples);
  printf ("Memory High Water Mark is %ld %s\n",HWM, HWMUnit);

  // Dump Execution Metrics in CSV
  
  FILE *fr=NULL;
  const char fName[]="Report.csv";
  fr=fopen(fName,"w");

  // report problem size

  ReportProblemSizeCSV(sx, sy, sz,
		       bord, st, 
		       fr);

  // report collected metrics

  ReportMetricsCSV(walltime, MSamples,
		   HWM, HWMUnit, fr);
  
  // report PAPI metrics

#ifdef PAPI
  ReportRawCountersCSV (values, fr);
#endif
  
  fclose(fr);

  fflush(stdout);

  // DRIVER_Finalize deallocate data, clean-up things etc 
  DRIVER_Finalize();

}

