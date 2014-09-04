/*
 * first flush model block 
 * Copyright: Xavier Flores-Alsina, Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 *            Krist V Gernaey, DTU. Denmark
 */


#define S_FUNCTION_NAME firstflush_ASM1

#include "simstruc.h"
#include <math.h>

#define ASM1_XINIT   ssGetArg(S,0)
#define SSPARS	    ssGetArg(S,1)

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 9);   /* number of continuous states                      (x[0] = TSS; x[1] = XI; x[2] = XS; x[3] = XBH; x[4] = XBHA; x[5] = XU; x[6] = XND;  x[7] = XD1; x[8] = XD2;                              */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states                                                        */
    ssSetNumInputs(        S, 10);  /* number of inputs                                 (u[0] = Q; u[1] = TSS, u[2-10] = XI, XS, XBH, XBA, XU, XND, XD1, XD2)       */
    ssSetNumOutputs(       S, 19);  /* number of outputs                                (y[0] = Q; y[1 -9] = TSS, XI, XS, XBH, XBA, XU, XND, XD1, XD2; y[10-18] = mass of TSS, XI, XS, XBH, XBA, XU, XND, XD1, XD2) */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag                                                          */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                                                           */
    ssSetNumSFcnParams(    S, 2);   /* number of input arguments                        (ASM1_PARS, SSPAR)             */
    ssSetNumRWork(         S, 0);   /* number of real work vector elements                                              */
    ssSetNumIWork(         S, 0);   /* number of integer work vector elements                                           */
    ssSetNumPWork(         S, 0);   /* number of pointer work vector elements                                           */
    
  
}

/*
 * mdlInitializeSampleTimes - initialize the sample times array
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

/*
 * mdlInitializeConditions - initialize the states
 */
static void mdlInitializeConditions(double *x0, SimStruct *S)
{
int i;

for (i = 0; i < 9; i++) {
   x0[i] = mxGetPr(ASM1_XINIT)[i];
}
}

/*
 * mdlOutputs - compute the outputs
 */

static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{

int i;
double M_Max, Q_lim, n, Ff, Ff_factor_out, Hill;

M_Max = mxGetPr(SSPARS)[0];
Q_lim = mxGetPr(SSPARS)[1];
n     = mxGetPr(SSPARS)[2];
Ff    = mxGetPr(SSPARS)[3];

for (i = 0; i < 7; i++) {
  if (x[i] < 0.0)
    x[i] = 0.0;
}

Ff_factor_out=(x[0]/M_Max);
Hill=pow(u[0]/10000,n)/(pow(Q_lim/10000,n)+pow(u[0]/10000,n));

  y[0]=u[0];                /* Flow rate */
  y[1]=u[1]*Ff_factor_out+Hill*x[0]*Ff/u[0];
  y[2]=u[2]*Ff_factor_out+Hill*x[1]*Ff/u[0];
  y[3]=u[3]*Ff_factor_out+Hill*x[2]*Ff/u[0];
  y[4]=u[4]*Ff_factor_out+Hill*x[3]*Ff/u[0];
  y[5]=u[5]*Ff_factor_out+Hill*x[4]*Ff/u[0];
  y[6]=u[6]*Ff_factor_out+Hill*x[5]*Ff/u[0];
  y[7]=u[7]*Ff_factor_out+Hill*x[6]*Ff/u[0];
  y[8]=u[8]+Hill*x[7]*Ff/u[0];
  y[9]=u[9]+Hill*x[8]*Ff/u[0];
  
  y[10]=x[0];                /* Internal state, mass of sediment (TSS) stored in the sewer (kg COD)*/
  y[11]=x[1]/1000;           /* Internal state, mass of sediment XI stored in the sewer (kg COD)*/
  y[12]=x[2]/1000;          /* Internal state, mass of sediment XS stored in the sewer (kg COD)*/
  y[13]=x[3]/1000;          /* Internal state, mass of sediment XBH stored in the sewer (kg COD)*/
  y[14]=x[4]/1000;          /* Internal state, mass of sediment XBA stored in the sewer (kg COD)*/
  y[15]=x[5]/1000;          /* Internal state, mass of sediment XP stored in the sewer (kg COD)*/
  y[16]=x[6]/1000;          /* Internal state, mass of sediment XND stored in the sewer (kg N)*/
  y[17]=x[7]/1000;          /* Internal state, mass of sediment XD1 stored in the sewer (kg )*/
  y[18]=x[8]/1000;          /* Internal state, mass of sediment XFD2 stored in the sewer (kg )*/
}

/*
 * mdlUpdate - perform action at major integration time step
 */

static void mdlUpdate(double *x, double *u, SimStruct *S, int tid)
{
}

/*
 * mdlDerivatives - compute the derivatives
 */
static void mdlDerivatives(double *dx, double *x, double *u, SimStruct *S, int tid)
{

double M_Max, Q_lim, n, Ff, Ff_factor, Hill;

M_Max = mxGetPr(SSPARS)[0];
Q_lim = mxGetPr(SSPARS)[1];
n     = mxGetPr(SSPARS)[2];
Ff    = mxGetPr(SSPARS)[3];

Ff_factor=(1-x[0]/M_Max);
Hill=pow(u[0]/10000,n)/(pow(Q_lim/10000,n)+pow(u[0]/10000,n));

dx[0] = u[0]*u[1]*Ff_factor-Hill*x[0]*Ff;  /* TSS balance */
dx[1] = u[0]*u[2]*Ff_factor-Hill*x[1]*Ff;  /* XI balance */
dx[2] = u[0]*u[3]*Ff_factor-Hill*x[2]*Ff;  /* XS balance */
dx[3] = u[0]*u[4]*Ff_factor-Hill*x[3]*Ff;  /* XBH balance */
dx[4] = u[0]*u[5]*Ff_factor-Hill*x[4]*Ff;  /* XBA balance */
dx[5] = u[0]*u[6]*Ff_factor-Hill*x[5]*Ff;  /* XP balance */
dx[6] = u[0]*u[7]*Ff_factor-Hill*x[6]*Ff;  /* XND balance */
dx[7] = u[0]*u[8]-Hill*x[7]*Ff;  /* XD1 balance */
dx[8] = u[0]*u[9]-Hill*x[8]*Ff;  /* XD2 balance */

}

/*
 * mdlTerminate - called when the simulation is terminated.
 */
static void mdlTerminate(SimStruct *S)
{
}

#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
