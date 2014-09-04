/*
 * Soil model block 
 * Copyright: Xavier Flores-Alsina, Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 *            Krist V Gernaey, DTU. Denmark
 */

#define S_FUNCTION_NAME unisoilmodel

#include "simstruc.h"
#include <math.h>

#define XINIT   ssGetArg(S,0)       /* definition of initial values of the soil model state variables (h1)              */
#define PAR	ssGetArg(S,1)           /* definition of the set of parmater values to include in the soil model            */

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 1);   /* number of continuous states                  (x[0] = h1)                         */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states                                                        */
    ssSetNumInputs(        S, 2);   /* number of inputs                             (u[0] = qin 1 - upstream infiltration and   u[1] = qin2 - rain generation)       */
    ssSetNumOutputs(       S, 4);   /* number of outputs                            (y[0] = qout1 - infiltration to the sewers, y[1] = overflow rate,           y[2] = quout2 - infiltration to aquifers and    y[3] = h1 level)*/
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag                                                          */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                                                           */
    ssSetNumSFcnParams(    S, 2);   /* number of input arguments                    (XINIT, PAR)                        */
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
   x0[0] = mxGetPr(XINIT)[0];     /* Initial h1 value */
}

/*
 * mdlOutputs - compute the outputs
 */

static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{

double HMAX, HINV, A, K, Kinf, Kdown;

HMAX  = mxGetPr(PAR)[0];   /* Maximum water level in tank */
HINV  = mxGetPr(PAR)[1];   /* Invert level in tank, i.e. water level corresponding to bottom of sewer pipes */
A     = mxGetPr(PAR)[2];   /* Surface area */
K     = mxGetPr(PAR)[3];   /* Parameter, permeability of soil for water penetration */
Kinf  = mxGetPr(PAR)[4];   /* Infiltration gain, can be an indication of the quality of the sewer system pipes */
Kdown = mxGetPr(PAR)[5];   /* Gain to adjust for flow rate to downstream aquifer */

if (x[0]<HINV)
    y[0]=0.0;
else
    y[0]=Kinf*sqrt(x[0]-HINV); /* Flow rate of infiltration water into sewer pipe */

if (u[0]<K*A)
    y[1]=0.0;
else
    y[1]=u[0]-K*A;         /* Overflow rate, i.e. rain water that is not penetrating in the soil */

y[2]=Kdown*x[0];           /* Flow rate of water to downstream aquifers */

y[3]=x[0];                 /* level in tank*/

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

double HMAX, HINV, A, K, Kinf, Kdown;
double qin1, qin2, qout1, qout2, Hinf;

HMAX  = mxGetPr(PAR)[0];   /* Maximum water level in tank */
HINV  = mxGetPr(PAR)[1];   /* Invert level in tank, i.e. water level corresponding to bottom of sewer pipes */
A     = mxGetPr(PAR)[2];   /* Surface area */
K     = mxGetPr(PAR)[3];   /* Parameter, permeability of soil for water penetration */
Kinf  = mxGetPr(PAR)[4];   /* Infiltration gain, can be an indication of the quality of the sewer system pipes */
Kdown = mxGetPr(PAR)[5];   /* Gain to adjust for flow rate to downstream aquifer */


// 1 condition qin1//
if (u[0]>K*A)	/* u[0] = rain flow rate as input to this model block, limited to a maximum (top soil permeability for rain) */
    qin1=K*A;
else
    qin1=u[0];

// 2 condition qin2//
if (x[0]>=HMAX)
    qout2 = 0;
else
    qout2= u[1]; /* u[1] = flow rate from other (upflow) aquifers, as input to this model block, limited to a maximum */

// 3 condition qout2//
if (x[0]>HINV)
    Hinf=x[0]-HINV;
else
    Hinf=0;

// 4 condition : avoid negative numbers //
if (x[0] < 0.0)
     x[0] = 0.0;
   else
     x[0]= x[0];



/* This is the only balance equation */

dx[0]=(1/A)*(qin1+u[1]-Kinf*sqrt(Hinf)-Kdown*x[0]);

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
