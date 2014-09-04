/*
 * asm1 fractionation model block 
 * Copyright: Xavier Flores-Alsina, Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 *            Krist V Gernaey, DTU. Denmark
 */

#define S_FUNCTION_NAME asmX_fractionation

#include "simstruc.h"
#include <math.h>


/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 0);   /* number of continuous states                                                      */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states                                                        */
    ssSetNumInputs(        S, 5);  /* number of inputs              (u[0] = Cli, u[1]= Ccj, u[2] = Csl, u[3] Csl_i, u[4] = flow */
    ssSetNumOutputs(       S, 4);  /* number of outputs                                                                */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag                                                          */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                                                           */
    ssSetNumSFcnParams(    S, 0);   /* number of input arguments                                                        */
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
}

/*
 * mdlOutputs - compute the outputs
 */

static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{

double Cli_flux, Ccj_flux, Csl_flux, Csl_i_flux;


Cli_flux = u[0];
Ccj_flux = u[1];
Csl_flux = u[2];
Csl_i_flux = u[3];


y[0] = Cli_flux / u[4];
y[1] = Ccj_flux / u[4];
y[2] = Csl_flux / u[4];
y[3] = Csl_i_flux / u[4];


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
