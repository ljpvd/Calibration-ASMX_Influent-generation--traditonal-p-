/*
 * asm1 combiner model block 
 * Copyright: Xavier Flores-Alsina, Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 *            Krist V Gernaey, DTU. Denmark
 */

#define S_FUNCTION_NAME asm1_combiner

#include "simstruc.h"


/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 0);   /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 42);  /* number of inputs  (u[0-20] = ASM1 state variables + TSS +  Q +  T + 3 S dummy states + 2 X dummy states; u[21-41] = ASM1 state variables + TSS +  Q +  T + 3 S dummy states + 2 X dummy states; */
    ssSetNumOutputs(       S, 21);  /* number of outputs (y[0-20] = ASM1 state variables + TSS +  Q +  T + 3 S dummy states + 2 X dummy states  */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 0);   /* number of input arguments             */
    ssSetNumRWork(         S, 0);   /* number of real work vector elements   */
    ssSetNumIWork(         S, 0);   /* number of integer work vector elements*/
    ssSetNumPWork(         S, 0);   /* number of pointer work vector elements*/
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
  y[0]=(u[0]*u[14] + u[21]*u[35])/(u[14]+u[35]);     /* ASM1-state: SI */
  y[1]=(u[1]*u[14] + u[22]*u[35])/(u[14]+u[35]);     /* ASM1-state: SS */
  y[2]=(u[2]*u[14] + u[23]*u[35])/(u[14]+u[35]);     /* ASM1-state: XI */
  y[3]=(u[3]*u[14] + u[24]*u[35])/(u[14]+u[35]);     /* ASM1-state: XS */
  y[4]=(u[4]*u[14] + u[25]*u[35])/(u[14]+u[35]);     /* ASM1-state: XBH */
  y[5]=(u[5]*u[14] + u[26]*u[35])/(u[14]+u[35]);     /* ASM1-state: XBA */
  y[6]=(u[6]*u[14] + u[27]*u[35])/(u[14]+u[35]);     /* ASM1-state: XP */
  y[7]=(u[7]*u[14] + u[28]*u[35])/(u[14]+u[35]);     /* ASM1-state: SO */
  y[8]=(u[8]*u[14] + u[29]*u[35])/(u[14]+u[35]);     /* ASM1-state: SNO */
  y[9]=(u[9]*u[14] + u[30]*u[35])/(u[14]+u[35]);     /* ASM1-state: SNH */
  y[10]=(u[10]*u[14] + u[31]*u[35])/(u[14]+u[35]);   /* ASM1-state: SND */
  y[11]=(u[11]*u[14] + u[32]*u[35])/(u[14]+u[35]);   /* ASM1-state: XND */
  y[12]=(u[12]*u[14] + u[33]*u[35])/(u[14]+u[35]);   /* ASM1-state: SALK */
  y[13]=(u[13]*u[14] + u[34]*u[35])/(u[14]+u[35]);   /* TSS */
  y[14]=(u[14]+u[35]);                               /* Q */
  y[15]=(u[15]*u[14] + u[36]*u[35])/(u[14]+u[35]);   /* T */
  y[16]=(u[16]*u[14] + u[37]*u[35])/(u[14]+u[35]);   /* S Dummy state 1 */
  y[17]=(u[17]*u[14] + u[38]*u[35])/(u[14]+u[35]);   /* S Dummy state 2 */
  y[18]=(u[18]*u[14] + u[39]*u[35])/(u[14]+u[35]);   /* S Dummy state 3 */
  y[19]=(u[19]*u[14] + u[40]*u[35])/(u[14]+u[35]);   /* X Dummy state 4 */
  y[20]=(u[20]*u[14] + u[41]*u[35])/(u[14]+u[35]);   /* X Dummy state 5 */

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
