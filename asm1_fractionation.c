/*
 * asm1 fractionation model block 
 * Copyright: Xavier Flores-Alsina, Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 *            Krist V Gernaey, DTU. Denmark
 */

#define S_FUNCTION_NAME asm1_fractionation

#include "simstruc.h"
#include <math.h>

#define ASM1_PARS   ssGetArg(S,0)       /* definition of the ASM1 kinetic and stoichiometric coefficients            */
#define ASM1_FRACTIONS   ssGetArg(S,1)  /* definition of the ASM1 kinetic and stoichiometric coefficients            */

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 0);   /* number of continuous states                                                      */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states                                                        */
    ssSetNumInputs(        S, 11);  /* number of inputs                            (u[0] = CODsol_HH; u[1] = CODpart_HH, u[2] = SNH_HH; u[3] = TKN_HH; u[4] = NOX_HH; ; u[5] = CODsol_Ind; u[6] = CODpart_Ind, u[7] = SNH_Ind; u[8] = TKN_Ind; u[9] = NOX_Ind;, u[10] = flow data;   )       */
    ssSetNumOutputs(       S, 13);  /* number of outputs                           (y[0] = SI; y[1] = SS; y[2] = XI; y[3] = XS; y[4] = XBH; y[5] = XBA; y[6] = XU; y[7] = SO; y[8] = SNO; y[9] = SNH; y [10] = SND; y[11] = XND; y[12] = SALK ) */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag                                                          */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                                                           */
    ssSetNumSFcnParams(    S, 2);   /* number of input arguments                    (ASM1_PARS, ASM1_FRACTIONS)         */
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

double i_XB, i_XP;
double SI_fr,SS_fr,XI_fr,XS_fr,XBH_fr,XBA_fr,XP_fr,SO_fr,SNO_fr,SNH_fr,SND_fr,XND_fr,SALK_fr;
double SI_cst,SS_cst,XI_cst,XS_cst,XBH_cst,XBA_cst,XP_cst,SO_cst,SNO_cst,SNH_cst,SND_cst,XND_cst,SALK_cst;
double CODsolflux, Norg;
    

/* load the stoichiometric parameters that are needed for the fractionation */
i_XB = mxGetPr(ASM1_PARS)[17];
i_XP = mxGetPr(ASM1_PARS)[18];

/* load the ASM1 component fractions */
SI_cst  =  mxGetPr(ASM1_FRACTIONS)[0];
XI_fr   =  mxGetPr(ASM1_FRACTIONS)[1];
XS_fr   =  mxGetPr(ASM1_FRACTIONS)[2];
XBH_fr  =  mxGetPr(ASM1_FRACTIONS)[3];
XBA_fr  =  mxGetPr(ASM1_FRACTIONS)[4];
XP_fr   =  mxGetPr(ASM1_FRACTIONS)[5];
SNO_fr  =  mxGetPr(ASM1_FRACTIONS)[6];
SNH_fr  =  mxGetPr(ASM1_FRACTIONS)[7];
SND_fr  =  mxGetPr(ASM1_FRACTIONS)[8];
XND_fr  =  mxGetPr(ASM1_FRACTIONS)[9];

CODsolflux=u[0]+u[5];
y[0] = CODsolflux/u[10];
y[1] = 0.0;
y[2] = XI_fr*(u[1]+u[6])/u[10];                     /* XI defined as a fraction of the particulate COD */
y[3] = XS_fr*(u[1]+u[6])/u[10];                     /* XS defined as a fraction of the particulate COD */
y[4] = XBH_fr*(u[1]+u[6])/u[10];                    /* XBH defined as a fraction of the particulate COD */
y[5] = XBA_fr*(u[1]+u[6])/u[10];                    /* XBA defined as a fraction of the particulate COD */
y[6] = XP_fr*(u[1]+u[6])/u[10];                     /* XP defined as a fraction of the particulate COD */
y[7] = 0;                                           /* SO is assumed to be zero */
y[8] = SNO_fr*(u[4]+u[9])/u[10];                   /* SNO defined as the result of combining the inorganic N fluxes */
y[9] = SNH_fr*(u[2]+u[7])/u[10];                    /* SNH defined as the result of combining the ammonium N fluxes */
    
Norg=(u[3]+u[8]-u[2]-u[7])/u[10]-i_XB*(y[4]+y[5])-i_XP*y[6];
                                                    /* Calculate the remaining organic N concentration, taking into account that some organic N is 
                                                     present in the active biomass and in decay products of biomass */
if (Norg>=0.0){
    y[10] = SND_fr*Norg;                            /* SND defined as a fraction of the remaining organic N */
    y[11] = XND_fr*Norg;                            /* XND defined as a fraction of the remaining organic N */
    }
else if (Norg<0.0){           
    y[10] = 0;
    y[11] = 0;
    }
y[12] = 7.0;                                        /* SALK is assumed to be constant and not limiting for the processes in the WWTP */   

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
