/*
 * sewer_asm1 model block 
 * Copyright: Xavier Flores-Alsina, Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 *            Krist V Gernaey, DTU. Denmark
 */

#define S_FUNCTION_NAME sewer_asm1

#include "simstruc.h"
#include <math.h>

#define ASM1_SEWERINIT ssGetArg(S,0)
#define ASM1_PARS ssGetArg(S,1)
#define SPAR ssGetArg(S,2)
#define SOSAT	ssGetArg(S,3)
#define MPOL_PAR ssGetArg(S,4)
#define KLa ssGetArg(S,5)



/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 19);  /* number of continuous states           (x[0-18] = ASM1 state variables +3 S dummy states + 2 X dummy states + h */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states */            
    ssSetNumInputs(        S, 22);  /* number of inputs (corresponding to outputs) (u[0-20] = ASM1 state variables + TSS +  Q +  T + 3 S dummy states + 2 X dummy states + h */
    ssSetNumOutputs(       S, 22);  /* number of outputs                           (y[0-20] = ASM1 state variables + TSS +  Q +  T + 3 S dummy states + 2 X dummy states + h */  
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 6);   /* number of input arguments             */
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
int i;

for (i = 0; i < 19; i++) {
   x0[i] = mxGetPr(ASM1_SEWERINIT)[i];
}
}

/*
 * mdlOutputs - compute the outputs
 */

static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{
  double X_I2TSS, X_S2TSS, X_BH2TSS, X_BA2TSS, X_P2TSS, A, C, Hmin, vol;
  
  int i;

  X_I2TSS = mxGetPr(ASM1_PARS)[19];
  X_S2TSS = mxGetPr(ASM1_PARS)[20];
  X_BH2TSS = mxGetPr(ASM1_PARS)[21];
  X_BA2TSS = mxGetPr(ASM1_PARS)[22];
  X_P2TSS = mxGetPr(ASM1_PARS)[23];
  A = mxGetPr(SPAR)[0];
  C = mxGetPr(SPAR)[1];
  Hmin = mxGetPr(SPAR)[2];
  
  
  if (x[18] < Hmin)
    x[18] = Hmin;
  
  vol = A*x[18];

  for (i = 0; i < 13; i++) {
      if (x[i] < 0.0)
      y[i] = 0.0;
      else
      y[i] = x[i]/vol;               /* ASM1-states: SI, SS, XI, XS, XBH, XBA, XP, SO, SNO, SNH, SND, XND, SALK */
  
  }

  y[13] = (X_I2TSS*x[2]+X_S2TSS*x[3]+X_BH2TSS*x[4]+X_BA2TSS*x[5]+X_P2TSS*x[6])/vol;   
                                     /* TSS */
  y[14] = C*pow((x[18] - Hmin),1.5); /* Q_out */
  y[15] = u[15];                     /* T */
  y[16] = x[13]/vol;                 /* Dummy state 1 */  
  y[17] = x[14]/vol;                 /* Dummy state 2 */
  y[18] = x[15]/vol;                 /* Dummy state 3 */
  y[19] = x[16]/vol;                 /* Dummy state 4 */
  y[20] = x[17]/vol;                 /* Dummy state 5 */
  y[21] = x[18];                         /* Level in tank */
  

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

double vol, A, C, Hmin, qin, qout;
double mu_H, K_S, K_OH, K_NO, b_H, mu_A, K_NH, K_OA, b_A, ny_g, k_a, k_h, K_X, ny_h;
double k_Des, ny_Dec, ny_Bio, KD_ox,  kDec_ox, kBio_ox,  KD_anox, kDec_anox, kBio_anox; 
double Y_H, Y_A, f_P, i_XB, i_XP;
double X_I2TSS, X_S2TSS, X_BH2TSS, X_BA2TSS, X_P2TSS;
double proc1, proc2, proc3, proc4, proc5, proc6, proc7, proc8, proc3x;
double proc9, proc10, proc11, proc12, proc13, proc14, proc15, proc16, proc17;
double reac1, reac2, reac3, reac4, reac5, reac6, reac7, reac8, reac9, reac10, reac11, reac12, reac13, reac16, reac17, reac18, reac19, reac20;
double SO_sat, SO_sat_temp, KLa_temp;
double xtemp[18];
int i;

A = mxGetPr(SPAR)[0];
C = mxGetPr(SPAR)[1];
Hmin = mxGetPr(SPAR)[2];

if (x[18] < Hmin)
    x[18] = Hmin; 

qin = u[14];
qout = C*pow((x[18] - Hmin),1.5);
vol = A*x[18];

///////////////////////////////

mu_H = mxGetPr(ASM1_PARS)[0];
K_S = mxGetPr(ASM1_PARS)[1];
K_OH = mxGetPr(ASM1_PARS)[2];
K_NO = mxGetPr(ASM1_PARS)[3];
b_H = mxGetPr(ASM1_PARS)[4];
mu_A = mxGetPr(ASM1_PARS)[5];
K_NH = mxGetPr(ASM1_PARS)[6];
K_OA = mxGetPr(ASM1_PARS)[7];
b_A = mxGetPr(ASM1_PARS)[8];
ny_g = mxGetPr(ASM1_PARS)[9];
k_a = mxGetPr(ASM1_PARS)[10];
k_h = mxGetPr(ASM1_PARS)[11];
K_X = mxGetPr(ASM1_PARS)[12];
ny_h = mxGetPr(ASM1_PARS)[13];
Y_H = mxGetPr(ASM1_PARS)[14];
Y_A = mxGetPr(ASM1_PARS)[15];
f_P = mxGetPr(ASM1_PARS)[16];
i_XB = mxGetPr(ASM1_PARS)[17];
i_XP = mxGetPr(ASM1_PARS)[18];

X_I2TSS = mxGetPr(ASM1_PARS)[19];
X_S2TSS = mxGetPr(ASM1_PARS)[20];
X_BH2TSS = mxGetPr(ASM1_PARS)[21];
X_BA2TSS = mxGetPr(ASM1_PARS)[22];
X_P2TSS = mxGetPr(ASM1_PARS)[23];

SO_sat_temp = mxGetPr(SOSAT)[0];

KLa_temp = mxGetPr(KLa)[0];

k_Des = mxGetPr(MPOL_PAR)[0];
ny_Dec = mxGetPr(MPOL_PAR)[1];
ny_Bio = mxGetPr(MPOL_PAR)[2];
KD_ox = mxGetPr(MPOL_PAR)[3];
kDec_ox = mxGetPr(MPOL_PAR)[4];
kBio_ox = mxGetPr(MPOL_PAR)[5];
KD_anox = mxGetPr(MPOL_PAR)[6];
kDec_anox = mxGetPr(MPOL_PAR)[7];
kBio_anox = mxGetPr(MPOL_PAR)[8];

for (i = 0; i < 18; i++) {
   if (x[i] < 0.0)
     xtemp[i] == 0.0;
   else
     xtemp[i] = x[i];
}
for (i = 0; i < 18; i++) {
   if (x[i] < 0.0)
     x[i] == 0.0;
   else
     x[i] = x[i];
}

  if (KLa_temp < 0.0)
        x[7] = fabs(KLa_temp);



// // proc1 = mu_H*(xtemp[1]/(K_S+xtemp[1]))*(xtemp[7]/(K_OH+xtemp[7]))*xtemp[4];
// // proc2 = mu_H*(xtemp[1]/(K_S+xtemp[1]))*(K_OH/(K_OH+xtemp[7]))*(xtemp[8]/(K_NO+xtemp[8]))*ny_g*xtemp[4];
// // proc3 = mu_A*(xtemp[9]/(K_NH+xtemp[9]))*(xtemp[7]/(K_OA+xtemp[7]))*xtemp[5];
// // proc4 = b_H*xtemp[4];
// // proc5 = b_A*xtemp[5];
// // proc6 = k_a*xtemp[10]*xtemp[4];
// // proc7 = k_h*((xtemp[3]/xtemp[4])/(K_X+(xtemp[3]/xtemp[4])))*((xtemp[7]/(K_OH+xtemp[7]))+ny_h*(K_OH/(K_OH+xtemp[7]))*(xtemp[8]/(K_NO+xtemp[8])))*xtemp[4];
// // proc8 = proc7*xtemp[11]/xtemp[3];
// 
proc9 =  k_Des*xtemp[16];
proc10 = k_Des* KD_ox *xtemp[13]*(xtemp[7]/(K_OH+xtemp[7]))*(X_I2TSS*xtemp[2]+X_S2TSS*xtemp[3]+X_BH2TSS*xtemp[4]+X_BA2TSS*xtemp[5]+X_P2TSS*xtemp[6]);
proc11 = kDec_ox * xtemp[14]*(K_S*ny_Dec /((K_S*ny_Dec)+xtemp[1]))*(xtemp[7]/(K_OH+xtemp[7]))*(X_I2TSS*xtemp[2]+X_S2TSS*xtemp[3]+X_BH2TSS*xtemp[4]+X_BA2TSS*xtemp[5]+X_P2TSS*xtemp[6]);
proc12 = kBio_ox * xtemp[13]*(K_S*ny_Bio /((K_S*ny_Bio)+xtemp[1]))*(xtemp[7]/(K_OH+xtemp[7]))*(X_I2TSS*xtemp[2]+X_S2TSS*xtemp[3]+X_BH2TSS*xtemp[4]+X_BA2TSS*xtemp[5]+X_P2TSS*xtemp[6]);
proc13 = k_Des* KD_anox *  xtemp[13] * (K_OH/(K_OH+xtemp[7]))*(X_I2TSS*xtemp[2]+X_S2TSS*xtemp[3]+X_BH2TSS*xtemp[4]+X_BA2TSS*xtemp[5]+X_P2TSS*xtemp[6]);
proc14 = kDec_anox* xtemp[14]*(K_S*ny_Dec /((K_S*ny_Dec)+xtemp[1]))* (K_OH/(K_OH+xtemp[7]))*(X_I2TSS*xtemp[2]+X_S2TSS*xtemp[3]+X_BH2TSS*xtemp[4]+X_BA2TSS*xtemp[5]+X_P2TSS*xtemp[6]);
proc15 = kBio_anox* xtemp[13]*(K_S*ny_Bio /((K_S*ny_Bio)+xtemp[1]))* (K_OH/(K_OH+xtemp[7]))*(X_I2TSS*xtemp[2]+X_S2TSS*xtemp[3]+X_BH2TSS*xtemp[4]+X_BA2TSS*xtemp[5]+X_P2TSS*xtemp[6]);

// 
// // reac1 = 0.0;
// // reac2 = (-proc1-proc2)/Y_H+proc7;
// // reac3 = 0.0;
// // reac4 = (1.0-f_P)*(proc4+proc5)-proc7;
// // reac5 = proc1+proc2-proc4;
// // reac6 = proc3-proc5;
// // reac7 = f_P*(proc4+proc5);
// // reac8 = -((1.0-Y_H)/Y_H)*proc1-((4.57-Y_A)/Y_A)*proc3;
// // reac9 = -((1.0-Y_H)/(2.86*Y_H))*proc2+proc3/Y_A;
// // reac10 = -i_XB*(proc1+proc2)-(i_XB+(1.0/Y_A))*proc3+proc6;
// // reac11 = -proc6+proc8;
// // reac12 = (i_XB-f_P*i_XP)*(proc4+proc5)-proc8;
// // reac13 = -i_XB/14.0*proc1+((1.0-Y_H)/(14.0*2.86*Y_H)-(i_XB/14.0))*proc2-((i_XB/14.0)+1.0/(7.0*Y_A))*proc3+proc6/14.0;
// 
// 
// reac16 = proc9 - proc10 + proc11 -proc12 -proc13 +proc14 - proc15;
// reac17 = -proc11 -proc14;
// reac18 = 0.0;
// reac19 = -proc9 + proc10 + proc13;
// reac20 = 0.0;

reac1 = 0.0;
reac2 = 0.0;
reac3 = 0.0;
reac4 = 0.0;
reac5 = 0.0;
reac6 = 0.0;
reac7 = 0.0;
reac8 = 0.0;
reac9 = 0.0;
reac10 = 0.0;
reac11 = 0.0;
reac12 = 0.0;
reac13 = 0.0;


reac16 = 0.0;
reac17 = 0.0;
reac18 = 0.0;
reac19 = 0.0;
reac20 = 0.0;

/* NOTE: the states are expressed in mass (not concentration) */
dx[0] = (qin*u[0]-qout*x[0]/vol) + reac1;    /* ASM1-state: SI */
dx[1] = (qin*u[1]-qout*x[1]/vol) + reac2;    /* ASM1-state: SS */
dx[2] = (qin*u[2]-qout*x[2]/vol) + reac3;    /* ASM1-state: XI */
dx[3] = (qin*u[3]-qout*x[3]/vol) + reac4;    /* ASM1-state: XS */
dx[4] = (qin*u[4]-qout*x[4]/vol) + reac5;    /* ASM1-state: XBH */
dx[5] = (qin*u[5]-qout*x[5]/vol) + reac6;    /* ASM1-state: XBA */
dx[6] = (qin*u[6]-qout*x[6]/vol) + reac7;    /* ASM1-state: XP */

if (KLa_temp < 0.0)
    dx[7] = 0.0;
else
    dx[7] = (qin*u[7]-qout*x[7]/vol) + reac8 + KLa_temp*(SO_sat_temp - x[7]);    /* ASM1-state: SO */
    
    
dx[7] = (qin*u[7]-qout*x[7]/vol) + reac8;
dx[8] = (qin*u[8]-qout*x[8]/vol) + reac9;    /* ASM1-state: SNO */
dx[9] = (qin*u[9]-qout*x[9]/vol) + reac10;    /* ASM1-state: SNH */
dx[10] = (qin*u[10]-qout*x[10]/vol) + reac11; /* ASM1-state: SND */
dx[11] = (qin*u[11]-qout*x[11]/vol) + reac12; /* ASM1-state: XND */
dx[12] = (qin*u[12]-qout*x[12]/vol) + reac13; /* ASM1-state: SALK */

dx[13] = (qin*u[16]-qout*x[13]/vol) + reac16; /* Cli */
dx[14] = (qin*u[17]-qout*x[14]/vol) + reac17; /* Ccj */
dx[15] = (qin*u[18]-qout*x[15]/vol) + reac18; /* S_D3 */
dx[16] = (qin*u[19]-qout*x[16]/vol) + reac19; /* Csl */
dx[17] = (qin*u[20]-qout*x[17]/vol) + reac20; /* X_D2 */

dx[18] = 1/A*(qin - qout);           /* h, level in tank */

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
