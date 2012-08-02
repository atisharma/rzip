/*
 * MATLAB Compiler: 2.0
 * Date: Thu Jun 21 14:04:38 2001
 * Arguments: "-t" "G2K" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __G2K_h
#define __G2K_h 1

#include "matlab.h"

extern mxArray * mlfG2K(mxArray * * gamma_opt, mxArray * G);
extern void mlxG2K(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#endif
