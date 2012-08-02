/*
 * MATLAB Compiler: 2.0
 * Date: Thu Jun 21 14:04:38 2001
 * Arguments: "-t" "G2K" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __care_h
#define __care_h 1

#include "matlab.h"

extern mxArray * mlfNCare(int nargout,
                          mxArray * * L,
                          mxArray * * G,
                          mxArray * * RR,
                          mxArray * A,
                          mxArray * B,
                          mxArray * Q,
                          ...);
extern mxArray * mlfCare(mxArray * * L,
                         mxArray * * G,
                         mxArray * * RR,
                         mxArray * A,
                         mxArray * B,
                         mxArray * Q,
                         ...);
extern void mlfVCare(mxArray * A, mxArray * B, mxArray * Q, ...);
extern void mlxCare(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#endif
