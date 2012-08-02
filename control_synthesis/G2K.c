/*
 * MATLAB Compiler: 2.0
 * Date: Thu Jun 21 14:04:38 2001
 * Arguments: "-t" "G2K" 
 */
#include "G2K.h"
#include "care.h"

/*
 * The function "MG2K" is the implementation version of the "G2K" M-function
 * from file "D:\matlab\MPC\G2K.m" (lines 1-64). It contains the actual
 * compiled code for that M-function. It is a static function and must only be
 * called from one of the interface functions, appearing below.
 */
/*
 * function [K, gamma_opt] = G2K(G);
 */
static mxArray * MG2K(mxArray * * gamma_opt, int nargout_, mxArray * G) {
    mxArray * K = mclGetUninitializedArray();
    mxArray * A = mclGetUninitializedArray();
    mxArray * Ac = mclGetUninitializedArray();
    mxArray * B = mclGetUninitializedArray();
    mxArray * Bc = mclGetUninitializedArray();
    mxArray * C = mclGetUninitializedArray();
    mxArray * Cc = mclGetUninitializedArray();
    mxArray * D = mclGetUninitializedArray();
    mxArray * Dc = mclGetUninitializedArray();
    mxArray * E = mclGetUninitializedArray();
    mxArray * EE = mclGetUninitializedArray();
    mxArray * F = mclGetUninitializedArray();
    mxArray * R = mclGetUninitializedArray();
    mxArray * S = mclGetUninitializedArray();
    mxArray * U = mclGetUninitializedArray();
    mxArray * V = mclGetUninitializedArray();
    mxArray * W = mclGetUninitializedArray();
    mxArray * X = mclGetUninitializedArray();
    mxArray * Y = mclGetUninitializedArray();
    mxArray * a = mclGetUninitializedArray();
    mxArray * a11 = mclGetUninitializedArray();
    mxArray * a12 = mclGetUninitializedArray();
    mxArray * a21 = mclGetUninitializedArray();
    mxArray * a22 = mclGetUninitializedArray();
    mxArray * b = mclGetUninitializedArray();
    mxArray * b1 = mclGetUninitializedArray();
    mxArray * b2 = mclGetUninitializedArray();
    mxArray * c = mclGetUninitializedArray();
    mxArray * c1 = mclGetUninitializedArray();
    mxArray * c2 = mclGetUninitializedArray();
    mxArray * d = mclGetUninitializedArray();
    mxArray * n = mclGetUninitializedArray();
    mxArray * o = mclGetUninitializedArray();
    mxArray * ss = mclUnassigned();
    mclValidateInputs("G2K", 1, &G);
    /*
     * 
     * A = G.A; B = G.b; C = G.C; D = G.D;
     */
    mlfAssign(&A, mlfIndexRef(G, ".A"));
    mlfAssign(&B, mlfIndexRef(G, ".b"));
    mlfAssign(&C, mlfIndexRef(G, ".C"));
    mlfAssign(&D, mlfIndexRef(G, ".D"));
    /*
     * 
     * % Calculate optimal controller
     * 
     * S = eye(size(D,2)) + D'*D;
     */
    mlfAssign(
      &S,
      mlfPlus(
        mlfEye(mlfSize(mclValueVarargout(), D, mlfScalar(2.0)), NULL),
        mlfMtimes(mlfCtranspose(D), D)));
    /*
     * R = eye(size(D,1)) + D*D';
     */
    mlfAssign(
      &R,
      mlfPlus(
        mlfEye(mlfSize(mclValueVarargout(), D, mlfScalar(1.0)), NULL),
        mlfMtimes(D, mlfCtranspose(D))));
    /*
     * 
     * [Y] = care((A-B/S*D'*C)', C', B/S*B', R);
     */
    mlfAssign(
      &Y,
      mlfNCare(
        1,
        NULL,
        NULL,
        NULL,
        mlfCtranspose(
          mlfMinus(
            A, mlfMtimes(mlfMtimes(mlfMrdivide(B, S), mlfCtranspose(D)), C))),
        mlfCtranspose(C),
        mlfMtimes(mlfMrdivide(B, S), mlfCtranspose(B)),
        R,
        NULL));
    /*
     * 
     * [X] = care((A-B/S*D'*C), B, C'/R*C, S);
     */
    mlfAssign(
      &X,
      mlfNCare(
        1,
        NULL,
        NULL,
        NULL,
        mlfMinus(
          A, mlfMtimes(mlfMtimes(mlfMrdivide(B, S), mlfCtranspose(D)), C)),
        B,
        mlfMtimes(mlfMrdivide(mlfCtranspose(C), R), C),
        S,
        NULL));
    /*
     * 
     * gamma_opt = sqrt(1 + max(eig(X*Y)));
     */
    mlfAssign(
      gamma_opt,
      mlfSqrt(
        mlfPlus(
          mlfScalar(1.0),
          mlfMax(NULL, mlfEig(NULL, mlfMtimes(X, Y), NULL), NULL, NULL))));
    /*
     * 
     * F = S\(D'*C + B'*X);
     */
    mlfAssign(
      &F,
      mlfMldivide(
        S,
        mlfPlus(
          mlfMtimes(mlfCtranspose(D), C), mlfMtimes(mlfCtranspose(B), X))));
    /*
     * W = (gamma_opt*gamma_opt - 1)*eye(size(X,2)) - Y*X;
     */
    mlfAssign(
      &W,
      mlfMinus(
        mlfMtimes(
          mlfMinus(mlfMtimes(*gamma_opt, *gamma_opt), mlfScalar(1.0)),
          mlfEye(mlfSize(mclValueVarargout(), X, mlfScalar(2.0)), NULL)),
        mlfMtimes(Y, X)));
    /*
     * 
     * Ac = W*(A-B*F) - gamma_opt*gamma_opt*(Y*C')*(C-D*F);
     */
    mlfAssign(
      &Ac,
      mlfMinus(
        mlfMtimes(W, mlfMinus(A, mlfMtimes(B, F))),
        mlfMtimes(
          mlfMtimes(
            mlfMtimes(*gamma_opt, *gamma_opt), mlfMtimes(Y, mlfCtranspose(C))),
          mlfMinus(C, mlfMtimes(D, F)))));
    /*
     * Bc = gamma_opt*gamma_opt*Y*C';
     */
    mlfAssign(
      &Bc,
      mlfMtimes(
        mlfMtimes(mlfMtimes(*gamma_opt, *gamma_opt), Y), mlfCtranspose(C)));
    /*
     * Cc = -B'*X;
     */
    mlfAssign(&Cc, mlfMtimes(mlfUminus(mlfCtranspose(B)), X));
    /*
     * Dc = -D';
     */
    mlfAssign(&Dc, mlfUminus(mlfCtranspose(D)));
    /*
     * 
     * % Get minimal realisation of controller
     * [U,EE,V] = svd(W);
     */
    mlfAssign(&U, mlfSvd(&EE, &V, W, NULL));
    /*
     * o = rank(EE);
     */
    mlfAssign(&o, mlfRank(EE, NULL));
    /*
     * n = size(EE,1) - o;
     */
    mlfAssign(
      &n, mlfMinus(mlfSize(mclValueVarargout(), EE, mlfScalar(1.0)), o));
    /*
     * E = EE(1:o,1:o);
     */
    mlfAssign(
      &E,
      mlfIndexRef(
        EE,
        "(?,?)",
        mlfColon(mlfScalar(1.0), o, NULL),
        mlfColon(mlfScalar(1.0), o, NULL)));
    /*
     * 
     * Ac = U'*Ac*V;
     */
    mlfAssign(&Ac, mlfMtimes(mlfMtimes(mlfCtranspose(U), Ac), V));
    /*
     * Bc = U'*Bc;
     */
    mlfAssign(&Bc, mlfMtimes(mlfCtranspose(U), Bc));
    /*
     * Cc = Cc*V;
     */
    mlfAssign(&Cc, mlfMtimes(Cc, V));
    /*
     * Dc = Dc;
     */
    mlfAssign(&Dc, Dc);
    /*
     * 
     * % partition
     * a11 = Ac(1:o,1:o);
     */
    mlfAssign(
      &a11,
      mlfIndexRef(
        Ac,
        "(?,?)",
        mlfColon(mlfScalar(1.0), o, NULL),
        mlfColon(mlfScalar(1.0), o, NULL)));
    /*
     * a12 = Ac(1:o,o+1:end);
     */
    mlfAssign(
      &a12,
      mlfIndexRef(
        Ac,
        "(?,?)",
        mlfColon(mlfScalar(1.0), o, NULL),
        mlfColon(
          mlfPlus(o, mlfScalar(1.0)),
          mlfEnd(Ac, mlfScalar(2), mlfScalar(2)),
          NULL)));
    /*
     * a21 = Ac(1+o:end,1:o);
     */
    mlfAssign(
      &a21,
      mlfIndexRef(
        Ac,
        "(?,?)",
        mlfColon(
          mlfPlus(mlfScalar(1.0), o),
          mlfEnd(Ac, mlfScalar(1), mlfScalar(2)),
          NULL),
        mlfColon(mlfScalar(1.0), o, NULL)));
    /*
     * a22 = Ac(1+o:end,1+o:end);
     */
    mlfAssign(
      &a22,
      mlfIndexRef(
        Ac,
        "(?,?)",
        mlfColon(
          mlfPlus(mlfScalar(1.0), o),
          mlfEnd(Ac, mlfScalar(1), mlfScalar(2)),
          NULL),
        mlfColon(
          mlfPlus(mlfScalar(1.0), o),
          mlfEnd(Ac, mlfScalar(2), mlfScalar(2)),
          NULL)));
    /*
     * b1 = Bc(1:o,:);
     */
    mlfAssign(
      &b1,
      mlfIndexRef(
        Bc, "(?,?)", mlfColon(mlfScalar(1.0), o, NULL), mlfCreateColonIndex()));
    /*
     * b2 = Bc(1+o:end,:);
     */
    mlfAssign(
      &b2,
      mlfIndexRef(
        Bc,
        "(?,?)",
        mlfColon(
          mlfPlus(mlfScalar(1.0), o),
          mlfEnd(Bc, mlfScalar(1), mlfScalar(2)),
          NULL),
        mlfCreateColonIndex()));
    /*
     * c1 = Cc(:,1:o);
     */
    mlfAssign(
      &c1,
      mlfIndexRef(
        Cc, "(?,?)", mlfCreateColonIndex(), mlfColon(mlfScalar(1.0), o, NULL)));
    /*
     * c2 = Cc(:,1+o:end);
     */
    mlfAssign(
      &c2,
      mlfIndexRef(
        Cc,
        "(?,?)",
        mlfCreateColonIndex(),
        mlfColon(
          mlfPlus(mlfScalar(1.0), o),
          mlfEnd(Cc, mlfScalar(2), mlfScalar(2)),
          NULL)));
    /*
     * d = Dc;
     */
    mlfAssign(&d, Dc);
    /*
     * 
     * % introduce transformation & put minimal realisation
     * a = E\(a11-a12/a22*a21);
     */
    mlfAssign(
      &a, mlfMldivide(E, mlfMinus(a11, mlfMtimes(mlfMrdivide(a12, a22), a21))));
    /*
     * b = E\(b1-a12/a22*b2);
     */
    mlfAssign(
      &b, mlfMldivide(E, mlfMinus(b1, mlfMtimes(mlfMrdivide(a12, a22), b2))));
    /*
     * c = c1-c2/a22*a21;
     */
    mlfAssign(&c, mlfMinus(c1, mlfMtimes(mlfMrdivide(c2, a22), a21)));
    /*
     * d = d-c2/a22*b2;
     */
    mlfAssign(&d, mlfMinus(d, mlfMtimes(mlfMrdivide(c2, a22), b2)));
    /*
     * 
     * K = ss(a,b,c,d);
     */
    mlfAssign(&K, mlfIndexRef(ss, "(?,?,?,?)", a, b, c, d));
    mclValidateOutputs("G2K", 2, nargout_, &K, gamma_opt);
    mxDestroyArray(A);
    mxDestroyArray(Ac);
    mxDestroyArray(B);
    mxDestroyArray(Bc);
    mxDestroyArray(C);
    mxDestroyArray(Cc);
    mxDestroyArray(D);
    mxDestroyArray(Dc);
    mxDestroyArray(E);
    mxDestroyArray(EE);
    mxDestroyArray(F);
    mxDestroyArray(R);
    mxDestroyArray(S);
    mxDestroyArray(U);
    mxDestroyArray(V);
    mxDestroyArray(W);
    mxDestroyArray(X);
    mxDestroyArray(Y);
    mxDestroyArray(a);
    mxDestroyArray(a11);
    mxDestroyArray(a12);
    mxDestroyArray(a21);
    mxDestroyArray(a22);
    mxDestroyArray(b);
    mxDestroyArray(b1);
    mxDestroyArray(b2);
    mxDestroyArray(c);
    mxDestroyArray(c1);
    mxDestroyArray(c2);
    mxDestroyArray(d);
    mxDestroyArray(n);
    mxDestroyArray(o);
    mxDestroyArray(ss);
    /*
     * 
     * 
     * 
     * 
     * 
     * 
     * 
     * 
     * 
     * 
     * 
     */
    return K;
}

/*
 * The function "mlfG2K" contains the normal interface for the "G2K" M-function
 * from file "D:\matlab\MPC\G2K.m" (lines 1-64). This function processes any
 * input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
mxArray * mlfG2K(mxArray * * gamma_opt, mxArray * G) {
    int nargout = 1;
    mxArray * K = mclGetUninitializedArray();
    mxArray * gamma_opt__ = mclGetUninitializedArray();
    mlfEnterNewContext(1, 1, gamma_opt, G);
    if (gamma_opt != NULL) {
        ++nargout;
    }
    K = MG2K(&gamma_opt__, nargout, G);
    mlfRestorePreviousContext(1, 1, gamma_opt, G);
    if (gamma_opt != NULL) {
        mclCopyOutputArg(gamma_opt, gamma_opt__);
    } else {
        mxDestroyArray(gamma_opt__);
    }
    return mlfReturnValue(K);
}

/*
 * The function "mlxG2K" contains the feval interface for the "G2K" M-function
 * from file "D:\matlab\MPC\G2K.m" (lines 1-64). The feval function calls the
 * implementation version of G2K through this function. This function processes
 * any input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
void mlxG2K(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[2];
    int i;
    if (nlhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: G2K Line: 1 Column: 0 The function \"G2K\""
            " was called with more than the declared number of outputs (2)"));
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: G2K Line: 1 Column: 0 The function \"G2K\""
            " was called with more than the declared number of inputs (1)"));
    }
    for (i = 0; i < 2; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 1 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 1; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 1, mprhs[0]);
    mplhs[0] = MG2K(&mplhs[1], nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 2 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 2; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}
