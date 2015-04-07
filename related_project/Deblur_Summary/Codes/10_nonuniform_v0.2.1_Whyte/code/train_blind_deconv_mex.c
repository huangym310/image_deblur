#include <mex.h>
#include <stddef.h>
#include <math.h>
#include "matrix.h"
#include "ow_homography.h"

/*
    Function to do the heavy computation for variational blind deblurring
        Single-threaded version

    Author:     Oliver Whyte <oliver.whyte@ens.fr>
    Date:       August 2010
    Copyright:  2010, Oliver Whyte
    Reference:  O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
    URL:        http://www.di.ens.fr/~whyte/deblurring/
*/

double sqr(double const x) { return x*x; }

/* Global variables */
double *Kblurry,*Ksharp;
double invKblurry[9];
int h_sharp, w_sharp, h_blurry, w_blurry;
double *theta_list;
double *Hcache;

int sharp_subim_l[2];
int blurry_subim_l[2];

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray const* mf1_ptr       = prhs[0]; /* first moment of q(f) */
    mxArray const* mf2_ptr       = prhs[1]; /* raw second moment of q(f) */
    mxArray const* mw1_ptr       = prhs[2]; /* first moment of q(w) */
    mxArray const* mw2_ptr       = prhs[3]; /* raw second moment of q(w) */
    mxArray const* g_ptr         = prhs[4]; /* blurry image */
    mxArray const* Ksharp_ptr    = prhs[5]; /* intrinsic calibration of sharp image */
    mxArray const* Kblurry_ptr   = prhs[6]; /* Intrinsic calibration of blurry image */
    mxArray const* theta_ptr     = prhs[7]; /* angles covered by camera kernel */
    mxArray const* obsmask_ptr   = prhs[8]; /* Mask of observed blurry pixels */
    double const mmu1            = mxGetScalar(prhs[9]); 
    double const mmu2            = mxGetScalar(prhs[10]);

    int const n_kernel = mxGetN(theta_ptr);

    double const *mf1   = mxGetPr(mf1_ptr);
    double const *mf2   = mxGetPr(mf2_ptr);
    double const *mw1   = mxGetPr(mw1_ptr);
    double const *mw2   = mxGetPr(mw2_ptr);
    double const *g = mxGetPr(g_ptr);
    double const *obsmask = mxGetPr(obsmask_ptr);

    int const *dims_sharp = mxGetDimensions(mf1_ptr);
    int const *dims_blurry = mxGetDimensions(g_ptr);
    int const n_sharp = mxGetNumberOfElements(mf1_ptr);
    int const n_blurry = mxGetNumberOfElements(g_ptr);

    mxArray* Ai_ptr;
    double* Ai ;
    double gbari;
    mxArray* exp_sqr_err_ptr ;
    double* exp_sqr_err;
    mxArray* f1f2_ptr;
    double* f1f2;
    mxArray* f2_ptr;
    double* f2;
    mxArray* w1w2_ptr;
    double* w1w2;
    mxArray* w2_ptr;
    double* w2;
    mxArray* mu1_ptr;
    double* mu1;
    mxArray* mu2_ptr;
    double* mu2;
    mxArray* rerror_ptr;
    double* rerror;
    mxArray* Hcache_ptr;
    int j, k;
    double erri;
    double* coeffs = (double*)mxMalloc(sizeof(double)*n_interp);
    mxArray *var_w_ptr;
    double *var_w;
    /* Calculate variance of each sharp pixel */
    mxArray *var_f_ptr;
    double *var_f;
    int xi, yi, subimage;
    int xinterp, yinterp;

    double term2i;
    double term3i;
    double term4i;
    double xj,yj;
    int xjfloor,yjfloor;
    double Bik = 0;

    int jj;
    bool inbounds;

    if(mxGetNumberOfDimensions(mf1_ptr) > 2) mexErrMsgTxt("Images must be grayscale");
    if(mxGetNumberOfDimensions(g_ptr) > 2) mexErrMsgTxt("Images must be grayscale");

    Ksharp = mxGetPr(Ksharp_ptr);
    Kblurry = mxGetPr(Kblurry_ptr);
    theta_list = mxGetPr(theta_ptr);

    h_sharp = dims_sharp[0];
    w_sharp = dims_sharp[1]/2;
    h_blurry = dims_blurry[0];
    w_blurry = dims_blurry[1]/2;
    
    /* Offsets for x and y gradients */
    sharp_subim_l[0] = 1;
    sharp_subim_l[1] = 1+w_sharp;
    blurry_subim_l[0] = 1;
    blurry_subim_l[1] = 1+w_blurry;

    /* Invert calibration matrix for blurry image */
    inv3(Kblurry,invKblurry);

    /* general */
    Ai_ptr = mxCreateNumericArray(2, dims_sharp, mxDOUBLE_CLASS, mxREAL);
    Ai      = mxGetPr(Ai_ptr);
    /* exp_sqr_err */
    exp_sqr_err_ptr = mxCreateNumericArray(2, dims_blurry, mxDOUBLE_CLASS, mxREAL);
    exp_sqr_err      = mxGetPr(exp_sqr_err_ptr);
    /* f1f2 */
    f1f2_ptr = mxCreateNumericArray(2, dims_sharp, mxDOUBLE_CLASS, mxREAL);
    f1f2      = mxGetPr(f1f2_ptr);
    /* f2 */
    f2_ptr = mxCreateNumericArray(2, dims_sharp, mxDOUBLE_CLASS, mxREAL);
    f2      = mxGetPr(f2_ptr);
    /* w1w2 */
    w1w2_ptr = mxCreateDoubleMatrix(n_kernel, 1, mxREAL);
    w1w2      = mxGetPr(w1w2_ptr);
    /* w2 */
    w2_ptr = mxCreateDoubleMatrix(n_kernel, 1, mxREAL);
    w2      = mxGetPr(w2_ptr);
    /* mu1=sum(sum(Dp.*(D-mD))); */
    mu1_ptr = mxCreateDoubleMatrix(1, 1, mxREAL);
    mu1      = mxGetPr(mu1_ptr);
    /* mu2 = data_points */
    mu2_ptr = mxCreateDoubleMatrix(1, 1, mxREAL);
    mu2      = mxGetPr(mu2_ptr);
    /* rerror = sum(exp_sqr_err) */
    rerror_ptr = mxCreateDoubleMatrix(1, 1, mxREAL);
    rerror      = mxGetPr(rerror_ptr);
    /* Hcache, used to cache homography matrices */
    Hcache_ptr = mxCreateDoubleMatrix(n_kernel, 9, mxREAL);
    Hcache              = mxGetPr(Hcache_ptr);

    /* Calculate variance of each kernel element */
    var_w_ptr = mxCreateDoubleMatrix(1,n_kernel,mxREAL);
    var_w = mxGetPr(var_w_ptr);
    for(k=0; k<n_kernel; ++k)
        var_w[k] = mw2[k] - sqr(mw1[k]);
    /* Calculate variance of each sharp pixel */
    var_f_ptr = mxCreateNumericArray(2, dims_sharp, mxDOUBLE_CLASS, mxREAL);
    var_f      = mxGetPr(var_f_ptr);
    for(j=0; j<n_sharp; ++j)
        var_f[j] = mf2[j] - sqr(mf1[j]);
    /* Calculate homographies and cache them */
    for(k=0; k<n_kernel; ++k)
        compute_homography_matrix(Ksharp, &theta_list[k*3], invKblurry, &Hcache[k*9]);

    /* Loop over all blurry image pixels:
    (xi,yi): 1-based row/column indices into blurry image */

    /* For each subimage, (x-gradients then y-gradients) */
    for(subimage=0; subimage<=1; ++subimage) {
        /* Loop over all blurry pixels i */
        for(xi=1; xi<=w_blurry; ++xi) {
            for(yi=1; yi<=h_blurry; ++yi) {
                int i = (xi-1+blurry_subim_l[subimage]-1)*h_blurry + yi - 1;
                if(obsmask[i])
                    ++(*mu2); /* Increment observation count */
                else
                    continue; /* Skip this blurry pixel */
                term2i = 0;
                term3i = 0;
                term4i = 0;
                /* Clear temporary variables */
                gbari = 0;
                /* Clear Ai */
                for(j=0; j<n_sharp; ++j)
                    Ai[j] = 0;
                /* Loop over all rotations k */
                for(k=0; k<n_kernel; ++k) {
                    Bik = 0;
                    /* Project blurry pixel into sharp image */
                    project_blurry_to_sharp(yi,xi,&Hcache[k*9],&xj,&yj);
                    /* Get interpolation coefficients */
                    interp_coeffs(xj,yj,&xjfloor,&yjfloor,coeffs);
                    /* Interpolate points */
                    for(jj=0; jj<n_interp; ++jj) {
                        yinterp = yjfloor + yoff[jj];
                        xinterp = xjfloor + xoff[jj];
                        inbounds = yinterp >=1 && yinterp <= h_sharp && xinterp >= 1 && xinterp <= w_sharp;
                        if(!inbounds)
                            continue;
                        j = (xinterp-1+sharp_subim_l[subimage]-1)*h_sharp + yinterp - 1;
                        Bik += coeffs[jj]*mf1[j];
                        term2i += sqr(coeffs[jj])*var_f[j]*var_w[k];
                        Ai[j] += coeffs[jj]*mw1[k];
                        f2[j] += sqr(coeffs[jj])*var_w[k];
                        w2[k] += sqr(coeffs[jj])*var_f[j];
                    }
                    /* w2 */
                    w2[k] += sqr(Bik);
                    /* Accumulate gbar[i] */
                    gbari += Bik*mw1[k];
                    /* exp_sqr_err */
                    term4i += sqr(Bik)*var_w[k];
                    /* Having calculated Bik, interpolate terms involving it */
                    for(jj=0; jj<n_interp; ++jj) {
                        yinterp = yjfloor + yoff[jj];
                        xinterp = xjfloor + xoff[jj];
                        inbounds = yinterp >=1 && yinterp <= h_sharp && xinterp >= 1 && xinterp <= w_sharp;
                        if(!inbounds)
                            continue;
                        j = (xinterp-1+sharp_subim_l[subimage]-1)*h_sharp + yinterp - 1;
                        f1f2[j] += coeffs[jj]*Bik*var_w[k];
                    }
                }
                /* Expected error */
                erri = obsmask[i]*(g[i] - gbari - mmu1);
                *mu1 += obsmask[i]*(g[i] - gbari);
                /* Accumulate terms for image parameters */
                for(j=0; j<n_sharp; ++j) {
                    if(Ai[j] > 0) {
                        term3i += sqr(Ai[j])*var_f[j];
                        f2[j] += sqr(Ai[j]);
                        f1f2[j] += Ai[j]*erri;
                    }
                }
                /* Loop over rotations again to calculate terms for kernel */
                for(k=0; k<n_kernel; ++k) {
                    Bik = 0;
                    /* Project blurry pixel into sharp image */
                    project_blurry_to_sharp(yi,xi,&Hcache[k*9],&xj,&yj);
                    /* Get interpolation coefficients */
                    interp_coeffs(xj,yj,&xjfloor,&yjfloor,coeffs);
                    /* Interpolate points */
                    for(jj=0; jj<n_interp; ++jj) {
                        yinterp = yjfloor + yoff[jj];
                        xinterp = xjfloor + xoff[jj];
                        inbounds = yinterp >=1 && yinterp <= h_sharp && xinterp >= 1 && xinterp <= w_sharp;
                        if(!inbounds)
                            continue;
                        j = (xinterp-1+sharp_subim_l[subimage]-1)*h_sharp + yinterp - 1;
                        Bik += coeffs[jj]*mf1[j];
                        w1w2[k] -= coeffs[jj]*Ai[j]*var_f[j];
                    }
                    /* w1w2 */
                    w1w2[k] += Bik*erri;
                }
                /* Add everything together for exp_sqr_err */
                exp_sqr_err[i] = sqr(erri) + term2i + term3i + term4i + mmu2 - sqr(mmu1);
                /* rerror is the sum of exp_sqr_err over all pixels */
                *rerror += exp_sqr_err[i];
            }
        }
    }
    /* Add parts together */
    for(j=0; j<n_sharp; ++j) {
        f1f2[j] += mf1[j]*f2[j];
    }
    for(k=0; k<n_kernel; ++k) {
        w1w2[k] += mw1[k]*w2[k];
    }
    plhs[0] = f1f2_ptr;
    plhs[1] = f2_ptr;
    plhs[2] = w1w2_ptr;
    plhs[3] = w2_ptr;
    plhs[4] = mu1_ptr;
    plhs[5] = mu2_ptr;
    plhs[6] = rerror_ptr;

    mxFree(coeffs);
}
  



