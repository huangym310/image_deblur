#include <mex.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "ow_mat3.h"

/*
    Function to mask the edge of a blurry image, affected by boundary effects
    
    Author:     Oliver Whyte <oliver.whyte@ens.fr>
    Date:       August 2010
    Copyright:  2010, Oliver Whyte
    Reference:  O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
    URL:        http://www.di.ens.fr/~whyte/deblurring/
*/

#define thetax(k) (theta[(k)*3])
#define thetay(k) (theta[(k)*3+1])
#define thetaz(k) (theta[(k)*3+2])

double round(double val)
{
    return floor(val + 0.5);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray const* imorig_ptr = prhs[0]; /* original / sharp image */
    mxArray const* dimswarp_ptr = prhs[1]; /* dimensions of warped / blurry image */
    mxArray const* Ko_ptr = prhs[2]; /* Intrinsic calibration of original / sharp image */
    mxArray const* Kw_ptr = prhs[3]; /* intrinsic calibration of warped / blurry image */
    mxArray const* theta_ptr = prhs[4]; /* angles covered by camera kernel */
    mxArray const* kernel_ptr = prhs[5]; /* kernel */
    bool const clamp_edges_to_zero = (bool const)mxGetScalar(prhs[6]); /* Set columns to zero if they involve edge pixels */
    bool const non_uniform = (bool const)mxGetScalar(prhs[7]); /* Non-uniform blur model? */
    mxArray const* theta_pre_ptr = prhs[8]; /* Registration of noisy image */

    int const nk = mxGetN(theta_ptr);

    double const *imorig = mxGetPr(imorig_ptr);
    double const *dimswarp = mxGetPr(dimswarp_ptr);
    double const *Ko = mxGetPr(Ko_ptr);
    double const *Kw = mxGetPr(Kw_ptr);
    double const *theta = mxGetPr(theta_ptr);
    double const *kernel = mxGetPr(kernel_ptr);
    double const *theta_pre = mxGetPr(theta_pre_ptr);
    /* double const *kernel = mxGetPr(kernel_ptr); */

    int const numdims = mxGetNumberOfDimensions(imorig_ptr);
    int const *imdimso = mxGetDimensions(imorig_ptr);
    int const heiorig = imdimso[0];
    int const widorig = imdimso[1];
    int const channels = (numdims==3) ? imdimso[2] : 1;
    /* int const heiorig = (int const)mxGetM(imorig_ptr); */
    /*   int const widorig = (int const)mxGetN(imorig_ptr); */
    int const heiwarp = (int const)dimswarp[0];
    int const widwarp = (int const)dimswarp[1];
    int const norig = heiorig*widorig;
    int const nwarp = heiwarp*widwarp;

    int imdimsw[3];

    double invKw[9];
    mxArray* imwarp_ptr;
    double *imwarp;

    int ro, co;
    int i, j, k;
    double denom, yo, xo, yoi, xoi;
    double H[9];
    double R[9];
    double f1, f2, f3, f4;

    int cw, rw;

    inv3(Kw,invKw);

    imdimsw[0] = heiwarp;
    imdimsw[1] = widwarp;
    imdimsw[2] = 1;
    
    imwarp_ptr = mxCreateNumericArray(3, imdimsw, mxDOUBLE_CLASS, mxREAL);
    imwarp = mxGetPr(imwarp_ptr);

    /* Check thetas */
    if(!non_uniform) {
        for(k=0; k<nk; ++k)
            if(fabs(round(thetay(k))-thetay(k)) > 1e-5 ||     \
            fabs(round(thetax(k))-thetax(k)) > 1e-5)
            mexErrMsgTxt("theta not integral for uniform blur");
        if(fabs(round(theta_pre[0])-theta_pre[0]) > 1e-5 || \
            fabs(round(theta_pre[1])-theta_pre[1]) > 1e-5)
            mexErrMsgTxt("theta_pre not integral for uniform blur");
    }

    /* Loop over all rotations */   
    if (non_uniform) {
        for(k=0; k<nk; ++k) {
            cp3(invKw,H);
            rot3(thetax(k),thetay(k),thetaz(k),R);
            mmip3(R,H);
            rot3(theta_pre[0],theta_pre[1],theta_pre[2],R);
            mmip3(R,H);
            mmip3(Ko,H);
            /* Loop over all warped image pixels:
            rw/cw: 1-based row/column indices into warped image */
            for(cw=1; cw<=widwarp; ++cw) {
                for(rw=1; rw<=heiwarp; ++rw) {
                    /* Project warped image pixel (xxwarp[cw],yywarp[rw]) into original image.
                    Position in original image is (xo,yo) */
                    denom = (H[2]*cw + H[5]*rw + H[8]) + 1e-16;
                    yo    = (H[1]*cw + H[4]*rw + H[7])/denom;
                    xo    = (H[0]*cw + H[3]*rw + H[6])/denom;
                    /* If very close to a whole pixel, round to that pixel */
                    if(fabs(xo-round(xo)) < 1e-5)
                        xo = round(xo);
                    if(fabs(yo-round(yo)) < 1e-5)
                        yo = round(yo);
                    /* If at the edge, make sure we're just inside the image */
                    if(yo==heiorig)
                        yo=heiorig-1e-5;
                    if(xo==widorig)
                        xo=widorig-1e-5;
                    /* i: 0-based linear index into warped image = sub2ind([heiwarp widwarp],rw,cw); */
                    i = (mwIndex)((cw-1)*heiwarp + rw - 1);
                    /* Check if (xo,yo) lies outside domain of original image */
                    if (yo<1 || yo>=heiorig || xo<1 || xo>=widorig) {
                        if (clamp_edges_to_zero)     
                            imwarp[i] = -1e10;
                        continue;
                    } else {
                        imwarp[i] += 1;
                    }
                }
            }
        }
    } else {
        mexErrMsgTxt("Only use this for non uniform blurs");
    }
    /* Set all flagged pixels to zero */
    for(i=0; i<nwarp; ++i)
        if(imwarp[i] <= -1e5)
        imwarp[i] = 0;
        /* Assign output */
    plhs[0] = imwarp_ptr;
}


