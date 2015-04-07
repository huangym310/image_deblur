#include <mex.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "ow_homography.h"

/*
    Function to apply a (uniform or non-uniform) blur kernel to an image

    Author:     Oliver Whyte <oliver.whyte@ens.fr>
    Date:       August 2010
    Copyright:  2010, Oliver Whyte
    Reference:  O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
    URL:        http://www.di.ens.fr/~whyte/deblurring/
*/

#define thetax(k) (theta_list[(k)*3])
#define thetay(k) (theta_list[(k)*3+1])
#define thetaz(k) (theta_list[(k)*3+2])

double round(double val)
{
    return floor(val + 0.5);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray const* imsharp_ptr     = prhs[0]; /* sharp image */
    mxArray const* dims_blurry_ptr = prhs[1]; /* dimensions of blurry image */
    mxArray const* Ksharp_ptr      = prhs[2]; /* Intrinsic calibration of sharp image */
    mxArray const* Kblurry_ptr     = prhs[3]; /* intrinsic calibration of blurry image */
    mxArray const* theta_list_ptr  = prhs[4]; /* Orientations covered by blur kernel */
    mxArray const* kernel_ptr      = prhs[5]; /* Kernel */
    bool const clamp_edges_to_zero = (bool const)mxGetScalar(prhs[6]); /* Set blurry pixels to zero if they involve sharp pixels outside the image? */
    bool const non_uniform         = (bool const)mxGetScalar(prhs[7]); /* Non-uniform blur model? */

    int const n_kernel = mxGetN(theta_list_ptr); /* Number of kernel elements */

    double const *imsharp     = mxGetPr(imsharp_ptr);
    double const *dims_blurry = mxGetPr(dims_blurry_ptr);
    double const *Ksharp      = mxGetPr(Ksharp_ptr);
    double const *Kblurry     = mxGetPr(Kblurry_ptr);
    double const *theta_list  = mxGetPr(theta_list_ptr);
    double const *kernel      = mxGetPr(kernel_ptr);

    int const numdims     = mxGetNumberOfDimensions(imsharp_ptr);
    int const *dims_sharp = mxGetDimensions(imsharp_ptr);
    int const h_sharp     = dims_sharp[0];
    int const w_sharp     = dims_sharp[1];
    int const channels    = (numdims==3) ? dims_sharp[2] : 1;
    int const h_blurry    = (int const)dims_blurry[0];
    int const w_blurry    = (int const)dims_blurry[1];
    int const n_sharp     = h_sharp*w_sharp;
    int const n_blurry    = h_blurry*w_blurry;

    int full_dims_blurry[3];

    mxArray* imblurry_ptr;
    double *imblurry;

    int xi, yi;
    double yj, xj;
    int yjfloor, xjfloor;
    double H[9];
    double R[9];
    double coeffs[4];
    double invKblurry[9];

    int i, j, c, k;

    /* Compute offsets for each channel, so we're not constantly computing c*height*width */
    int* chan_offset_sharp = (int *)mxMalloc(sizeof(int) * channels);
    int* chan_offset_blurry = (int *)mxMalloc(sizeof(int) * channels);
    for(c=0; c<channels; ++c) {
        chan_offset_sharp[c]  = c*h_sharp*w_sharp;
        chan_offset_blurry[c] = c*h_blurry*w_blurry;
    }

    full_dims_blurry[0] = h_blurry;
    full_dims_blurry[1] = w_blurry;
    full_dims_blurry[2] = channels;
    imblurry_ptr = mxCreateNumericArray(3, full_dims_blurry, mxDOUBLE_CLASS, mxREAL);
    imblurry = mxGetPr(imblurry_ptr);

    /* Check thetas -- for uniform blur theta must contain only integers */
    if(!non_uniform) {
        for(k=0; k<n_kernel; ++k)
            if(fabs(round(thetay(k))-thetay(k)) > 1e-10 ||     \
               fabs(round(thetax(k))-thetax(k)) > 1e-10)
               mexErrMsgTxt("theta not integral for uniform blur");
    }

    inv3(Kblurry,invKblurry);
    /* Loop over all rotations */
    for(k=0; k<n_kernel; ++k) {
        if (kernel[k]==0) continue;
        /* Loop over all blurry image pixels:
            (xi,yi): 1-based 2D position in blurry image */
        for(xi=1; xi<=w_blurry; ++xi) {
            for(yi=1; yi<=h_blurry; ++yi) {
                /* i: 0-based linear index into blurry image */
                i = (mwIndex)((xi-1)*h_blurry + yi - 1);
                /* Project blurry image pixel (xi,yi) into sharp image.
                    1-based 2D position in sharp image is (xj,yj) */
                if (non_uniform) {
                    compute_homography_matrix(Ksharp, &theta_list[k*3], invKblurry, H);
                    project_blurry_to_sharp(yi, xi, H, &xj, &yj);
                } else {
                    xj = xi + round(thetax(k));
                    yj = yi + round(thetay(k));
                }
                /* Check if (xj,yj) lies outside domain of sharp image */
                if (yj<1 || yj>=h_sharp || xj<1 || xj>=w_sharp) {
                    /* Flag the pixel as having been influenced by edge effects */
                    if (clamp_edges_to_zero)
                        for(c=0; c<channels; ++c)
                            imblurry[i+chan_offset_blurry[c]] = -1e10;
                    continue;
                }
                /* Get floor(xj,yj) and interpolation coefficients */
                if(non_uniform) {
                    interp_coeffs(xj,yj,&xjfloor,&yjfloor,coeffs);
                } else {
                    xjfloor = (int)xj;
                    yjfloor = (int)yj;
                }
                /* j: 0-based linear index into sharp image */
                j = (mwIndex)((xjfloor-1)*h_sharp + yjfloor - 1);
                /* Do interpolation */
                if(non_uniform) {
                    for(c=0; c<channels; ++c)
                        imblurry[i+chan_offset_blurry[c]] += kernel[k]*(coeffs[0]*imsharp[j+chan_offset_sharp[c]] + \
                                                                        coeffs[1]*imsharp[j+1+chan_offset_sharp[c]] + \
                                                                        coeffs[2]*imsharp[j+h_sharp+chan_offset_sharp[c]] + \
                                                                        coeffs[3]*imsharp[j+h_sharp+1+chan_offset_sharp[c]]);
                } else {
                    for(c=0; c<channels; ++c)
                        imblurry[i+chan_offset_blurry[c]] += kernel[k]*imsharp[j+chan_offset_sharp[c]];
                }
            }
        }
    }
    /* Set all flagged pixels to zero */
    for(i=0; i<n_blurry; ++i)
        for(c=0; c<channels; ++c)
            if(imblurry[i+chan_offset_blurry[c]] <= -1e5)
                imblurry[i+chan_offset_blurry[c]] = 0;
    /* Assign output */
    plhs[0] = imblurry_ptr;

    mxFree(chan_offset_sharp);
    mxFree(chan_offset_blurry);
}


