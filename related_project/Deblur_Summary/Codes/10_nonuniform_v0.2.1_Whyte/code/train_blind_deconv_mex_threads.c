#include <pthread.h>
#include <errno.h>
/* #include <unistd.h> */ /* for system information */
#include <stddef.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"
#include "ow_homography.h"

/*
    Function to do the heavy computation for variational blind deblurring
        Multi-threaded version, requires libpthread. i.e. mex with option -lpthread

    Author:     Oliver Whyte <oliver.whyte@ens.fr>
    Date:       August 2010
    Copyright:  2010, Oliver Whyte
    Reference:  O. Whyte, J. Sivic, A. Zisserman, and J. Ponce. ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
    URL:        http://www.di.ens.fr/~whyte/deblurring/
*/


#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

double sqr(double const x) { return x*x; }

/* Global variables -- read-only */
double *Kblurry,*Ksharp;
double invKblurry[9];
int h_sharp, w_sharp, n_sharp, h_blurry, w_blurry, n_blurry, n_kernel;
double *theta_list;
double *Hcache;
double *mf1, *mf2, *mw1, *mw2, *g, *obsmask, mmu1, mmu2;
int sharp_subim_l[2];
int blurry_subim_l[2];
double *var_w, *var_f;

/* threads global variables */
struct thread_data{
    int  thread_id, top, bottom, left, right, subimage;
    double *w1w2, *w2, *f1f2, *f2, *mu1, *mu2, *Ai, *exp_sqr_err, *Bi, *coeffs;
    int *xjfloor, *yjfloor;
};

void *thread_func(void* threadarg) {
    /* Get struct argument */
    struct thread_data *my_data;
    my_data = (struct thread_data *) threadarg;
    int thread_id = my_data->thread_id;
    int top = my_data->top;
    int bottom = my_data->bottom;
    int left = my_data->left;
    int right = my_data->right;
    int subimage = my_data->subimage;
    double* Ai = my_data->Ai;
    double* w1w2 = my_data->w1w2;
    double* w2 = my_data->w2;
    double* f1f2 = my_data->f1f2;
    double* f2 = my_data->f2;
    double* mu1 = my_data->mu1;
    double* mu2 = my_data->mu2;
    double* exp_sqr_err = my_data->exp_sqr_err;
    double* Bi = my_data->Bi;
    int *xjfloor = my_data->xjfloor;
    int *yjfloor = my_data->yjfloor;
    double *coeffs = my_data->coeffs;

    double erri;
    double sqrcoeffs[n_interp];

    /* Loop over all blurry image pixels:
        (xi,yi): 1-based 2D position in blurry image */
    int xi, yi;
    int i, j, k, jj, kjj;
    double sqrBik, sqrAij, w2inc, term2i, term3i, term4i;
    double xj, yj;
    int xinterp, yinterp;
    bool inbounds;
    /* Loop over all blurry pixels i */
    for(xi=left; xi<=right; ++xi) {
        for(yi=top; yi<=bottom; ++yi) {
            i = (xi-1+blurry_subim_l[subimage]-1)*h_blurry + yi - 1;
            if(obsmask[i])
                ++(*mu2); /* Increment observation count */
            else
                continue; /* Skip this blurry pixel */
            term2i = 0;
            term3i = 0;
            term4i = 0;
            /* Clear temporary variables */
            double gbari = 0;
            /* Clear Ai */
            for(j=0; j<n_sharp; ++j)
                Ai[j] = 0;
            /* Loop over all orientations k */
            for(k=0; k<n_kernel; ++k) {
                Bi[k] = 0;
                /* Project blurry pixel into sharp image */
                project_blurry_to_sharp(yi,xi,&Hcache[k*9],&xj,&yj);
                /* Get interpolation coefficients */
                interp_coeffs(xj,yj,&xjfloor[k],&yjfloor[k],&coeffs[k*n_interp]);
                for(jj=0; jj<n_interp; ++jj) sqrcoeffs[jj] = sqr(coeffs[k*n_interp+jj]);
                /* Interpolate points */
                for(jj=0; jj<n_interp; ++jj) {
                    yinterp = yjfloor[k] + yoff[jj];
                    xinterp = xjfloor[k] + xoff[jj];
                    inbounds = yinterp >= 1 && yinterp <= h_sharp && xinterp >= 1 && xinterp <= w_sharp;
                    if(!inbounds)
                        continue;
                    j       = (xinterp-1+sharp_subim_l[subimage]-1)*h_sharp + yinterp - 1;
                    kjj     = k*n_interp+jj;
                    Bi[k]  += coeffs[kjj]*mf1[j];
                    Ai[j]  += coeffs[kjj]*mw1[k];
                    f2[j]  += sqrcoeffs[jj]*var_w[k];
                    w2inc   = sqrcoeffs[jj]*var_f[j];
                    w2[k]  += w2inc;
                    term2i += w2inc*var_w[k]; /* += sqrcoeffs[jj]*var_f[j]*var_w[k]; */
                }
                sqrBik = Bi[k]*Bi[k];
                /* w2 */
                w2[k] += sqrBik;
                /* Accumulate gbar[i] */
                    gbari += Bi[k]*mw1[k];
                /* exp_sqr_err */
                term4i += sqrBik*var_w[k];
                /* Having calculated Bik, interpolate terms involving it */
                for(jj=0; jj<n_interp; ++jj) {
                    yinterp = yjfloor[k] + yoff[jj];
                    xinterp = xjfloor[k] + xoff[jj];
                    inbounds = yinterp >= 1 && yinterp <= h_sharp && xinterp >= 1 && xinterp <= w_sharp;
                    if(!inbounds)
                        continue;
                    j = (xinterp-1+sharp_subim_l[subimage]-1)*h_sharp + yinterp - 1;
                    f1f2[j] += coeffs[k*n_interp+jj]*Bi[k]*var_w[k];
                }
            }
            /* Expected error */
            erri = obsmask[i]*(g[i] - gbari - mmu1);
            *mu1 += obsmask[i]*(g[i] - gbari);
            /* Accumulate terms for image parameters */
            for(j=0; j<n_sharp; ++j) {
                if(Ai[j] > 0) {
                    sqrAij = Ai[j]*Ai[j];
                    term3i += sqrAij*var_f[j];
                    f2[j] += sqrAij;
                    f1f2[j] += Ai[j]*erri;
                }
            }
            /* Loop over rotations again to calculate terms for kernel */
            for(k=0; k<n_kernel; ++k) {
                /* Interpolate points */
                for(jj=0; jj<n_interp; ++jj) {
                    yinterp = yjfloor[k] + yoff[jj];
                    xinterp = xjfloor[k] + xoff[jj];
                    inbounds = yinterp >= 1 && yinterp <= h_sharp && xinterp >= 1 && xinterp <= w_sharp;
                    if(!inbounds)
                        continue;
                    j = (xinterp-1+sharp_subim_l[subimage]-1)*h_sharp + yinterp - 1;
                    w1w2[k] -= coeffs[k*n_interp+jj]*Ai[j]*var_f[j];
                }
                /* w1w2 */
                w1w2[k] += Bi[k]*erri;
            }
            /* Add everything together for exp_sqr_err */
            exp_sqr_err[i] = sqr(erri) + term2i + term3i + term4i + mmu2 - sqr(mmu1);
        }
    }
    pthread_exit(NULL);
}

void errorCheck(const char *fnname, const int rc) {
    if(rc) {
        mexPrintf(fnname);
        switch(rc) {
            case ESRCH:
                mexErrMsgTxt(" failed with return code ESRCH");
                break;
            case EDEADLK:
                mexErrMsgTxt(" failed with return code EDEADLK");
                break;
            case EINVAL:
                mexErrMsgTxt(" failed with return code EINVAL");
                break;
            default:
                mexErrMsgTxt(" failed with unknown return code");
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray const* mf1_ptr        = prhs[0]; /* first moment of q(f) */
    mxArray const* mf2_ptr        = prhs[1]; /* raw second moment of q(f) */
    mxArray const* mw1_ptr        = prhs[2]; /* first moment of q(w) */
    mxArray const* mw2_ptr        = prhs[3]; /* raw second moment of q(w) */
    mxArray const* g_ptr          = prhs[4]; /* blurry image */
    mxArray const* Ksharp_ptr     = prhs[5]; /* intrinsic calibration of sharp image */
    mxArray const* Kblurry_ptr    = prhs[6]; /* Intrinsic calibration of blurry image */
    mxArray const* theta_list_ptr = prhs[7]; /* angles covered by camera kernel */
    mxArray const* obsmask_ptr    = prhs[8]; /* Mask of observed blurry pixels */
    mmu1                          = mxGetScalar(prhs[9]); 
    mmu2                          = mxGetScalar(prhs[10]);
    int num_threads               = (nrhs > 11) ? (int)mxGetScalar(prhs[11]) : 2;
        
    n_kernel = mxGetN(theta_list_ptr);

    mf1 = mxGetPr(mf1_ptr);
    mf2 = mxGetPr(mf2_ptr);
    mw1 = mxGetPr(mw1_ptr);
    mw2 = mxGetPr(mw2_ptr);
    g = mxGetPr(g_ptr);
    Ksharp = mxGetPr(Ksharp_ptr);
    Kblurry = mxGetPr(Kblurry_ptr);
    theta_list = mxGetPr(theta_list_ptr);
    obsmask = mxGetPr(obsmask_ptr);
    
    if(mxGetNumberOfDimensions(mf1_ptr) > 2) mexErrMsgTxt("Images must be grayscale");
    if(mxGetNumberOfDimensions(g_ptr) > 2) mexErrMsgTxt("Images must be grayscale");
    int const *dims_sharp = mxGetDimensions(mf1_ptr);
    h_sharp = dims_sharp[0];
    w_sharp = dims_sharp[1]/2;
    int const *dims_blurry = mxGetDimensions(g_ptr);
    h_blurry = dims_blurry[0];
    w_blurry = dims_blurry[1]/2;
    n_sharp = mxGetNumberOfElements(mf1_ptr);
    n_blurry = mxGetNumberOfElements(g_ptr);
    
    /* X and Y gradients make up 2 "subimages" in the arrays, which need to be handled separately */
    sharp_subim_l[0] = 1;
    sharp_subim_l[1] = 1+w_sharp;
    blurry_subim_l[0] = 1;
    blurry_subim_l[1] = 1+w_blurry;
    
    /* Find number of processors available */
    if(num_threads == 1) {
        mexErrMsgTxt("num_threads must be greater than 1");
/*      mexErrMsgIdAndTxt("deblur:numthreads","num_threads must be greater than 1. Your system reports having %d processors.",sysconf(_SC_NPROCESSORS_ONLN)); */
    }
    /* Number of threads must be an even number */
    num_threads = min(ceil(num_threads/2), w_blurry)*2;
    
    /* Invert calibration matrix for blurry image */
    inv3(Kblurry,invKblurry);

    int t;
    mxArray *Ai_ptr[num_threads];
    double *Ai[num_threads];
    mxArray *w1w2_ptr[num_threads];
    double *w1w2[num_threads];
    mxArray *w2_ptr[num_threads];
    double *w2[num_threads];
    mxArray *f1f2_ptr[num_threads];
    double *f1f2[num_threads];
    mxArray *f2_ptr[num_threads];
    double *f2[num_threads];
    mxArray *mu1_ptr[num_threads];
    double *mu1[num_threads];
    mxArray *mu2_ptr[num_threads];
    double *mu2[num_threads];
    mxArray *exp_sqr_err_ptr[num_threads];
    double *exp_sqr_err[num_threads];
    double *Bi[num_threads];
    int *xjfloor[num_threads];
    int *yjfloor[num_threads];
    double *coeffs[num_threads];
    for(t=0; t<num_threads; ++t) {
        /* general */
        Ai_ptr[t]          = mxCreateNumericArray(2, dims_sharp, mxDOUBLE_CLASS, mxREAL);
        Ai[t]              = mxGetPr(Ai_ptr[t]);
        /* exp_sqr_err */
        exp_sqr_err_ptr[t] = mxCreateNumericArray(2, dims_blurry, mxDOUBLE_CLASS, mxREAL);
        exp_sqr_err[t]     = mxGetPr(exp_sqr_err_ptr[t]);
        /* f1f2 */
        f1f2_ptr[t]        = mxCreateNumericArray(2, dims_sharp, mxDOUBLE_CLASS, mxREAL);
        f1f2[t]            = mxGetPr(f1f2_ptr[t]);
        /* f2 */
        f2_ptr[t]          = mxCreateNumericArray(2, dims_sharp, mxDOUBLE_CLASS, mxREAL);
        f2[t]              = mxGetPr(f2_ptr[t]);
        /* w1w2 */
        w1w2_ptr[t]        = mxCreateDoubleMatrix(n_kernel, 1, mxREAL);
        w1w2[t]            = mxGetPr(w1w2_ptr[t]);
        /* w2 */
        w2_ptr[t]          = mxCreateDoubleMatrix(n_kernel, 1, mxREAL);
        w2[t]              = mxGetPr(w2_ptr[t]);
        /* mu1             = sum(sum(Dp.*(D-mD))); */
        mu1_ptr[t]         = mxCreateDoubleMatrix(1, 1, mxREAL);
        mu1[t]             = mxGetPr(mu1_ptr[t]);
        /* mu2             = data_points */
        mu2_ptr[t]         = mxCreateDoubleMatrix(1, 1, mxREAL);
        mu2[t]             = mxGetPr(mu2_ptr[t]);
        /* Temporary storage for threads */
        Bi[t]              = (double*) mxCalloc(n_kernel, sizeof(double));
        xjfloor[t]         = (int*) mxCalloc(n_kernel, sizeof(int));
        yjfloor[t]         = (int*) mxCalloc(n_kernel, sizeof(int));
        coeffs[t]          = (double*) mxCalloc(n_kernel*n_interp, sizeof(double));
    }
    mxArray* rerror_ptr = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *rerror = mxGetPr(rerror_ptr);
    /* Hcache, used to cache homography matrices */
    mxArray* Hcache_ptr = mxCreateDoubleMatrix(9, n_kernel, mxREAL);
    Hcache = mxGetPr(Hcache_ptr);

    int i, j, k;
    /* Calculate variance of each kernel element */
    mxArray *var_w_ptr = mxCreateDoubleMatrix(1, n_kernel, mxREAL);
    var_w = mxGetPr(var_w_ptr);
    for(k=0; k<n_kernel; ++k)
        var_w[k] = mw2[k] - sqr(mw1[k]);
    /* Calculate variance of each sharp pixel */
    mxArray *var_f_ptr = mxCreateNumericArray(2, dims_sharp, mxDOUBLE_CLASS, mxREAL);
    var_f = mxGetPr(var_f_ptr);
    for(j=0; j<n_sharp; ++j)
        var_f[j] = mf2[j] - sqr(mf1[j]);
    /* Calculate homographies and cache them */
    for(k=0; k<n_kernel; ++k)
        compute_homography_matrix(Ksharp, &theta_list[k*3], invKblurry, &Hcache[k*9]);
    /* Setup multiple threads */
    pthread_t thread[num_threads];
    struct thread_data thread_data_array[num_threads];
    void *status;
    int rc=0;
    pthread_attr_t attr;

    /* Initialize and set thread detached attribute so that we can join the thread */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    /* Divide each subimage of blurry image into num_threads/2 slices */
    int thread_width = (int)ceil(((double)w_blurry)/((double)num_threads/2));
    
    /* Launch num_threads/2 threads per subimage */
    for(t=0; t<num_threads; ++t) {
        int subimage = (t < num_threads/2)? 0 : 1;
        thread_data_array[t].thread_id = t;
        thread_data_array[t].top = 1;
        thread_data_array[t].bottom = h_blurry;
        thread_data_array[t].left = t*thread_width+1-blurry_subim_l[subimage]+1;
        thread_data_array[t].right = min((t+1)*thread_width-blurry_subim_l[subimage]+1,w_blurry);
        thread_data_array[t].subimage = subimage;
        /* Assign each thread a set of arrays for its outputs */
        thread_data_array[t].Ai = Ai[t];
        thread_data_array[t].w1w2 = w1w2[t];
        thread_data_array[t].w2 = w2[t];
        thread_data_array[t].f1f2 = f1f2[t];
        thread_data_array[t].f2 = f2[t];
        thread_data_array[t].mu1 = mu1[t];
        thread_data_array[t].mu2 = mu2[t];
        thread_data_array[t].exp_sqr_err = exp_sqr_err[t];
        thread_data_array[t].Bi = Bi[t];
        thread_data_array[t].xjfloor = xjfloor[t];
        thread_data_array[t].yjfloor = yjfloor[t];
        thread_data_array[t].coeffs = coeffs[t];
        rc = pthread_create(&thread[t], &attr, thread_func, (void *) &thread_data_array[t]);
        errorCheck("pthread_create()", rc);
    }

    /* Free attribute and wait for the other threads */
    pthread_attr_destroy(&attr);

    /* Join threads to wait until they're all finished */
    for(t=0; t<num_threads; ++t) {
        rc = pthread_join(thread[t], &status);
        errorCheck("pthread_join()", rc);
        if (status != NULL) {
            mexPrintf("thread %ld exited with status %ld",t,(long)status);
            mexErrMsgTxt("\n");
        }
    }
    
    /* Combine threads' results */

    /* Add threads' results together */
    for(t=1; t<num_threads; ++t) {
        for(j=0; j<n_sharp; ++j) {
            (f1f2[0])[j] += (f1f2[t])[j];
            (f2[0])[j] += (f2[t])[j];
        }
        for(k=0; k<n_kernel; ++k) {
            (w1w2[0])[k] += (w1w2[t])[k];
            (w2[0])[k] += (w2[t])[k];
        }
        *(mu1[0]) += *(mu1[t]);
        *(mu2[0]) += *(mu2[t]);
    }
    
    /* Final addition, f1f2 += mf1.*f2 */
    for(j=0; j<n_sharp; ++j) {
        (f1f2[0])[j] += mf1[j]*(f2[0])[j];
    }
    for(k=0; k<n_kernel; ++k) {
        (w1w2[0])[k] += mw1[k]*(w2[0])[k];
    }
    *rerror = 0;
    for(t=0; t<num_threads; ++t) {
        for(i=0; i<n_blurry; ++i) {
            /* rerror is the sum of exp_sqr_err over all pixels */
            *rerror += (exp_sqr_err[t])[i];
        }
    }
    plhs[0] = f1f2_ptr[0];
    plhs[1] = f2_ptr[0];
    plhs[2] = w1w2_ptr[0];
    plhs[3] = w2_ptr[0];
    plhs[4] = mu1_ptr[0];
    plhs[5] = mu2_ptr[0];
    plhs[6] = rerror_ptr;
/*  pthread_exit(NULL); */
}
  



