Title:      Matlab code for "Non-uniform Deblurring for Shaken Images"  
Author:     Oliver Whyte <oliver.whyte@ens.fr>  
Version:    0.2.1
Date:       September 30, 2010  
Copyright:  2010, Oliver Whyte


Matlab code for "Non-uniform Deblurring for Shaken Images"
==========================================================

This package contains code to perform blind deblurring of non-uniform /
spatially-varying blur caused by camera shake, using the model described in
[][#Whyte10], applied within the algorithm described in [][#Fergus06] and
[][#Miskin00]. Please cite these three papers if using this code in an 
academic publication.

Please send bug reports to <oliver.whyte@ens.fr>

[#Whyte10]: O. Whyte, J. Sivic, A. Zisserman, and J. Ponce.
``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.

[#Fergus06]: R. Fergus, B. Singh, A. Hertzmann, S. T. Roweis, and W.
T. Freeman. ``Removing camera shake from a single photograph''. ACM
Trans. Graphics (Proc. SIGGRAPH 2006), 2006.

[#Miskin00]: J. W. Miskin and D. J. C. MacKay. ``Ensemble Learning
for Blind Image Separation and Deconvolution''. In Advances in
Independent Component Analysis, 2000.


## 1 Installing ##

* Compile the mex files. In Matlab, `cd` to the `code` directory, and run
  `mex_all`
* Add the `code` directory to your Matlab path

### 1.1 Compatibility ###

* From version 0.2, the mex files should compile and run on Windows. This has
  not yet been extensively tested however.
* Multi-threaded version requires libpthread to be available (tested on 64-bit
  Linux Ubuntu 9.10 with Matlab R2007b & R2009a, and 32-bit Mac OS X 10.5 with
  Matlab R2006b)
    * For Windows, you might try 
      [pthreads-win32](http://sourceware.org/pthreads-win32/), although I have
      no experience with this myself. If you have any success with this, 
      please do let me know


## 2 Running ##

* In Matlab, `cd` to `results` directory
* Ensure that the `code` directory is on the Matlab path
* Run one of the scripts called `deblur_...`
* To create scripts for your own images, the existing scripts should provide a
  good start
    * A basic `deblur_...` script has a few components:
        1. Set `CONFIG_FNAME` for the name of the directory where the outputs 
           will go
        2. Set `NON_UNIFORM = 1 / 0` for non-uniform / uniform blur model
        3. Run an `image_settings_...` script, which contains per-image
           settings
        4. Run the `default_config` script, which sets a lot of parameters to
           their default values
        5. Run the `deblur` script, which runs the actual algorithm
* Note that the algorithm can take a long time to run with the non-uniform
  blur model, e.g. the included `pantheon` example takes 3 1/2 hours on my 
  workstation, with `NUM_THREADS = 4` (see below)

### 2.1 Differences with the code from Fergus et al. ###

* The code is not compatible with the original codes from Miskin & MacKay or
  Fergus et al., as many of the functions have been restructured and take
  different arguments. However, you should get similar, if not identical,
  results to that of Fergus et al., when using the uniform blur model

### 2.2 A few significant options ###

Please see the original [readme from Fergus et al.'s code](Fergus_readme.txt) 
for an explanation of more options, and more information on the algorithm in 
general.

`PRESCALE`
:   Factor to downsample original blurry image.

`BLUR_KERNEL_SIZE`
:   Size of blur in blurry image after downsampling by `PRESCALE`.

`blur_{x,y,z}_lims = [min, max]`
:   Maximum extent of blur due to contributions of <<theta_X>>, <<theta_Y>>,
    and <<theta_Z>>. Usually derived from `BLUR_KERNEL_SIZE`.

`focal_length_35mm`
:   Camera's focal length in 35mm equivalent. This is sometimes available
    directly in the EXIF tags of an image. If not, the actual focal length
    should be in the EXIF tags. With this, you will need to know the size of
    the camera's sensor in mm (you might try looking on
    <http://www.dpreview.com/reviews/specs.asp>). Then,   
    `focal_length_35mm = focal_length / sensor_width_in_mm * 36`.

`AXIS = [xmin, xmax, ymin, ymax]`
:   Region of (downsampled) blurry image to use for kernel estimation.

`FIRST_INIT_MODE_BLUR = 'xbar' / 'ybar' / 'zbar'`
:   Initialisation of kernel, depending on approximate shape of blur. One of 3
    axes: `xbar` for approximately vertical blur, `ybar` for horizontal, or
    `zbar` for in-plane rotation.

`DISPLAY_EACH_ITERATION = true / false`
:   After each iteration, show the latent image and kernel, and plot the value
    of the cost function over all iterations.

`SAVE_EACH_ITERATION = true / false`
:   Save images of the latent image and kernel after each iteration. Note that
    this option will cause hundreds of (small) image files to be saved in the
    results directory.

`PLOT_LINESEARCH = true / false`
:   At each iteration, show value of cost function at different points along
    search direction.

`NUM_THREADS`
:   On systems where libpthread is available, the non-uniform blur model can
    be run multi-threaded for speed. On other systems, `NUM_THREADS` must be
    set to 1.


## 3 Acknowledgements ##

We would like to thank Rob Fergus for making his code available 
online (<http://cs.nyu.edu/~fergus/research/deblur.html>) and for kindly 
allowing us to release our modified version of his code, as well as James 
Miskin and David MacKay for also making their original code available 
online (<http://www.inference.phy.cam.ac.uk/jwm1003/>).

Thanks to Xianyong Fang for testing and debugging the mex files on Windows.


## 4 License ##

The code in this package is based on the code provided by Fergus et al., and 
as such several parts are subject to their own license. For functions 
originally distributed by Fergus et al., please see the original [readme from 
Fergus et al.'s code](Fergus_readme.txt) for details. These functions are 
marked with the header:

	Author: Rob Fergus (or other)
	Version: 1.0, distribution code.
	Project: Removing Camera Shake from a Single Image, SIGGRAPH 2006 paper
	Copyright 2006, Massachusetts Institute of Technology
	
For other functions marked with the header:

	Author:		Oliver Whyte <oliver.whyte@ens.fr>
	Date:		August 2010
	Copyright:	2010, Oliver Whyte
	Reference:	O. Whyte, J. Sivic, A. Zisserman, and J. Ponce.
	   ``Non-uniform Deblurring for Shaken Images''. In Proc. CVPR, 2010.
	URL:		http://www.di.ens.fr/~whyte/deblurring/

the following license applies:

* * *

Copyright (c) 2010, Oliver Whyte

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the ``Software''), to 
deal in the Software without restriction, including without limitation the 
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED ``AS IS'', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

* * *


## 5 Changelog ##

### Version 0.2.1 (30-Sep-2010) ###

- Modified mex files to allow compilation on Windows

### Version 0.2 (23-Sep-2010) ###

- Added Xianyong Fang's modifications to mex files for Windows
- Corrected sign of theta in fiddle_lucy3_rot.m

### Version 0.1 (30-Aug-2010) ###

- Initial release

