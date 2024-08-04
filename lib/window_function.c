/*******************************************************************************
	window_function.c -- Window Function
	
	Author: Hironobu Inatsuka
*******************************************************************************/
#include <ruby.h>
//#include <ruby/internal/memory.h> // ALLOC_N()
//#include <ruby/internal/intern/array.h> // rb_ary_new(), rb_ary_store()
#include "ruby/wave/globals.h"
#include "internal/algorithm/wf.h"

#ifndef HAVE_CYL_BESSEL_I0
double cyl_bessel_i0(double);
# include "missing/cyl_bessel_i0.c"
#endif


/*
 *  module Wave::WindowFunction
 *  
 *  - Overview -
 *  Wave::WindowFunction is modularizes the window functions that used at the C level.
 *  The implementation is a discrete type often used in waveform filtering.
 *  DSP programming is multi-threaded, and the implementation like this is called "user-level" as opposed to kernel-level.
 *  The implementation is more focused on test suites as a collection of algorithms than on execution speed; 
 *  This module is unlikely to be used directly in DSP applications.
 *  
 *  A callback method is used for array generation.
 *  The number of window function arrays is no less than 100 per generates at once. 
 *  If the master frequency is 96kHz, you will generate twice that number of arrays;
 *  This problem is solved by incorporating callback function into the iteration to speed up processing.
 *  
 *  The design philosophy is based on Dr. Naofumi Aoki.
 *  
 *  By the way, the continuous window function uses the x-axis as a variable with 0 as the median value, 
 *  and the domain is do: $-\frac{1}{2}\leq{x}\leq\frac{1}{2}$
 *  The discrete type has a different approach in that it limits the domain of x to "$0 \leq x \leq 1$"
 *  and uses a discrete signal with the number of arrays as quartiles.
 *  
 *    ```Ruby
 *    def cont_hann(x)
 *      if (-0.5 <= x && x <= 0.5)
 *        0.5 + 0.5 * Math.cos(2 * Math::PI * x)
 *      else
 *        0.0
 *      end
 *    end
 *    
 *    len = 5
 *    len.times.map{|x| x = (1.0 * x / len) - (2.0 / len); cont_hann(x)}
 *    # => [0.09549150281252633,
 *    # =>  0.6545084971874737,
 *    # =>  1.0,
 *    # =>  0.6545084971874738,
 *    # =>  0.09549150281252633]
 *    
 *    Wave::WindowFunction.hann(len)
 *    # => [0.09549150281252627,
 *    # =>  0.6545084971874737,
 *    # =>  1.0,
 *    # =>  0.6545084971874737,
 *    # =>  0.09549150281252633]
 *    ```
 *  
 *  - About algorithm -
 *  This implementation employs an iterative process due to the discrete nature of the algorithm.
 *  The iterator has the modified discrete cosine transform (so-called MDCT, which is used in audio processing), and the one-dimensional rule.
 *  "one-dimensional" is so named in contrast to the discrete cosine transform which has poles.
 *  Note that inverse functions also exist for window functions (MDCT windows also have them), 
 *  however, since this module is a subclass of waveform, it does not contain these inverse algorithms.
 *  
 */



static inline VALUE
rb_wf_ary_new(void (*func)(double, long, double *), long len, double param)
{
	VALUE ary = rb_ary_new2(len);
	double *w = ALLOC_N(double, len);
	
	func(param, len, w);
	
	for (volatile long n = 0; n < len; n++)
		rb_ary_store(ary, n, DBL2NUM(w[n]));
	
	return ary;
}

/*******************************************************************************
	Dirichlet Window / Rectangular Window
*******************************************************************************/
#include "internal/solver/window_function/rectangular.h"

static void
wf_cb_rectangular(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_rectangular_expr,
		0.,
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.dirichlet(len) -> [*Float]
 *    Wave::WindowFunction.rectangular(len) -> [*Float]
 *  
 *  Returns an array of discrete rectangular windows. Specify the number of arrays with len.
 *  In Europe, it is well known as the "Dirichlet window" and is one of the most commonly used window functions. 
 *  It is always a scalar quantity of 1.0.
 *  
 *    Wave::WindowFunction.rectangular(5)
 *    # => [1.0, 1.0, 1.0, 1.0, 1.0]
 */
static VALUE
wf_rectangular(VALUE unused_obj, VALUE len)
{
	return rb_wf_ary_new(wf_cb_rectangular, NUM2LONG(len), 0.);
}


/*******************************************************************************
	Hamming Window / Generalized Hamming Window
*******************************************************************************/
#include "internal/solver/window_function/hamming.h"

static void
wf_cb_hamming(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_hamming_expr, 
		0.,
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

#include "internal/solver/window_function/generalized_hamming.h"

static void
wf_cb_generalized_hamming(double alpha, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_generalized_hamming_expr, 
		wf_generalized_hamming_calc_param(alpha),
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.hamming(len) -> [*Float]
 *    Wave::WindowFunction.hamming(len, alpha) -> [*Float]
 *  
 *  Returns an array of discrete Hamming windows. Specify the number of arrays with len.
 *  Hamming window is one of the commonly used window functions.
 *  A discrete Hamming window is usually defined as:
 *  $ w(x)=\frac{25}{46}-\frac{21}{46}\cos(2\pi{x}), 0 \leq x \leq 1 $
 *  
 *  When $\alpha$ is assigned as the second argument, returns an array of discrete generalized Hamming windows for the number of arrays len.
 *  The generalized Hamming window is a generalization of the Hann window and the Hamming window, 
 *  and usually takes the real parameter $\alpha$ as the domain $\frac{1}{2}\leq\alpha 1$.
 *  
 *  If $\alpha$ is a value outside this domain, a RangeError exception will occur. 
 *  By the way, ambiguity is a strength in electrical mathematics, so it is not appropriate for exception handling to work. 
 *  This occurs because the value obtained is far from the expected value.
 *  
 *  A discrete generalized Hamming window is usually defined as:
 *  $ w(x)=\alpha-(1-\alpha)\cos(2\pi{x}), 0 \leq x \leq 1 $
 *  
 *  
 *    Wave::WindowFunction.hamming(5)
 *    # => [0.174144415611437,
 *    # =>  0.684551236562476,
 *    # =>  1.0,
 *    # =>  0.684551236562476,
 *    # =>  0.17414441561143706]
 *    Wave::WindowFunction.hamming(5, 25/46r)
 *    # => [0.17414441561143695,
 *    # =>  0.684551236562476,
 *    # =>  1.0,
 *    # =>  0.684551236562476,
 *    # =>  0.17414441561143695]
 *    
 *    Wave::WindowFunction.hamming(5, 1)
 *    # => [1.0, 1.0, 1.0, 1.0, 1.0]
 *    Wave::WindowFunction.rectangular(5)
 *    # => [1.0, 1.0, 1.0, 1.0, 1.0]
 *    
 *    Wave::WindowFunction.hamming(5, 0) # => RangeError
 */
static VALUE
wf_hamming(int argc, VALUE *argv, VALUE unused_obj)
{
	VALUE len, param = Qnil;
	rb_scan_args(argc, argv, "11", &len, &param);
	if (argc == 1)
	{
		return rb_wf_ary_new(wf_cb_hamming, NUM2LONG(len), 0.);
	}
	else
	{
		return rb_wf_ary_new(wf_cb_generalized_hamming, NUM2LONG(len), NUM2DBL(param));
	}
}


/*******************************************************************************
	Hann window / Parameterized Hann window
*******************************************************************************/
#include "internal/solver/window_function/hann.h"

static void
wf_cb_hann(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_hann_expr,
		0.,
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.hann(len) -> [*Float]
 *    Wave::WindowFunction.hanning(len) -> [*Float]
 *    Wave::WindowFunction.hann(len, alpha) -> [*Float]
 *    Wave::WindowFunction.hanning(len, alpha) -> [*Float]
 *  
 *  Returns an array of discrete Hann windows. Specify the number of arrays with len.
 *  The Hann window is one of the commonly used window functions.
 *  This function is also called the 'Hanning window' after the parameter-modified Hamming window.
 *  A discrete Hann window is usually defined as:
 *  $ w(x)=\frac{1}{2}-\frac{1}{2}\cos(2\pi{x}), 0 \leq x \leq 1 $
 *  Where, the coefficient $\alpha=\frac{1}{2}$ is related to the order $1-\alpha$ on the cosine term.
 *  
 *  If the second argument alpha is given, the return value is a parameterized Hann window.
 *  $\alpha$ is a real parameter whose domain is $\frac{1}{2}\leq\alpha 1$.
 *  
 *  If $\alpha$ is a value outside this domain, a RangeError exception will occur. 
 *  By the way, ambiguity is a strength in electrical mathematics, so it is not appropriate for exception handling to work. 
 *  This occurs because the value obtained is far from the expected value.
 *  
 *    Wave::WindowFunction.hann(5)
 *    # => [0.09549150281252627,
 *    # =>  0.6545084971874737,
 *    # =>  1.0,
 *    # =>  0.6545084971874737,
 *    # =>  0.09549150281252633]
 *    Wave::WindowFunction.hann(5, 0.5)
 *    # => [0.09549150281252627,
 *    # =>  0.6545084971874737,
 *    # =>  1.0,
 *    # =>  0.6545084971874737,
 *    # =>  0.09549150281252627]
 *    
 *    Wave::WindowFunction.hann(5, 1)
 *    # => [1.0, 1.0, 1.0, 1.0, 1.0]
 *    Wave::WindowFunction.rectangular(5)
 *    # => [1.0, 1.0, 1.0, 1.0, 1.0]
 *    
 *    Wave::WindowFunction.hann(5, 0) # => RangeError
 */
static VALUE
wf_hann(int argc, VALUE *argv, VALUE unused_obj)
{
	VALUE len, param = Qnil;
	rb_scan_args(argc, argv, "11", &len, &param);
	if (argc == 1)
	{
		return rb_wf_ary_new(wf_cb_hann, NUM2LONG(len), 0.);
	}
	else
	{
		return rb_wf_ary_new(wf_cb_generalized_hamming, NUM2LONG(len), NUM2DBL(param));
	}
}


/*******************************************************************************
	Bartlett Window
*******************************************************************************/
#include "internal/solver/window_function/bartlett.h"

static void
wf_cb_bartlett(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_bartlett_expr,
		0.,
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.bartlett(len) -> [*Float]
 *  
 *  Returns an array of discrete Bartlett windows. Specify the number of arrays with len.
 *  The Bartlett window, also known as the triangular window, is a window function that often appears in reference books.
 *  The definition can be expressed as follows:
 *  $w(x)=1-2\left| x-\frac{1}{2} \right|, 0 \leq x \leq 1$
 
 *    Wave::WindowFunction.bartlett(5)
 *    # => [0.19999999999999996, 
 *    # =>  0.6,
 *    # =>  1.0, 
 *    # =>  0.6, 
 *    # =>  0.19999999999999996]
 */
static VALUE
wf_bartlett(VALUE unused_obj, VALUE len)
{
	return rb_wf_ary_new(wf_cb_bartlett, NUM2LONG(len), 0.);
}

/*******************************************************************************
	Blackman Window
*******************************************************************************/
#include "internal/solver/window_function/blackman.h"

static void
wf_cb_blackman(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_blackman_expr,
		0.,
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.blackman(len) -> [*Float]
 *  
 *  Returns an array of discrete Blackman windows. Specify the number of arrays with len.
 *  The Blackman window is a commonly used window function.
 *  The definition can be expressed as follows:
 *  $w(x)=0.42-0.5\cos(2\pi x)+0.08\cos(4\pi x), 0 \leq x \leq 1$
 *  
 *    Wave::WindowFunction.blackman(5)
 *    # => [0.040212862362522056,
 *    # =>  0.5097871376374778,
 *    # =>  1.0,
 *    # =>  0.5097871376374778,
 *    # =>  0.040212862362522056]
 */
static VALUE
wf_blackman(VALUE unused_obj, VALUE len)
{
	return rb_wf_ary_new(wf_cb_blackman, NUM2LONG(len), 0.);
}

/*******************************************************************************
	Gaussian Window
*******************************************************************************/
#include "internal/solver/window_function/gaussian.h"

static void
wf_cb_gaussian(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_gaussian_expr, 
		0.,
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

#include "internal/solver/window_function/gaussian_with_param.h"

static void
wf_cb_gaussian_with_param(double sigma, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_gaussian_with_param_expr, 
		wf_gaussian_calc_param(sigma),
		WFIF_ITER_1D,
		WFIF_KURT,
		WFIF_NOCNTL,
		WFIF_KURT
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.gaussian(len) -> [*Float]
 *    Wave::WindowFunction.gaussian(len, sigma) -> [*Float]
 
 *  Returns an array of discrete Gaussian windows. Specify the number of arrays with len.
 *  In general, a discrete Gaussian window satisfies the following equation:
 *  $w(x)=e^{-((-1+2x)^2/8\sigma^2)}, 0 \leq{x}\leq 1$
 *  Where $\sigma$ is the standard deviation. 
 *  If the standard deviation is $3/10$, the following equation is created.
 *  $w(x) = w(x, 3/10)$
 
 *    Wave::WindowFunction.gaussian(5)
 *    # => [0.4111122905071874,
 *    # =>  0.8007374029168081,
 *    # =>  1.0,
 *    # =>  0.8007374029168082,
 *    # =>  0.4111122905071874]
 *    Wave::WindowFunction.gaussian(5, 3/10r)
 *    # => [0.41111229050718734,
 *    # =>  0.8007374029168081,
 *    # =>  1.0,
 *    # =>  0.8007374029168082,
 *    # =>  0.41111229050718734]
 */
static VALUE
wf_gaussian(int argc, VALUE *argv, VALUE unused_obj)
{
	VALUE len, param = Qnil;
	rb_scan_args(argc, argv, "11", &len, &param);
	if (argc == 1)
	{
		return rb_wf_ary_new(wf_cb_gaussian, NUM2LONG(len), 0.);
	}
	else
	{
		return rb_wf_ary_new(wf_cb_gaussian_with_param, NUM2LONG(len), NUM2DBL(param));
	}
}

/*******************************************************************************
	Kaiser Window
*******************************************************************************/
#include "internal/solver/window_function/kaiser.h"

static void
wf_cb_kaiser(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_kaiser_expr, 
		0.,
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

#include "internal/solver/window_function/kaiser_with_param.h"

static void
wf_cb_kaiser_with_param(double alpha, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_kaiser_with_param_expr, 
		alpha,
		WFIF_ITER_1D,
		WFIF_KURT,
		WFIF_KURT,
		WFIF_RECT
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.kaiser(len) -> [*Float]
 *    Wave::WindowFunction.kaiser(len, alpha) -> [*Float]
 *  
 *  Returns an array of discrete Kaiser windows. Specify the number of arrays with len.
 *  The Kaiser window, also known as the Kaiser-Bessel window, is commonly used in finite impulse response {FIR} filter design and spectral analysis.
 *  The discrete Kaiser window is given by
 *  $ w(x)=\frac{I_0(\alpha 2 \sqrt{-(x-1)x}}{I_0(\alpha)} $
 *  We obtain.
 *  Where $I_n(x)$ is the first kind of modified Bessel function, n is the zeroth order, and $\alpha$ is the shape parameter.
 *  The following are equivalent:
 *  $ w(x) = w(x, 3) $
 *  
 *    Wave::WindowFunction.kaiser(5)
 *    # => [0.4076303841265242,
 *    # =>  0.8184078580166961,
 *    # =>  1.0,
 *    # =>  0.8184078580166961,
 *    # =>  0.4076303841265242]
 *    Wave::WindowFunction.kaiser(5, 3)
 *    # => [0.4076303841265242,
 *    # =>  0.8184078580166961,
 *    # =>  1.0,
 *    # =>  0.8184078580166961,
 *    # =>  0.4076303841265242]
 */
static VALUE
wf_kaiser(int argc, VALUE *argv, VALUE unused_obj)
{
	VALUE len, param = Qnil;
	rb_scan_args(argc, argv, "11", &len, &param);
	if (argc == 1)
	{
		return rb_wf_ary_new(wf_cb_kaiser, NUM2LONG(len), 0.);
	}
	else
	{
		return rb_wf_ary_new(wf_cb_kaiser_with_param, NUM2LONG(len), NUM2DBL(param));
	}
}

/*******************************************************************************
	Bartlett - Hann Window
*******************************************************************************/
#include "internal/solver/window_function/bartlett_hann.h"

static void
wf_cb_bartlett_hann(double unused_obj, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_bartlett_hann_expr, 
		0.,
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.bartlett_hann(len) -> [*Float]
 *  
 *  Returns an array of discrete (modified) Bartlett-Hann windows. Specify the number of arrays with len.
 *  The modified Bartlett-Hann window is
 *  $ w(x)=0.62-0.48 |x - 0.5| + 0.38 \cos(2\pi(x - 0.5)) , 0 \leq x \leq 1 $
 *  is defined as. where $|x|$ is the absolute value of x.
 *  
 *    Wave::WindowFunction.bartlett_hann(5)
 *    # => [0.12057354213751997,
 *    # =>  0.6414264578624801,
 *    # =>  1.0,
 *    # =>  0.6414264578624801,
 *    # =>  0.12057354213751997]
 */
static VALUE
wf_bartlett_hann(VALUE unused_obj, VALUE len)
{
	return rb_wf_ary_new(wf_cb_bartlett_hann, NUM2LONG(len), 0.);
}

/*******************************************************************************
	Blackman-Harris window
*******************************************************************************/
#include "internal/solver/window_function/blackman_harris.h"

static void
wf_cb_blackman_harris(double unused_obj, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_blackman_harris_expr, 
		0.,
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.blackman_harris(len) -> [*Float]
 *  
 *  Returns an array of discrete Blackman-Harris windows. Specify the number of arrays with len.
 *  In general, the Blackman-Harris window is given by the following formula
 *  $ w(x)=a_0-a_1 cos(2\pi x) + a_2 cos(4\pi x) - a_3 cos(6\pi x), 0 \leq x \leq 1 $
 *  defined by it as a minimum of four terms.
 *  Where, $a_n$ is the coefficient of the smallest four terms, and the mean and median of the following are $\frac{1}{4}$.
 *  $ \begin{array}{rcl} 
 *     a_0 & = & \frac{35875}{100000} \\ 
 *     a_1 & = & \frac{48829}{100000} \\ 
 *     a_2 & = & \frac{14128}{100000} \\ 
 *     a_3 & = & \frac{1168}{100000} 
 *    \end{array} $
 *  
 *    Wave::WindowFunction.blackman_harris(5)
 *    # => [0.010982331276248888,
 *    # =>  0.3858926687237511,
 *    # =>  1.0,
 *    # =>  0.3858926687237511,
 *    # =>  0.010982331276248888]
 */
static VALUE
wf_blackman_harris(VALUE unused_obj, VALUE len)
{
	return rb_wf_ary_new(wf_cb_blackman_harris, NUM2LONG(len), 0.);
}

/*******************************************************************************
	Nuttall Window
*******************************************************************************/
#include "internal/solver/window_function/nuttall.h"

static void
wf_cb_nuttall(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_nuttall_expr, 
		0.,
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.nuttall(len) -> [*Float]
 *  
 *  Returns an array of discrete Nuttall windows. Specify the number of arrays with len.
 *  The Nuttall window is a four-term symmetric Blackman-Harris window as
 *  $ w(x)=a_0-a_1 cos(2\pi x) + a_2 cos(4\pi x) - a_3 cos(6\pi x), 0 \leq x \leq 1 $
 *  is the L point according to Nuttall's definition.
 *  Where, $a_n$ is a coefficient and is set to 4 or less.
 *  $ \begin{array}{rcl} 
 *     a_0 & = & \frac{88942}{250000} \\ 
 *     a_1 & = & \frac{121849}{250000} \\ 
 *     a_2 & = & \frac{36058}{250000} \\ 
 *     a_3 & = & \frac{3151}{250000} 
 *    \end{array} $
 *  
 *    Wave::WindowFunction.nuttall(5)
 *    # => [0.009921342339417317,
 *    # =>  0.37949865766058255,
 *    # =>  1.0,
 *    # =>  0.37949865766058255,
 *    # =>  0.009921342339417317]
 */
static VALUE
wf_nuttall(VALUE unused_obj, VALUE len)
{
	return rb_wf_ary_new(wf_cb_nuttall, NUM2LONG(len), 0.);
}

/*******************************************************************************
	Blackman-Nutall window
*******************************************************************************/
#include "internal/solver/window_function/blackman_nuttall.h"

static void
wf_cb_blackman_nuttall(double unused_obj, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_blackman_nuttall_expr, 
		0.,
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.blackman_nuttall(len) -> [*Float]
 *  
 *  Returns an array of discrete Blackman-Nuttal windows. Specify the number of arrays with len.
 *  In general, the Blackman-Nuttal window is expressed by the following formula:
 *  $ w(x)=a_0-a_1 cos(2\pi x) + a_2 cos(4\pi x) - a_3 cos(6\pi x), 0 \leq x \leq 1 $
 *  is defined as.
 *  Where, $a_n$ are the coefficients of the smallest four terms and have the following values:
 *  $ \begin{array}{rcl} 
 *     a_0 & = & \frac{3635819}{10000000} \\ 
 *     a_1 & = & \frac{4891775}{10000000} \\ 
 *     a_2 & = & \frac{1365995}{10000000} \\ 
 *     a_3 & = & \frac{106411}{10000000} 
 *    \end{array} $
 *  
 *    Wave::WindowFunction.blackman_nuttall(5)
 *    # => [0.013328836896113066,
 *    # =>  0.3956259131038869,
 *    # =>  1.0,
 *    # =>  0.3956259131038869,
 *    # =>  0.013328836896113066]
 */
static VALUE
wf_blackman_nuttall(VALUE unused_obj, VALUE len)
{
	return rb_wf_ary_new(wf_cb_blackman_nuttall, NUM2LONG(len), 0.);
}

/*******************************************************************************
	Flat Top Window
*******************************************************************************/
#include "internal/solver/window_function/flat_top.h"

static void
wf_cb_flat_top(double alpha, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_flat_top_expr, 
		0.,
		WFIF_ITER_1D,
		WFIF_NOCNTL,
		WFIF_NOCNTL,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.flat_top(len) -> [*Float]
 *  
 *  Returns an array of discrete flat-top windows. Specify the number of arrays with len.
 *  In general, flat top windows are defined by the following formula
 *  $w(x)=a_0-a_1\cos(2\pi x)+a_2\cos(4\pi x)-a_3\cos(6\pi x)+a_4\cos(8\pi x), 0\leq x\leq1$
 *  is defined as.
 *  Where, $a_n$ are coefficients, and their values ​​are as follows:
 *  $ \begin{array}{rcl} 
 *     a_0 & = & \frac{215578947}{1000000000} \\ 
 *     a_1 & = & \frac{416631580}{1000000000} \\ 
 *     a_2 & = & \frac{277263158}{1000000000} \\ 
 *     a_3 & = & \frac{83578947}{1000000000} \\ 
 *     a_4 & = & \frac{6947368}{1000000000}
 *    \end{array} $
 *  
 *    Wave::WindowFunction.flat_top(5)
 *    # => [-0.015597277660432994,
 *    # =>  0.054544645160432864,
 *    # =>  1.0,
 *    # =>  0.054544645160432864,
 *    # =>  -0.015597277660432994]
 */
static VALUE
wf_flat_top(VALUE unused_obj, VALUE len)
{
	return rb_wf_ary_new(wf_cb_flat_top, NUM2LONG(len), 0.);
}



/*******************************************************************************
	KBD Window (Kayser-Bessel derived window)
*******************************************************************************/
#include "internal/solver/window_function/kbd_with_param.h"

static void
wf_cb_kbd_with_param(double alpha, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_kbd_with_param_expr, 
		alpha,
		WFIF_ITER_MDCT,
		WFIF_RECT,
		WFIF_RECT,
		WFIF_NOCNTL
	};
	wf_iter_cb(wfif, len, w);
}

/*
 *  call-seq:
 *    Wave::WindowFunction.kbd(x, alpha) -> [*Float]
 *    Wave::WindowFunction.kaiser_bessel_derived(x, alpha) -> [*Float]
 *  
 *  Returns an array of discrete KBD windows. Specify the number of arrays with len.
 *  KBD window is an acronym for Kaiser-Bessel-derived window, 
 *  a variation of the Kaiser window designed for use in the modified discrete cosine transform (MDCT).
 *  
 *    Wave::WindowFunction.kbd(5, 3)
 *    # => [0.4114947429371883,
 *    # =>  0.9996957233074878,
 *    # =>  1.0,
 *    # =>  0.9996957233074878,
 *    # =>  0.4114947429371883]
 *    # 
 */
static VALUE
wf_kbd(int argc, VALUE *argv, VALUE unused_obj)
{
	VALUE len, param;
	rb_scan_args(argc, argv, "20", &len, &param);
	
	return rb_wf_ary_new(wf_cb_kbd_with_param, NUM2LONG(len), NUM2DBL(param));
}


/******************************************************************************/

// Entry Point
void
InitVM_WindowFunction(void)
{
	rb_define_module_function(rb_mWaveWindowFunction, "rectangular", wf_rectangular, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "dirichlet", wf_rectangular, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "hann", wf_hann, -1);
	rb_define_module_function(rb_mWaveWindowFunction, "hanning", wf_hann, -1);
	rb_define_module_function(rb_mWaveWindowFunction, "hamming", wf_hamming, -1);
	rb_define_module_function(rb_mWaveWindowFunction, "bartlett", wf_bartlett, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "blackman", wf_blackman, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "gaussian", wf_gaussian, -1);
	rb_define_module_function(rb_mWaveWindowFunction, "kaiser", wf_kaiser, -1);
	rb_define_module_function(rb_mWaveWindowFunction, "bartlett_hann", wf_bartlett_hann, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "nuttall", wf_nuttall, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "blackman_harris", wf_blackman_harris, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "blackman_nuttall", wf_blackman_nuttall, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "flat_top", wf_flat_top, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "kbd", wf_kbd, -1);
	rb_define_module_function(rb_mWaveWindowFunction, "kaiser_bessel_derived", wf_kbd, -1);
}

/*******************************************************************************
	for C API
*******************************************************************************/

static inline void
wf_iter_make_rect(long N, double w[])
{
	for (volatile long n = 0; n < N; n++)
		w[n] = 1.;
}

static inline void
wf_iter_make_kurt(long N, double w[])
{
	if (N % 2 == 0)
	{
		for (volatile long n = 0; n < (N/2); n++)
		{
			if (n == 0)
			{
				w[n] = 0.;
				continue;
			}
			w[n] = 0.;
			w[N-n] = 0.;
		}
		w[N/2] = 1.;
	}
	else
	{
		for (volatile long n = 0; n < (N/2); n++)
		{
			w[n] = 0;
			w[N-1-n] = 0;
		}
		w[N/2] = 1.;
	}
}

static inline void
wf_iter_cb_sp(enum WFIF_SP_EVAL_TYPE handle, long N, double w[])
{
	switch (handle) {
	case WFIF_RECT:
		wf_iter_make_rect(N, w);
		break;
	case WFIF_KURT:
		wf_iter_make_kurt(N, w);
		break;
	case WFIF_NOCNTL:
	default:
		break;
	}
}

static inline enum WFIF_SP_EVAL_TYPE
wf_iter_errhdl(wf_iterfunc_t wfif)
{
	enum WFIF_SP_EVAL_TYPE handle = WFIF_NOCNTL;
	
	if (wfif.handle_param_nan != WFIF_NOCNTL && isnan(wfif.param))
		handle = wfif.handle_param_nan;
	else if (wfif.handle_param_inf != WFIF_NOCNTL && isinf(wfif.param))
		handle = wfif.handle_param_inf;
	else if (wfif.handle_param_zero != WFIF_NOCNTL && (wfif.param == 0))
		handle = wfif.handle_param_zero;
	
	return handle;
}

static inline void
wf_iter_rule_1d(wf_iterfunc_t wfif, long N, double w[])
{
	if (N % 2 == 0)
	{
		for (volatile long n = 0; n < (N/2); n++)
		{
			volatile const double value = wfif.iterfunc(n, N, wfif.param);
			if (n == 0)
			{
				w[n] = value;
				continue;
			}
			w[n] = value;
			w[N-n] = value;
		}
		w[N/2] = 1.;
	}
	else
	{
		for (volatile long n = 0; n < (N/2); n++)
		{
			volatile const double value = wfif.iterfunc(n+0.5, N, wfif.param);
			w[n] = value;
			w[N-1-n] = value;
		}
		w[N/2] = 1.;
	}
}

static inline void
wf_iter_rule_mdct(wf_iterfunc_t wfif, long N, double w[])
{
	double sum = 0.;
	
	if (N % 2 == 0)
	{
		for (volatile long n = 0; n < (N/2); n++)
		{
			volatile const double value = wfif.iterfunc(n, N, wfif.param);
			sum += value;
			w[n] = sum;
		}
		sum += wfif.iterfunc(N/2, N, wfif.param);
		for (volatile long n = 0; n < (N/2); n++)
		{
			RUBY_ASSERT(signbit(w[n]));
			w[n] = isinf(w[n]) ? 1. : sqrt(w[n]/sum);
			w[N-1-n] = w[n];
		}
	}
	else
	{
		for (volatile long n = 0; n < (N/2); n++)
		{
			volatile const double value = wfif.iterfunc(n+0.5, N, wfif.param);
			sum += value;
			w[n] = sum;
		}
		sum += wfif.iterfunc(N/2.0, N, wfif.param);
		for (volatile long n = 0; n < (N/2); n++)
		{
			RUBY_ASSERT(signbit(w[n]));
			w[n] = isinf(w[n]) ? 1. : sqrt(w[n]/sum);
			w[N-1-n] = w[n];
		}
		w[N/2] = 1.;
	}
}

void
wf_iter_cb(wf_iterfunc_t wfif, long N, double w[])
{
	RUBY_ASSERT(wfif.iterfunc == NULL);
	
	enum WFIF_SP_EVAL_TYPE handle = wf_iter_errhdl(wfif);
	
	if (handle != WFIF_NOCNTL)
	{
		wf_iter_cb_sp(handle, N, w);
	}
	else
	{
		switch (wfif.iter_rule) {
		case WFIF_ITER_1D:
			wf_iter_rule_1d(wfif, N, w);
			break;
		case WFIF_ITER_MDCT:
			wf_iter_rule_mdct(wfif, N, w);
			break;
		default:
			break;
		}
	}
}
