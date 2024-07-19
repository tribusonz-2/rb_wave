/*******************************************************************************
	window_function.c -- Window Function
	
	Author: Hironobu Inatsuka
*******************************************************************************/
#include <ruby.h>
#include "include/wave.h"
#include "internal/algorithm/wf.h"

#ifndef HAVE_CYL_BESSEL_I0
double cyl_bessel_i0(double);
# include "missing/cyl_bessel_i0.c"
#endif


/*
 *  module Wave::WindowFunction
 *  
 *  Wave::WindowFunctionモジュールは，波形フィルタリングでよく用いられる離散型の窓関数をRubyのユーザレベルから叩けるようにしたフロントエンドである．
 *  ユーザレベル実装は実行速度よりはアルゴリズム集としてのテストスィートの色が強い．デジタルフィルタリングではdouble型のスカラ型をアロケートして，役目を終えれば使い捨てるのが実際である．
 *  モジュールを使うのはFFTの窓掛けとして周波数特性の分析に使うのが最もで，デジタルフィルタを開発するならば，C APIを用いるなどして，コアクラス開発者との合流が望ましい．
 *  連続型の窓関数は0を中央値とするx軸を変数とし，定義域を$-\frac{1}{2}\leq{x}\leq\frac{1}{2}$とする．
 *  離散型はxの定義域を$0 \leq x \leq 1$と限定し，配列数を四分位とした離散信号を用いるという異なったアプローチを持っている．
 *  
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
 */


/*******************************************************************************
	ハン窓
*******************************************************************************/
#include "internal/solver/window_function/hann.h"

static void
wf_cb_hann(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_hann_eval,
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
 *  
 *  離散型ハン窓の配列を返す．lenで配列数を指定する．
 *  ハン窓はよく使われる窓関数の一つである．
 *  離散型のハン窓は通常以下で定義される:
 *  $ w(x)=\frac{1}{2}-\frac{1}{2}\cos(2\pi{x}), 0 \leq x \leq 1 $
 *  ここで，係数$\alpha=\frac{1}{2}$は余弦の項の次数$1-\alpha$におかれる．
 *  係数はパラメタ化することができ，その実効範囲は$\frac{1}{2}\leq\alpha\leq\1$である．
 *  ハン窓を改良したハミング窓の係数は$\frac{25}{46}$であり，範囲内である．
 *  
 *    Wave::WindowFunction.hann(5)
 *    # => [0.09549150281252627,
 *    # =>  0.6545084971874737,
 *    # =>  1.0,
 *    # =>  0.6545084971874737,
 *    # =>  0.09549150281252633]
 */
static VALUE
wf_hann(VALUE unused_obj, VALUE len)
{
	return rb_wf_iter(wf_cb_hann, NUM2LONG(len), 0.);
}


/*******************************************************************************
	ハミング窓
*******************************************************************************/
#include "internal/solver/window_function/hamming.h"

static void
wf_cb_hamming(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_hamming_eval, 
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
 *    Wave::WindowFunction.hamming(len) -> [*Float]
 *  
 *  離散型ハミング窓の配列を返す．lenで配列数を指定する．
 *  ハミング窓はよく使われる窓関数の一つである．
 *  離散型のハミング窓は通常以下で定義される:
 *  $ w(x)=\frac{25}{46}-\frac{21}{46}\cos(2\pi{x}), 0 \leq x \leq 1 $
 *  
 *     Wave::WindowFunction.hamming(5)
 *     # => [0.174144415611437,
 *     # =>  0.684551236562476,
 *     # =>  1.0,
 *     # =>  0.684551236562476,
 *     # =>  0.17414441561143706]
 */
static VALUE
wf_hamming(VALUE unused_obj, VALUE len)
{
	return rb_wf_iter(wf_cb_hamming, NUM2LONG(len), 0.);
}


/*******************************************************************************
	ガウス窓
*******************************************************************************/
#include "internal/solver/window_function/gaussian.h"

static void
wf_cb_gaussian(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_gaussian_eval, 
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
		wf_gaussian_with_param_eval, 
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
 
 *  離散型ガウス窓の配列を返す．lenで配列数を指定する．
 *  一般に，離散型のガウス窓は以下の式を満たす．
 *  $w(x)=e^{-((-1+2x)^2/8\sigma^2)}, 0 \leq{x}\leq 1$
 *  ここで，$\sigma$は標準偏差である．標準偏差は$3/10$としたとき，以下を等式とする．
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
		return rb_wf_iter(wf_cb_gaussian, NUM2LONG(len), 0.);
	}
	else
	{
		return rb_wf_iter(wf_cb_gaussian_with_param, NUM2LONG(len), NUM2DBL(param));
	}
}

/*******************************************************************************
	カイザー窓
*******************************************************************************/
#include "internal/solver/window_function/kaiser.h"

static void
wf_cb_kaiser(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_kaiser_eval, 
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
		wf_kaiser_with_param_eval, 
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
 *    Wave::WindowFunction.kaiser(x) -> [*Float]
 *    Wave::WindowFunction.kaiser(x, alpha) -> [*Float]
 *  
 *  離散型カイザー窓の配列を返す．lenで配列数を指定する．
 *  カイザー窓はカイザー・ベッセル窓としても知られ，一般に有限インパルス応答{FIR}フィルタ設計とスペクトル分析に使用される．
 *  離散型のカイザー窓は次式
 *  $ w(x)=\frac{I_0(\alpha 2 \sqrt{-(x-1)x}}{I_0(\alpha)} $
 *  を得る．
 *  ただし，$ I_n(x) $は第一種変形ベッセル関数，nはゼロ次であり．$\alpha$は形状パラメタである．
 *  以下は等価である．
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
		return rb_wf_iter(wf_cb_kaiser, NUM2LONG(len), 0.);
	}
	else
	{
		return rb_wf_iter(wf_cb_kaiser_with_param, NUM2LONG(len), NUM2DBL(param));
	}
}


#define TEST
#ifdef TEST
static VALUE
math_cyl_bessel_i0(VALUE unused_obj, VALUE x)
{
	return DBL2NUM(cyl_bessel_i0(NUM2DBL(x)));
}
#endif

/******************************************************************************/

void
InitVM_WindowFunction(void)
{
	rb_define_module_function(rb_mWaveWindowFunction, "hann", wf_hann, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "hanning", wf_hann, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "hamming", wf_hamming, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "gaussian", wf_gaussian, -1);
	rb_define_module_function(rb_mWaveWindowFunction, "kaiser", wf_kaiser, -1);

#ifdef TEST
	rb_define_module_function(rb_mMath, "cyl_bessel_i0", math_cyl_bessel_i0, 1);
#endif
}

