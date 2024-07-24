/*******************************************************************************
	window_function.c -- Window Function
	
	Author: Hironobu Inatsuka
*******************************************************************************/
#include <ruby.h>
//#include <ruby/internal/memory.h> // ALLOC_N()
//#include <ruby/internal/intern/array.h> // rb_ary_new(), rb_ary_store()
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
 *  モジュールを使うのはFFTの窓掛けとして周波数特性の分析が最もで，デジタルフィルタを開発するならば，C APIを用いるなどして，コアクラス開発者との合流が望ましい．
 *  
 *  窓関数の配列は一度の生成につき100は下らない．マスター周波数が96kHzなら，その倍の配列数を生成することになる．
 *  このためコールバック方式を採用し，イテレーションに組み込むことで，高速な生成を実現している．
 *  設計思想は青木直史博士に基づいている．
 *  
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
	ディリクレ窓 / 矩形窓(レクタンギュラ窓)
*******************************************************************************/
#include "internal/solver/window_function/rectangular.h"

static void
wf_cb_rectangular(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_rectangular_eval,
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
 *  離散型ディリクレ窓の配列を返す．lenで配列数を指定する．
 *  ディリクレ窓は矩形窓(レクタンギュラ窓)としてよく知られ，よく使われる窓関数の一つである．常に1.0のスカラー量となる．
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
	ハミング窓 / 一般化ハミング窓
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

#include "internal/solver/window_function/generalized_hamming.h"

static void
wf_cb_generalized_hamming(double alpha, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_generalized_hamming_eval, 
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
 *  離散型ハミング窓の配列を返す．lenで配列数を指定する．
 *  ハミング窓はよく使われる窓関数の一つである．
 *  離散型のハミング窓は通常以下で定義される:
 *  $ w(x)=\frac{25}{46}-\frac{21}{46}\cos(2\pi{x}), 0 \leq x \leq 1 $
 *  
 *  第二引数に$\alpha$が当てられたときには配列数len分の離散型一般化ハミング窓の配列を返す．
 *  一般化ハミング窓は．ハン窓とハミング窓を一般化したものであり，通常は実数パラメタ$\alpha$を定義域$\frac{1}{2}\leq\alpha 1$とするものである．
 *  \alphaはこの定義域外の値のときにはRangeErrorの例外が発生する．電気数学では曖昧さも強みなため例外処理が働くのは珍しい．これは得られる値が期待とはかけ離れているためである．
 *  離散型の一般化ハミング窓は通常以下で定義される:
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
	ハン窓 / パラメタ化されたハン窓
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
 *    Wave::WindowFunction.hann(len, alpha) -> [*Float]
 *    Wave::WindowFunction.hanning(len, alpha) -> [*Float]
 *  
 *  離散型ハン窓の配列を返す．lenで配列数を指定する．
 *  ハン窓はよく使われる窓関数の一つである．
 *  この関数はパラメタ修正されたハミング窓に因んでハニング窓とも呼ばれる．
 *  離散型のハン窓は通常以下で定義される:
 *  $ w(x)=\frac{1}{2}-\frac{1}{2}\cos(2\pi{x}), 0 \leq x \leq 1 $
 *  ここで，係数$\alpha=\frac{1}{2}$は余弦の項の次数$1-\alpha$の関係がある．
 *  
 *  第二引数alphaが当てられた場合，パラメタ化されたハン窓をlen分の配列を確保し返す．
 *  $\alpha$は実数パラメタであり，定義域を$\frac{1}{2}\leq\alpha 1$とするものである．
 *  \alphaはこの定義域外の値のときにはRangeErrorの例外が発生する．電気数学では曖昧さも強みなため例外処理が働くのは珍しい．これは得られる値が期待とはかけ離れているためである． *  
 
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
	バートレット窓
*******************************************************************************/
#include "internal/solver/window_function/bartlett.h"

static void
wf_cb_bartlett(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_bartlett_eval,
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
 *  離散型バートレット窓の配列を返す．lenで配列数を指定する．
 *  バートレット窓は三角窓とも呼ばれ，よくリファレンスに出てくる窓関数である．
 *  定義式は以下で表せる．
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
	ブラックマン窓
*******************************************************************************/
#include "internal/solver/window_function/blackman.h"

static void
wf_cb_blackman(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = {
		wf_blackman_eval,
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
 *  離散型ブラックマン窓の配列を返す．lenで配列数を指定する．
 *  ブラックマン窓はよく使われる窓関数である．
 *  定義式は以下で表せる．
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
		return rb_wf_ary_new(wf_cb_gaussian, NUM2LONG(len), 0.);
	}
	else
	{
		return rb_wf_ary_new(wf_cb_gaussian_with_param, NUM2LONG(len), NUM2DBL(param));
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
		return rb_wf_ary_new(wf_cb_kaiser, NUM2LONG(len), 0.);
	}
	else
	{
		return rb_wf_ary_new(wf_cb_kaiser_with_param, NUM2LONG(len), NUM2DBL(param));
	}
}

/*******************************************************************************
	バートレット・ハン窓
*******************************************************************************/
#include "internal/solver/window_function/bartlett_hann.h"

static void
wf_cb_bartlett_hann(double unused_obj, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_bartlett_hann_eval, 
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
 *  離散型の修正バートレット・ハン窓の配列を返す．lenで配列数を指定する．
 *  修正バートレット・ハン窓は
 *  $ w(x)=0.62-0.48 |x - 0.5| + 0.38 \cos(2\pi(x - 0.5)) , 0 \leq x \leq 1 $
 *  で定義される．ここで$|x|$は絶対値である．
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
	ブラックマン・ハリス窓
*******************************************************************************/
#include "internal/solver/window_function/blackman_harris.h"

static void
wf_cb_blackman_harris(double unused_obj, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_blackman_harris_eval, 
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
 *  離散型ブラックマン・ハリス窓の配列を返す．lenで配列数を指定する．
 *  一般にブラックマン・ハリス窓は以下の式
 *  $ w(x)=a_0-a_1 cos(2\pi x) + a_2 cos(4\pi x) - a_3 cos(6\pi x), 0 \leq x \leq 1 $
 *  の最小4項で定義される．
 *  ここで$a_n$は最小4項の係数であり，以下の平均値・中央値が$\frac{1}{4}$な値を持つ．
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
	ナットール窓
*******************************************************************************/
#include "internal/solver/window_function/nuttall.h"

static void
wf_cb_nuttall(double unused_param, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_nuttall_eval, 
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
 *  離散型ナットール窓の配列を返す．lenで配列数を指定する．
 *  ナットール窓は以下の4項対称ブラックマン・ハリス窓
 *  $ w(x)=a_0-a_1 cos(2\pi x) + a_2 cos(4\pi x) - a_3 cos(6\pi x), 0 \leq x \leq 1 $
 *  をナットールの定義によるL点としたものである．
 *  ここで$a_n$は係数であり，4点を以下とする．
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
	ブラックマン・ナットール窓
*******************************************************************************/
#include "internal/solver/window_function/blackman_nuttall.h"

static void
wf_cb_blackman_nuttall(double unused_obj, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_blackman_nuttall_eval, 
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
 *  離散型ブラックマン・ナットール窓の配列を返す．lenで配列数を指定する．
 *  一般にブラックマン・ナットール窓は以下の式
 *  $ w(x)=a_0-a_1 cos(2\pi x) + a_2 cos(4\pi x) - a_3 cos(6\pi x), 0 \leq x \leq 1 $
 *  で定義される．
 *  ここで$a_n$は最小4項の係数であり，以下の値を持つ．
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
	フラットトップ窓
*******************************************************************************/
#include "internal/solver/window_function/flat_top.h"

static void
wf_cb_flat_top(double alpha, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_flat_top_eval, 
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
 *  離散型フラットトップ窓の配列を返す．lenで配列数を指定する．
 *  一般にフラットトップ窓は以下の式
 *  $w(x)=a_0-a_1\cos(2\pi x)+a_2\cos(4\pi x)-a_3\cos(6\pi x)+a_4\cos(8\pi x), 0\leq x\leq1$
 *  で定義される．
 *  ここで$a_n$は係数であり，値は各々以下である．
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
	KBD窓 (カイザー・ベッセル派生窓)
*******************************************************************************/
#include "internal/solver/window_function/kbd_with_param.h"

static void
wf_cb_kbd_with_param(double alpha, long len, double w[])
{
	wf_iterfunc_t wfif = { 
		wf_kbd_with_param_eval, 
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
 *  離散型KBD窓の配列を返す．lenで配列数を指定する．
 *  KBD窓はカイザー・ベッセル・派生窓の頭文字語であり，カイザー窓を修正離散コサイン変換(MDCT)での使用に設計したものである．
 *  
 *    Wave::WindowFunction.kbd(5, 3)
 *    # => [0.4114947429371883,
 *    # =>  0.9996957233074878,
 *    # =>  1.0,
 *    # =>  0.9996957233074878,
 *    # =>  0.4114947429371883]
 */
static VALUE
wf_kbd(int argc, VALUE *argv, VALUE unused_obj)
{
	VALUE len, param;
	rb_scan_args(argc, argv, "20", &len, &param);
	
	return rb_wf_ary_new(wf_cb_kbd_with_param, NUM2LONG(len), NUM2DBL(param));
}


/******************************************************************************/

// エントリポイント
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

// レクタンギュラの配列を作る
static inline void
wf_iter_make_rect(long N, double w[])
{
	for (volatile long n = 0; n < N; n++)
		w[n] = 1.;
}

// 中央値が1，他が0の配列を作る
static inline void
wf_iter_make_kurt(long N, double w[])
{
	if (N % 2 == 0) /* サイズが偶数のとき */
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
	else /* サイズが奇数のとき */
	{
		for (volatile long n = 0; n < (N/2); n++)
		{
			w[n] = 0;
			w[N-1-n] = 0;
		}
		w[N/2] = 1.;
	}
}

// 特殊な窓関数配列を生成する
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

// エラーハンドリング
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

// 一元イテレータ・ルール
static inline void
wf_iter_rule_1d(wf_iterfunc_t wfif, long N, double w[])
{
	if (N % 2 == 0) /* サイズが偶数のとき */
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
	else /* サイズが奇数のとき */
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

// MDCTイテレータ・ルール(畳み込み和ルーチン)
static inline void
wf_iter_rule_mdct(wf_iterfunc_t wfif, long N, double w[])
{
	double sum = 0.;
	
	if (N % 2 == 0) /* サイズが偶数のとき */
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
		w[N/2] = 1.;
	}
	else /* サイズが奇数のとき */
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
