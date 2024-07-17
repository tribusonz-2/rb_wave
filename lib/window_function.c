/*******************************************************************************
	window_function.c -- Window Function
	
	Author: Hironobu Inatsuka
*******************************************************************************/
#include <ruby.h>
#include "include/wave.h"
#include "internal/algorithm/wf.h"

/*
 *  module Wave::WindowFunction
 *  
 *  Wave::WindowFunctionモジュールは，波形フィルタリングでよく用いられる離散型の窓関数をRubyのユーザレベルから叩けるようにしたフロントエンドである．
 *  ユーザレベル実装は実行速度よりはアルゴリズム集としてのテストスィートの色が強い．デジタルフィルタリングではdouble型のスカラ型をアロケートして，役目を終えれば使い捨てるのが実際である．
 *  モジュールを使うのはFFTの窓掛けとして周波数特性の分析に使うのが最もで，デジタルフィルタを開発するならば，C APIを用いるなどして，コアクラス開発者との合流が望ましい．
 *  連続型の窓関数は0を中央値とするx軸を変数とする．離散型はxのドメインを$0 \leq x \leq 1$と限定し，配列数を四分位とした離散信号を用いるという異なったアプローチを持っている．
 */

typedef struct {
	double (*func_even)(long, long, double);
	double (*func_odd)(long, long, double);
	double param;
} wf_iterfunc_t;

static inline VALUE
rb_wf_iter_cb(wf_iterfunc_t wfif, long len)
{
	VALUE ary = rb_ary_new2(len);
	
	if (len % 2 == 0) /* サイズが偶数のとき */
	{
		for (volatile long i = 0; i < len; i++)
		{
			volatile const double value = wfif.func_even(i, len, wfif.param);
			rb_ary_store(ary, i, DBL2NUM(value));
		}
	}
	else /* サイズが奇数のとき */
	{
		for (volatile long i = 0; i < len; i++)
		{
			volatile const double value = wfif.func_odd(i, len, wfif.param);
			rb_ary_store(ary, i, DBL2NUM(value));
		}
	}
	return ary;
}

static inline VALUE
rb_wf_iter_cb_with_zeroarg(wf_iterfunc_t wfif, long len)
{
	VALUE ary = rb_ary_new2(len);
	
	if (isnan(wfif.param) || wfif.param == 0)
	{
		if (len % 2 == 0) /* サイズが偶数のとき */
		{
			for (volatile long i = 0; i < len; i++)
			{
				volatile const double value = ((-1 + 2. * i / len) == 0);
				rb_ary_store(ary, i, DBL2NUM(value));
			}
		}
		else /* サイズが奇数のとき */
		{
			for (volatile long i = 0; i < len; i++)
			{
				volatile const double value = ((-1 + 2 * (i + 0.5) / len) == 0);
				rb_ary_store(ary, i, DBL2NUM(value));
			}
		}
	}
	else
	{
		if (len % 2 == 0) /* サイズが偶数のとき */
		{
			for (volatile long i = 0; i < len; i++)
			{
				volatile const double value = wfif.func_even(i, len, wfif.param);
				rb_ary_store(ary, i, DBL2NUM(value));
			}
		}
		else /* サイズが奇数のとき */
		{
			for (volatile long i = 0; i < len; i++)
			{
				volatile const double value = wfif.func_odd(i, len, wfif.param);
				rb_ary_store(ary, i, DBL2NUM(value));
			}
		}
	}
	return ary;
}


#include "internal/solver/window_function/hann.h"

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
	wf_iterfunc_t wfif = { wf_hann_even, wf_hann_odd, 0. };
	return rb_wf_iter_cb(wfif, NUM2LONG(len));
}


#include "internal/solver/window_function/hamming.h"

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
	wf_iterfunc_t wfif = { wf_hamming_even, wf_hamming_odd, 0. };
	return rb_wf_iter_cb(wfif, NUM2LONG(len));
}

#include "internal/solver/window_function/gaussian.h"
#include "internal/solver/window_function/gaussian_with_param.h"

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
		wf_iterfunc_t wfif = { wf_gaussian_even, wf_gaussian_odd, 0. };
		return rb_wf_iter_cb(wfif, NUM2LONG(len));
	}
	else
	{
		wf_iterfunc_t wfif = { 
			wf_gaussian_with_param_even, 
			wf_gaussian_with_param_odd, 
			wf_gaussian_calc_param(NUM2DBL(param)) 
		};
		return rb_wf_iter_cb_with_zeroarg(wfif, NUM2LONG(len));
	}
}


void
InitVM_WindowFunction(void)
{
	rb_define_module_function(rb_mWaveWindowFunction, "hann", wf_hann, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "hanning", wf_hann, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "hamming", wf_hamming, 1);
	rb_define_module_function(rb_mWaveWindowFunction, "gaussian", wf_gaussian, -1);
}

