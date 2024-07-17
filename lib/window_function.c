/*******************************************************************************
	window_function.c -- Window Function
	
	Author: Hironobu Inatsuka
*******************************************************************************/
#include <ruby.h>
#include "include/wave.h"
#include "internal/algorithm/wf.h"

typedef struct {
	double (*func_even)(long, long, double);
	double (*func_odd)(long, long, double);
	double param;
} wf_iterfunc_t;

static inline VALUE
rb_wf_iter_cb(wf_iterfunc_t wfif, long len)
{
	VALUE ary = rb_ary_new2(len);
	
	if (len % 2 == 0) /* �T�C�Y�������̂Ƃ� */
	{
		for (volatile long i = 0; i < len; i++)
		{
			volatile const double value = wfif.func_even(i, len, wfif.param);
			rb_ary_store(ary, i, DBL2NUM(value));
		}
	}
	else /* �T�C�Y����̂Ƃ� */
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
		if (len % 2 == 0) /* �T�C�Y�������̂Ƃ� */
		{
			for (volatile long i = 0; i < len; i++)
			{
				volatile const double value = ((-1 + 2. * i / len) == 0);
				rb_ary_store(ary, i, DBL2NUM(value));
			}
		}
		else /* �T�C�Y����̂Ƃ� */
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
		if (len % 2 == 0) /* �T�C�Y�������̂Ƃ� */
		{
			for (volatile long i = 0; i < len; i++)
			{
				volatile const double value = wfif.func_even(i, len, wfif.param);
				rb_ary_store(ary, i, DBL2NUM(value));
			}
		}
		else /* �T�C�Y����̂Ƃ� */
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
 *  ���U�^�n�����̔z���Ԃ��Dlen�Ŕz�񐔂��w�肷��D
 *  �n�����͂悭�g���鑋�֐��̈�ł���D
 *  ���U�^�̃n�����͒ʏ�ȉ��Œ�`�����:
 *  $ w(x)=\frac{1}{2}-\frac{1}{2}\cos(2\pi{x}), 0 \leq x \leq 1 $
 *  �����ŁC�W��$\alpha=\frac{1}{2}$�͗]���̍��̎���$1-\alpha$�ɂ������D
 *  �W���̓p�����^�����邱�Ƃ��ł��C���̎����͈͂�$\frac{1}{2}\leq\alpha\leq\1$�ł���D
 *  �n���������ǂ����n�~���O���̌W����$\frac{25}{46}$�ł���C�͈͓��ł���D
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
 *  ���U�^�n�~���O���̔z���Ԃ��Dlen�Ŕz�񐔂��w�肷��D
 *  �n�~���O���͂悭�g���鑋�֐��̈�ł���D
 *  ���U�^�̃n�~���O���͒ʏ�ȉ��Œ�`�����:
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
 
 *  ���U�^�K�E�X���̔z���Ԃ��Dlen�Ŕz�񐔂��w�肷��D
 *  ��ʂɁC���U�^�̃K�E�X���͈ȉ��̎��𖞂����D
 *  $w(x)=e^{-((-1+2x)^2/8\sigma^2)}$
 *  �����ŁC$\sigma$�͕W���΍��ł���D�W���΍���$3/10$�Ƃ����Ƃ��C�ȉ��𓙎��Ƃ���D
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

