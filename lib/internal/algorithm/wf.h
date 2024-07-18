#ifndef RB_WAVE_ALGO_WF_H_INCLUDED
#define RB_WAVE_ALGO_WF_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

// The design philosophy is based on Dr. Naofumi Aoki

typedef struct {
	double (*func)(double, long, double);
	double param;
} wf_iterfunc_t;

// 一元イテレータ
// 使用: #hann, #hamming
static inline void
wf_iter_cb(wf_iterfunc_t wfif, long N, double w[])
{
	if (N % 2 == 0) /* サイズが偶数のとき */
	{
		for (volatile long n = 0; n < (N/2); n++)
		{
			volatile const double value = wfif.func(n, N, wfif.param);
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
			volatile const double value = wfif.func(n+0.5, N, wfif.param);
			w[n] = value;
			w[N-1-n] = value;
		}
		w[N/2] = 1.;
	}
}

static inline void
wf_iter_paramzero(wf_iterfunc_t unused_obj, long N, double w[])
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

// 一元イテレータ (パラメータが0のとき特殊ルーチン設定済み仕様)
// 使用: #gaussian
static inline void
wf_iter_cb_with_paramzero(wf_iterfunc_t wfif, long N, double w[])
{
	if ((wfif.param != wfif.param) || wfif.param == 0)
		wf_iter_paramzero(wfif, N, w);
	else
		wf_iter_cb(wfif, N, w);
}


#include <ruby/internal/memory.h> // ALLOC_N()
#include <ruby/internal/intern/array.h> // rb_ary_new()

static inline VALUE
rb_wf_iter(void (*func)(double, long, double *), long len, double param)
{
	VALUE ary = rb_ary_new2(len);
	double *w = ALLOC_N(double, len);
	
	func(param, len, w);
	
	for (volatile long n = 0; n < len; n++)
		rb_ary_store(ary, n, DBL2NUM(w[n]));
	
	return ary;
}


#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_ALGO_WF_H_INCLUDED */
