#ifndef RB_WAVE_ALGO_WF_H_INCLUDED
#define RB_WAVE_ALGO_WF_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

// The design philosophy is based on Dr. Naofumi Aoki

#include <ruby/assert.h> // RUBY_ASSERT

enum WFIF_ITER_TYPE {
	WFIF_ITER_1D,   // 1-dimensional Iterator
	WFIF_ITER_MDCT1 // Modified-DCT with convolusional product
} ;

enum WFIF_SP_EVAL_TYPE {
	WFIF_NOCNTL,  // No control
	WFIF_RECT,    // Rectangular
	WFIF_KURT,    // if center: 1.0, otherwise: 0.0
} ;


typedef struct {
	double (*iterfunc)(double, long, double);
	double param;
	enum WFIF_ITER_TYPE handle_iter;
	enum WFIF_SP_EVAL_TYPE handle_param_nan;
	enum WFIF_SP_EVAL_TYPE handle_param_inf;
	enum WFIF_SP_EVAL_TYPE handle_param_zero;
} wf_iterfunc_t;


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

// 一元イテレータ
static inline void
wf_iter_cb_1d(wf_iterfunc_t wfif, long N, double w[])
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


// 一元イテレータ
static inline void
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
		switch (wfif.handle_iter) {
		case WFIF_ITER_1D:
			wf_iter_cb_1d(wfif, N, w);
			break;
		case WFIF_ITER_MDCT1:
		default:
			break;
		}
	}
}


#include <ruby/internal/memory.h> // ALLOC_N()
#include <ruby/internal/intern/array.h> // rb_ary_new(), rb_ary_store()

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
