#ifndef RB_WAVE_ALGO_WF_H_INCLUDED
#define RB_WAVE_ALGO_WF_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

enum WFIF_ITER_RULE {
	WFIF_ITER_1D,  // 1-dimensional Iterator
	WFIF_ITER_MDCT // Modified-DCT with convolutional summation
} ;

enum WFIF_SP_EVAL_TYPE {
	WFIF_NOCNTL,  // No control
	WFIF_RECT,    // Rectangular
	WFIF_KURT,    // if center: 1.0, otherwise: 0.0
} ;


typedef struct {
	double (*iterfunc)(double, long, double);
	double param;
	enum WFIF_ITER_RULE iter_rule;
	enum WFIF_SP_EVAL_TYPE handle_param_nan;
	enum WFIF_SP_EVAL_TYPE handle_param_inf;
	enum WFIF_SP_EVAL_TYPE handle_param_zero;
} wf_iterfunc_t;


void wf_iter_cb(wf_iterfunc_t wfif, long N, double w[]);

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
