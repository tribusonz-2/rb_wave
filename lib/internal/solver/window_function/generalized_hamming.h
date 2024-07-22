#ifndef RB_WAVE_WF_GENERALIZED_HAMMING_H_INCLUDED
#define RB_WAVE_WF_GENERALIZED_HAMMING_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_generalized_hamming_calc_param(double alpha)
{
	if (0.5 <= alpha && alpha <= 1.0)
		return alpha;
	else
		rb_raise(rb_eRangeError, "parameter `alpha' is out of domain");
}

static inline double
wf_generalized_hamming_eval(double n, long N, double alpha)
{
	return alpha - (1 - alpha) * cos(2.0 * M_PI * n / N);
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_GENERALIZED_HAMMING_H_INCLUDED */
