#ifndef RB_WAVE_WF_HANN_H_INCLUDED
#define RB_WAVE_WF_HANN_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_hann_expr(double n, long N, double unused_param)
{
	return 0.5 - 0.5 * cos(2.0 * M_PI * n / N);
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_HANN_H_INCLUDED */
