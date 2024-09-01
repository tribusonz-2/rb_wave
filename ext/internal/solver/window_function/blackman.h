#ifndef RB_WAVE_WF_BLACKMAN_H_INCLUDED
#define RB_WAVE_WF_BLACKMAN_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_blackman_expr(double n, long N, double unused_param)
{
	return 0.42 - 
		0.5 * cos(2 * M_PI * n / N) + 
		0.08 * cos(4 * M_PI * n / N);
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_BLACKMAN_H_INCLUDED */
