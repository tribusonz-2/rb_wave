#ifndef RB_WAVE_WF_BARTLETT_H_INCLUDED
#define RB_WAVE_WF_BARTLETT_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_bartlett_expr(double n, long N, double unused_param)
{
	return 1 - 2 * fabs(n / N - 0.5);
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_BARTLETT_H_INCLUDED */
