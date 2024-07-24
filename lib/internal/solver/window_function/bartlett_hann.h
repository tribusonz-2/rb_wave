#ifndef RB_WAVE_WF_BARTLETT_HANN_H_INCLUDED
#define RB_WAVE_WF_BARTLETT_HANN_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_bartlett_hann_eval(double n, long N, double unused_param)
{
	return 0.62 - 0.48 * fabs(n / N - 0.5) + 
	       0.38 * cos(2 * M_PI * (n / N - 0.5));
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_BARTLETT_HANN_H_INCLUDED */
