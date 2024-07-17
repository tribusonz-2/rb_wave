#ifndef RB_WAVE_WF_HANN_H_INCLUDED
#define RB_WAVE_WF_HANN_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

// The design philosophy is based on Dr. Naofumi Aoki

static inline double
wf_hann_even(long n, long N, double unused_param)
{
	return 0.5 - 0.5 * cos(2.0 * M_PI * n / N);
}

static inline double
wf_hann_odd(long n, long N, double unused_param)
{
	return 0.5 - 0.5 * cos(2.0 * M_PI * (n + 0.5) / N);
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_HANN_H_INCLUDED */
