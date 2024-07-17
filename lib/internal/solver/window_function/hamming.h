#ifndef RB_WAVE_WF_HAMMING_H_INCLUDED
#define RB_WAVE_WF_HAMMING_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

// The design philosophy is based on Dr. Naofumi Aoki

static inline double
wf_hamming_even(long n, long N, double unused_param)
{
	return 25./46. - 21./46. * cos(2.0 * M_PI * n / N);
}

static inline double
wf_hamming_odd(long n, long N, double unused_param)
{
	return 25./46. - 21./46. * cos(2.0 * M_PI * (n + 0.5) / N);
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_HAMMING_H_INCLUDED */
