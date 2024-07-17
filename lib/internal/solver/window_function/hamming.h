#ifndef RB_WAVE_WF_HAMMING_H_INCLUDED
#define RB_WAVE_WF_HAMMING_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_hamming_eval(double n, long N, double unused_param)
{
	static const double coef1 = 25./46.;
	static const double coef2 = 21./46.;
	
	return coef1 - coef2 * cos(2.0 * M_PI * n / N);
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_HAMMING_H_INCLUDED */
