#ifndef RB_WAVE_WF_BLACKMAN_NUTTALL_H_INCLUDED
#define RB_WAVE_WF_BLACKMAN_NUTTALL_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_blackman_nuttall_expr(double n, long N, double unused_param)
{
	static const double coef[4] = 
	{
		3635819 / 10000000.,
		4891775 / 10000000.,
		1365995 / 10000000.,
		 106411 / 10000000.
	};
	
	return coef[0] - 
	       coef[1] * cos(2 * M_PI * n / N) +
	       coef[2] * cos(4 * M_PI * n / N) -
	       coef[3] * cos(6 * M_PI * n / N);
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_BLACKMAN_NUTTALL_H_INCLUDED */
