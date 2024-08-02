#ifndef RB_WAVE_WF_NUTALL_H_INCLUDED
#define RB_WAVE_WF_NUTALL_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_nuttall_expr(double n, long N, double unused_param)
{
	static const double coef[4] = 
	{
		 88942 / 250000.,
		121849 / 250000.,
		 36058 / 250000.,
		  3151 / 250000.
	};
	
	return coef[0] -
	       coef[1] * cos(2 * M_PI * n / N) +
	       coef[2] * cos(4 * M_PI * n / N) -
	       coef[3] * cos(6 * M_PI * n / N);
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_NUTALL_H_INCLUDED */
