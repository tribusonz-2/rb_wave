#ifndef RB_WAVE_WF_FLAT_TOP_H_INCLUDED
#define RB_WAVE_WF_FLAT_TOP_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_flat_top_expr(double n, long N, double unused_param)
{
	static const double coef[5] = 
	{
		215578947 / 1000000000.0,
		416631580 / 1000000000.0,
		277263158 / 1000000000.0,
		 83578947 / 1000000000.0,
		  6947368 / 1000000000.0
	};
	
	return coef[0] - 
	       coef[1] * cos(2 * M_PI * n / N) +
	       coef[2] * cos(4 * M_PI * n / N) -
	       coef[3] * cos(6 * M_PI * n / N) +
	       coef[4] * cos(8 * M_PI * n / N);
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_FLAT_TOP_H_INCLUDED */
