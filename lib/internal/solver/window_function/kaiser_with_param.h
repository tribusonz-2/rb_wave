#ifndef RB_WAVE_WF_KAISER_WITH_PARAM_H_INCLUDED
#define RB_WAVE_WF_KAISER_WITH_PARAM_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_kaiser_calc_param(double alpha)
{
	return cyl_bessel_i0(alpha);
}

static inline double
wf_kaiser_with_param_eval(double n, long N, double alpha)
{
	const double x = n / N, denom = cyl_bessel_i0(alpha);
	
	if (isinf(denom))
		return x == 0.5;
	/* else */
		return cyl_bessel_i0(alpha * 2 * sqrt(-(x - 1) * x)) / denom;
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_KAISER_WITH_PARAM_H_INCLUDED */
