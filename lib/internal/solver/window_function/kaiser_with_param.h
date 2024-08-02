#ifndef RB_WAVE_WF_KAISER_WITH_PARAM_H_INCLUDED
#define RB_WAVE_WF_KAISER_WITH_PARAM_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_kaiser_with_param_expr(double n, long N, double alpha)
{
	const double x = n / N;
	static double prev_alpha, denom;
	static bool check = false, reach_inf = false;
	
	if (!check || prev_alpha != alpha)
	{
		prev_alpha = alpha;
		denom = cyl_bessel_i0(prev_alpha);
		reach_inf = isinf(denom) ? true : false;
		if (!check)  check = true;
	}
	
	if (reach_inf)
		return x == 0.5;
	/* else */
		return cyl_bessel_i0(alpha * 2 * sqrt(-(x - 1) * x)) / denom;
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_KAISER_WITH_PARAM_H_INCLUDED */
