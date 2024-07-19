#ifndef RB_WAVE_WF_KAISER_H_INCLUDED
#define RB_WAVE_WF_KAISER_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_kaiser_eval(double n, long N, double unused_param)
{
	static double i0_3 = 0.;
	const double x = n / N;
	
	if (i0_3 == 0)
		i0_3 = cyl_bessel_i0(3);
	
	return cyl_bessel_i0(6 * sqrt(-(x - 1) * x)) / i0_3;
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_KAISER_H_INCLUDED */
