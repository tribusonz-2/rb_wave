#ifndef RB_WAVE_WF_KBD_WITH_PARAM_H_INCLUDED
#define RB_WAVE_WF_KBD_WITH_PARAM_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_kbd_with_param_eval(double n, long N, double alpha)
{
	const double x = n / N;
	double t1;
	
	t1 = 4.0 * x - 1.0;
	return cyl_bessel_i0(M_PI * alpha * sqrt(1.0 - t1 * t1));
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_KBD_WITH_PARAM_H_INCLUDED */
