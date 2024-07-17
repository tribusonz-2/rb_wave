#ifndef RB_WAVE_WF_GAUSSIAN_H_INCLUDED
#define RB_WAVE_WF_GAUSSIAN_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_gaussian_eval(double n, long N, double unused_param)
{
	static const double t1 = -(25./18.);
	const double t2 = (-1 + 2. * n / N);
	return exp(t1 * t2 * t2);
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_GAUSSIAN_H_INCLUDED */
