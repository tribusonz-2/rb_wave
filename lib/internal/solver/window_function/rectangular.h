#ifndef RB_WAVE_WF_RECTANGULAR_H_INCLUDED
#define RB_WAVE_WF_RECTANGULAR_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

static inline double
wf_rectangular_eval(double n, long N, double unused_param)
{
	return 1.;
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_RECTANGULAR_H_INCLUDED */
