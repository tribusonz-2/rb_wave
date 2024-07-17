#ifndef RB_WAVE_ALGO_WF_H_INCLUDED
#define RB_WAVE_ALGO_WF_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

// The design philosophy is based on Dr. Naofumi Aoki

static inline void
wf_iterfunc(
  double (*func_even)(long, long, double),
  double (*func_odd)(long, long, double),
  long N, double w[])
{
	if (N % 2 == 0) /* サイズが偶数のとき */
	{
		for (volatile long n = 0; n < N; n++)
		{
			w[n] = func_even(n, N, 0.);
		}
	}
	else /* サイズが奇数のとき */
	{
		for (volatile long n = 0; n < N; n++)
		{
			w[n] = func_odd(n, N, 0.);
		}
	}
}

#include <ruby.h> // VALUE, rb_ary_new(), rb_ary_store()




#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_ALGO_WF_H_INCLUDED */
