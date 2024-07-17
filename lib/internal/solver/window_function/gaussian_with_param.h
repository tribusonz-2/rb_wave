#ifndef RB_WAVE_WF_GAUSSIAN_WITH_PARAM_H_INCLUDED
#define RB_WAVE_WF_GAUSSIAN_WITH_PARAM_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

/*
 * �p�����^�����ꂽ�K�E�X���̗��U�^�́C�ȉ��̂悤�ɒ�`�ł���D
 * 
 *   ```
 *   def gaussian(n, _N, sigma)
 *     t1 = _N % 2 == 0 ?
 *       -1 + 2.0 * n / _N :
 *       -1 + 2.0 * (n + 0.5) / _N
 *     t2 = 8 * sigma * sigma
 *     Math.exp(-(t1 * t1 / t2))
 *   end
 *   ```
 * 
 * �������Cn��N�ɂ��Ď��̊֌W
 * $0 \leq n \leq N$
 * �𖞂����D
 *
 * �����ŁC$\sigma$�͕W���΍��ł���C�p�����^�ɂ������D
 * $\sigma$�́A�v�Z�̂Ƃ�������킸������ɋ߂Â��ƁC���̖�����𕪕�Ƃ����W���ɂȂ�C���͕K��1�ɂȂ�D
 * �t�ɁC$\sigma$��0�ɋ߂��悤���ƁC�v�Z�̓[�����Z�ɂȂ蓚�̑��o�Ƃ܂ł͂����Ȃ��D����͋K�i�ł���IEEE754�{���x�ŏ���ۂĂȂ��ɂ��D
 * 
 *   ```
 *   def calc_stdev(sigma)
 *     8 * sigma * sigma
 *   end
 *   
 *   calc_stdev(1e-160) # => 8.0e-320
 *   calc_stdev(1e-170) # => 0.0
 *   ```
 * 
 * ���̂��߁C���U�^��x=0�ȕϐ��̂Ƃ��ɂ�$0.0/0.0$���v�Z���Ă��܂��CNaN{Not a Number}�𑗏o���Ă��܂��D
 * 
 *   ```
 *   len = 6
 *   len.times{|x| p gaussian(x, len, 0)}
 *   # => 0.0
 *   # => 0.0
 *   # => 0.0
 *   # => NaN
 *   # => 0.0
 *   # => 0.0
 *   ```
 * 
 * ����𗘗p���āC�W���΍���0�̂Ƃ��ɂ͈ȉ��̂悤�ɐ��v�Z�փX�C�b�`����D
 * $ w(x, 0) = \left\{ \begin{array}{cl} x \neq 0 & : 0 \\ x = 0 & : 1 \end{array} \right. $
 * 
 */
static inline double
wf_gaussian_calc_param(double sigma)
{
	return 8 * sigma * sigma;
}


static inline double
wf_gaussian_with_param_even(long n, long N, double t2)
{
	double t1 = -1 + 2. * n / N;
	return exp(-(t1 * t1 / t2));
}

static inline double
wf_gaussian_with_param_odd(long n, long N, double t2)
{
	double t1 = -1 + 2 * (n + 0.5) / N;
	return exp(-(t1 * t1 / t2));
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_GAUSSIAN_WITH_PARAM_H_INCLUDED */
