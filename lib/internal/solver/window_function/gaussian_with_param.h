#ifndef RB_WAVE_WF_GAUSSIAN_WITH_PARAM_H_INCLUDED
#define RB_WAVE_WF_GAUSSIAN_WITH_PARAM_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

/*
 * パラメタ化されたガウス窓の離散型は，以下のように定義できる．
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
 * ただし，nはNについて次の関係
 * $0 \leq n \leq N$
 * を満たす．
 *
 * ここで，$\sigma$は標準偏差であり，パラメタにおかれる．
 * $\sigma$は、係数としての寄与のとき$\mathcal{R}\to\pm\infty$に対し1位の極を持つ．このときの$\frac{1}{8 \sigma^2}$の$\sigma\in \mathbb{R}$での値を1とする．
 * この計算は数値根を0とおく．したがって，$\sigma=0$のときは項をゼロ除算するものであり，電気数学としての解は
 * $w(n, N) = \left\{ \begin{array}{cl} n/N = \frac{1}{2} & : \ 1 \\ otherwise & : \ 0 \end{array} \right.$
 * の式を得る．
 * また$\sigma$がゼロに近いとき，数値計算ではゼロと等価である場合がある．規格であるIEEE754倍精度で情報を保てないのが原因である．
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
 * このため，${x}\neq{0}$を満たしているにも関わらず$0.0/0.0$を計算してしまい，NaN{Not a Number}を送出してしまう．
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
 */
static inline double
wf_gaussian_calc_param(double sigma)
{
	return 8 * sigma * sigma;
}

static inline double
wf_gaussian_with_param_eval(double n, long N, double t2)
{
	double t1 = -1 + 2. * n / N;
	return exp(-(t1 * t1 / t2));
}

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_WF_GAUSSIAN_WITH_PARAM_H_INCLUDED */
