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
 * $\sigma$は、計算のとき正負問わず無限大に近づくと，正の無限大を分母とした係数になり，解は必ず1になる．
 * 逆に，$\sigma$が0に近いようだと，計算はゼロ除算になり答の送出とまではいかない．これは規格であるIEEE754倍精度で情報を保てないによる．
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
 * このため，離散型でx=0な変数のときには$0.0/0.0$を計算してしまい，NaN{Not a Number}を送出してしまう．
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
 * これを利用して，標準偏差が0のときには以下のように専門計算へスイッチする．
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
