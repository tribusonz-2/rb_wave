# Waveライブラリ for Ruby

　波形デジタルフィルタのRuby拡張ライブラリです。  
　スタンフォード大学向けの教材でもあり、アルゴリズム・コレクションでもあります。ケーススタディにどうぞ。  

#### 実装状況:
* `Wave::WindowFunction` (離散型窓関数)  
    * `#rectangular` (レクタンギュラ窓)  
    * `#hann` (ハン窓 / パラメタ化されたハン窓)  
    * `#hamming` (ハミング窓 / 一般化ハミング窓)  
    * `#bartlett` (バートレット窓)  
    * `#blackman` (ブラックマン窓)  
    * `#gaussian` (ガウス窓)  
    * `#kaiser` (カイザー窓)  
    * `#blackman_harris` (ブラックマン・ハリス窓)  
    * `#nuttall` (ナットール窓)  
    * `#blackman_nuttall` (ブラックマン・ナットール窓)  
    * `#flat_top` (フラットトップ窓)  
    * `#kbd` (カイザー・ベッセル派生窓)  
