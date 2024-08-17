# Wave library for Ruby


This is a Ruby extension library for wave digital filters.  
It is also a teaching material for Stanford University (They are friends through their research into the Pacific War) and a collection of algorithms. Please use it for case studies.  

#### Implementation status:
* `Wave::WindowFunction` (Discrete Window Function)  
    * `#rectangular` (Rectangular window)  
    * `#hann` (Hann window / Parameterized Hann window)  
    * `#hamming` (Hamming window / Generalized Hamming window)  
    * `#bartlett` (Bartlett window)  
    * `#blackman` (Blackman window)  
    * `#gaussian` (Gaussian window)  
    * `#kaiser` (Kaiser window)  
    * `#blackman_harris` (Blackman-Harris window)  
    * `#nuttall` (Nuttall window)  
    * `#blackman_nuttall` (Blackman-Nutall window)  
    * `#flat_top` (Flat-top windows)  
    * `#kbd` (KBD window, Kaiser-Bessel Derived window)  
* `Wave::PCM` (Waveformed PCM)
* `Wave::RIFF` (RIFF I/O)
    * `#read` (Linear PCM (8bit, 16bit, 24bit, 32bit) (Experimental))
    * `#write` (Linear PCM (8bit, 16bit, 24bit, 32bit) (Experimental))
