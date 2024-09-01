#ifndef RB_WAVE_GLOBALS_H_INCLUDED
#define RB_WAVE_GLOBALS_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include <ruby/internal/value.h> // VALUE
#include "ruby/ext_extern.h"

RUBY_EXT_EXTERN VALUE rb_mWave;
RUBY_EXT_EXTERN VALUE rb_cWavePCM;
RUBY_EXT_EXTERN VALUE rb_mWaveFFT;
RUBY_EXT_EXTERN VALUE rb_mWaveWindowFunction;
RUBY_EXT_EXTERN VALUE rb_cWaveRIFF;
RUBY_EXT_EXTERN VALUE rb_eWaveSemanticError;

#if defined(__cplusplus)
}
#endif

#endif /* RB_WAVE_GLOBALS_H_INCLUDED */
