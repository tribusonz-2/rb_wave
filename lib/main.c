/*******************************************************************************
	main.c -- EntryPoint for Wave
	
	Author: Hironobu Inatsuka
*******************************************************************************/
#include <ruby.h>
#define USE_GLOBAL_VARIABLE
#include "ruby/wave/globals.h"

void InitVM_WindowFunction(void);

void
Init_wave(void)
{
	rb_mWave = rb_define_module("Wave");
	rb_cWavePCM = rb_define_class_under(rb_mWave, "PCM", rb_cObject);
	rb_mWaveFFT = rb_define_module_under(rb_mWave, "FFT");
	rb_mWaveWindowFunction = rb_define_module_under(rb_mWave, "WindowFunction");
	
	InitVM(WindowFunction);
}
