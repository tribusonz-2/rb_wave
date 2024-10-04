/*******************************************************************************
	riff.c - 
	
	$author$
	
	@license: MIT Licence
	
*******************************************************************************/
#include <ruby.h>
#include <ruby/io.h>
#include "ruby/wave/globals.h"
#include "ruby/wave/pcm.h"
#include "internal/riffchunk.h"
#include <stdint.h>

#define SupportedVersion "1.0.0"

struct RIFF {
	VALUE f; // File
	uint16_t format_type; // wave_format_type
	uint16_t bits; // bits_per_sample (Dependence by 'wave_format_type')
	uint32_t sampling_rate; // blocks_per_second
	uint16_t channels; // Number of channels
	uint16_t samples_per_block;
} ;



static uint16_t
dbyte2u16le(VALUE bytes)
{
	unsigned char *ptr = (unsigned char *)StringValuePtr(bytes);
	uint16_t size = ptr[0] | ptr[1] << 8;
	return size;
}

static uint32_t
qbyte2u32le(VALUE bytes)
{
	unsigned char *ptr = (unsigned char *)StringValuePtr(bytes);
	uint32_t size = ptr[0] | ptr[1] << 8 | ptr[2] << 16 | ptr[3] << 24;
	return size;
}


void
pcm_read_8bit(unsigned char buf[], double s[])
{
	char data = (buf[0] - 0x80);
	*s = data / (double)0x80;
}


void
pcm_read_16bit(unsigned char buf[], double s[])
{
	int16_t data = (int16_t)(buf[0] | buf[1] << 8);
	*s = data / (double)0x8000;
}

void
pcm_read_24bit(unsigned char buf[], double s[])
{
	 int data = (buf[0] | buf[1] << 8 | buf[2] << 16);
	 if (data & 0x800000)  data -= 0x1000000;
	*s = data / (double)0x800000;
}

void
pcm_read_32bit(unsigned char buf[], double s[])
{
	int32_t data = (int32_t)(buf[0] | buf[1] << 8 | buf[2] << 16 | buf[3] << 24);
	*s = data / (double)0x80000000;
}

static inline void
must_be_nonzero_error(const char *memb)
{
	rb_raise(rb_eWaveSemanticError, "'%s' must be non-zero", memb);
}

static void
io_readpartial(VALUE io, VALUE io_buf, long len)
{
	static ID readpartial;
	if (!readpartial)
		readpartial = rb_intern_const("readpartial");
	rb_funcall(io, readpartial, 2, LONG2FIX(len), io_buf);
}


static inline VALUE
wave_read_linear_pcm(char *file_name)
{
	const int BUFFER_SIZE = 0x1000;
	VALUE io = rb_file_open(file_name, "rb");
	VALUE io_buf = rb_str_new(0,0);
	
	uint32_t riff_chunk_size;
	uint32_t fmt_chunk_size;
	uint16_t wave_format_type;
	uint16_t channels;
	uint32_t samples_per_sec;
	uint32_t bytes_per_sec;
	uint16_t block_size;
	uint16_t bits_per_sample;
	uint32_t data_chunk_size;
	
	VALUE pcm_ary;
	double **mat;
	long length;
	void (*func)(unsigned char *, double *);
	long idx;
	int buffer_size;
	uint16_t samples_per_block = 1;
	
	// RIFF chunk
	io_readpartial(io, io_buf, 4);
	if (!RTEST(rb_str_equal(io_buf, rb_str_new_cstr("RIFF"))))
		rb_raise(rb_eWaveSemanticError, "unknown RIFF chunk ID: %"PRIsVALUE"", io_buf);
	
	io_readpartial(io, io_buf, 4);
	riff_chunk_size = qbyte2u32le(io_buf);
	
	io_readpartial(io, io_buf, 4);
	if (!RTEST(rb_str_equal(io_buf, rb_str_new_cstr("WAVE"))))
		rb_raise(rb_eWaveSemanticError, "unknown file format type: %"PRIsVALUE"", io_buf);

	// format chunk
	io_readpartial(io, io_buf, 4);
	if (!RTEST(rb_str_equal(io_buf, rb_str_new_cstr("fmt "))))
		rb_raise(rb_eWaveSemanticError, "no format chunk");
	
	io_readpartial(io, io_buf, 4);
	fmt_chunk_size = qbyte2u32le(io_buf);
	
	io_readpartial(io, io_buf, 2);

	wave_format_type = dbyte2u16le(io_buf);
	if (wave_format_type != 1)
		rb_raise(rb_eWaveSemanticError, "not a linear PCM");
	
	io_readpartial(io, io_buf, 2);
	channels = dbyte2u16le(io_buf);
	if (!channels)
		must_be_nonzero_error("channels");
	
	io_readpartial(io, io_buf, 4);
	samples_per_sec = qbyte2u32le(io_buf);
	if (!samples_per_sec)
		must_be_nonzero_error("samples_per_sec");
	
	io_readpartial(io, io_buf, 4);
	bytes_per_sec = qbyte2u32le(io_buf);
	if (!bytes_per_sec)
		must_be_nonzero_error("bytes_per_sec");
	
	io_readpartial(io, io_buf, 2);
	block_size = dbyte2u16le(io_buf);
	if (!block_size)
		must_be_nonzero_error("block_size");
	
	io_readpartial(io, io_buf, 2);
	bits_per_sample = dbyte2u16le(io_buf);
	if (!bits_per_sample)
		must_be_nonzero_error("bits_per_sample");
	
	if ((bits_per_sample / 8 * channels) != block_size)
		rb_raise(rb_eWaveSemanticError, "'block_size' mismatch");
	
	if ((samples_per_sec * block_size) != bytes_per_sec)
		rb_raise(rb_eWaveSemanticError, "'bytes_per_sec' mismatch");
	
	// data chunk
	io_readpartial(io, io_buf, 4);
	if (!RTEST(rb_str_equal(io_buf, rb_str_new_cstr("data"))))
		rb_raise(rb_eWaveSemanticError, "no data chunk");
	
	io_readpartial(io, io_buf, 4);
	data_chunk_size = qbyte2u32le(io_buf);
	
	if ((data_chunk_size % block_size) != 0)
		rb_raise(rb_eWaveSemanticError, "'data_chunk_size' is not a multiple of 'block_size'");
	
	switch (bits_per_sample) {
	case 8:  func = pcm_read_8bit;  break;
	case 16: func = pcm_read_16bit; break;
	case 24: func = pcm_read_24bit; break;
	case 32: func = pcm_read_32bit; break;
	default: rb_raise(rb_eWaveSemanticError, 
		"unrecognized (or unsupported) bits per sample: %d (for wave format type: %d)", 
		bits_per_sample, wave_format_type);
		break;
	}
	
	length = data_chunk_size / block_size;
	pcm_ary = rb_ary_new2(channels);
	for (long i = 0; i < channels; i++)
	{
		rb_ary_store(pcm_ary, i, rb_pcm_new(length, samples_per_sec));
	}
	
	mat = ALLOCA_N(double*, channels);
	for (long i = 0; i < channels; i++)
	{
		VALUE obj = rb_ary_entry(pcm_ary, i);
		mat[i] = WaveformDataPtr(obj);
	}
	
	buffer_size = BUFFER_SIZE / block_size * block_size;
	
	idx = 0;
	for (long data_offset = 0; data_offset < data_chunk_size; data_offset += buffer_size)
	{
		if ((1. * data_offset + buffer_size) > data_chunk_size)
			buffer_size = data_chunk_size - data_offset;
		io_readpartial(io, io_buf, buffer_size);
		unsigned char *buf_ptr = (unsigned char *)RSTRING_PTR(io_buf);
		for ( ; ; )
		{
			for (long i = 0; i < channels; i++)
			{
				double *s_ptr = mat[i];
				func(buf_ptr+(i*block_size/channels), s_ptr+idx);
			}
			idx += samples_per_block;
			buf_ptr += block_size;
			
			if (buf_ptr == (unsigned char *)RSTRING_END(io_buf))
				break;
		}
	}
	rb_str_resize(io_buf, 0);
	rb_io_close(io);
	return pcm_ary;
}

static VALUE
test_wave_read_linear_pcm(VALUE unused_obj, VALUE fname)
{
	return wave_read_linear_pcm(StringValuePtr(fname));
}

static bool
ary_all_pcm_p(VALUE ary)
{
	if (TYPE(ary) != T_ARRAY)
		rb_raise(rb_eTypeError, "not an %"PRIsVALUE" includes %"PRIsVALUE"", rb_cArray, rb_cWavePCM);
	for (long i = 0; i < RARRAY_LEN(ary); i++)
	{
		VALUE elem = rb_ary_entry(ary, i);
		if (CLASS_OF(elem) != rb_cWavePCM)
			return false;
	}
	return true;
}


static inline double
fclip(double x, double min, double max)
{
	return fmin(fmax(min, x), max);
}

static inline double
wave_normalize(double x, double min, double max, double rate)
{
	if (x != x)  x = 0;
	return fclip(x * rate, min, max);
}

void
pcm_write_8bit(unsigned char buf[], double s[])
{
	double digitize = wave_normalize(s[0], INT8_MIN, INT8_MAX, 0x80);
	
	*buf = (unsigned char)(digitize + 0x80);
}

void
pcm_write_16bit(unsigned char buf[], double s[])
{
	double digitize = wave_normalize(s[0], INT16_MIN, INT16_MAX, 0x8000);
	int16_t bytes = (int16_t)digitize;
	
	buf[0] = bytes & 0xFF;
	buf[1] = (bytes >> 8) & 0xFF;
}

void
pcm_write_24bit(unsigned char buf[], double s[])
{
	double digitize = wave_normalize(s[0], -0x800000, 0x7FFFFF, 0x800000);
	int32_t bytes = (int32_t)digitize;

	buf[0] = bytes & 0xFF;
	buf[1] = (bytes >> 8) & 0xFF;
	buf[2] = (bytes >> 16) & 0xFF;
}

void
pcm_write_32bit(unsigned char buf[], double s[])
{
	double digitize = wave_normalize(s[0], INT32_MIN, INT32_MAX, 0x80000000);
	int32_t bytes = (int32_t)digitize;

	buf[0] = bytes & 0xFF;
	buf[1] = (bytes >> 8) & 0xFF;
	buf[2] = (bytes >> 16) & 0xFF;
	buf[3] = (bytes >> 24) & 0xFF;
}

static void
io_writepartial(VALUE io, VALUE buf)
{
	if (rb_io_bufwrite(io, (unsigned char *)StringValuePtr(buf), RSTRING_LEN(buf)) == -1)
		rb_raise(rb_eIOError, "write failure");
}

static VALUE
rb_str_cat_uintle(VALUE str, uint32_t value, size_t sz)
{
	static char s[4];
	rb_integer_pack(ULL2NUM(value), s, sz, 8, 0, 
		INTEGER_PACK_LITTLE_ENDIAN | INTEGER_PACK_2COMP);
	return rb_str_buf_cat(str, s, sz);
}

static VALUE
rb_str_buf_z_new(long len)
{
	VALUE bin = rb_str_buf_new(len);
	char *ptr = RSTRING_PTR(bin);
	rb_str_resize(bin, len);
	MEMZERO(ptr, char, len);
	
	return bin;
}

static VALUE
rb_str_buf_z_resize(VALUE bin, long len)
{
	char *ptr = RSTRING_PTR(bin);
	if (RSTRING_LEN(bin) == len)
	{
		MEMZERO(ptr, char, len);
	}
	else if (RSTRING_LEN(bin) < len)
	{
		rb_str_resize(bin, len);
		MEMZERO(ptr, char, len);
	}
	else if (RSTRING_LEN(bin) > len)
	{
		rb_str_resize(bin, len);
	}

	return bin;
}


static inline VALUE
wave_write_linear_pcm(VALUE pcm_ary, int16_t bits, char *file_name)
{
	const int BUFFER_SIZE = 0x1000;
	
	if (!ary_all_pcm_p(pcm_ary))
		rb_raise(rb_eArgError, "not a %"PRIsVALUE"", rb_cWavePCM);
	
	VALUE io = rb_file_open(file_name, "wb");
	VALUE io_buf;
	
	char riff_chunk_ID[4];
	uint32_t riff_chunk_size;
	char file_format_type[4];
	char fmt_chunk_ID[4];
	uint32_t fmt_chunk_size;
	uint16_t wave_format_type;
	uint16_t channels;
	uint32_t samples_per_sec;
	uint32_t bytes_per_sec;
	uint16_t block_size;
	uint16_t bits_per_sample;
	char data_chunk_ID[4];
	uint32_t data_chunk_size;
	
	double **mat;
	long length;
	void (*func)(unsigned char *, double *);
	long idx;
	int buffer_size;
	uint16_t samples_per_block = 1;
	
	riff_chunk_ID[0] = 'R';
	riff_chunk_ID[1] = 'I';
	riff_chunk_ID[2] = 'F';
	riff_chunk_ID[3] = 'F';

	file_format_type[0] = 'W';
	file_format_type[1] = 'A';
	file_format_type[2] = 'V';
	file_format_type[3] = 'E';
	
	fmt_chunk_ID[0] = 'f';
	fmt_chunk_ID[1] = 'm';
	fmt_chunk_ID[2] = 't';
	fmt_chunk_ID[3] = ' ';
	fmt_chunk_size = 16;
	wave_format_type = 1;

	if (RARRAY_LEN(pcm_ary) > UINT16_MAX)
		rb_raise(rb_eRangeError, "too many PCM classes");
	channels = (uint16_t)RARRAY_LEN(pcm_ary);

	samples_per_sec = 0;
	length = 0;
	mat = ALLOCA_N(double*, channels);
	for (long i = 0; i < channels; i++)
	{
		VALUE obj = rb_ary_entry(pcm_ary, i);
		mat[i] = WaveformDataPtr(obj);
		if (!samples_per_sec)
			samples_per_sec = rb_pcm_fs(obj);
		else
			if (samples_per_sec != rb_pcm_fs(obj))
				rb_raise(rb_eRuntimeError, 
				"Exporting each channel's the different sampling frequency is not supported yet");
		if (!length)
			length = rb_pcm_len(obj);
		else
			if (length != rb_pcm_len(obj))
				rb_raise(rb_eRuntimeError, 
				"Exporting each channel's the different length is not supported yet");
	}
	

	bits_per_sample = bits;
	switch (bits_per_sample) {
	case 8:  func = pcm_write_8bit;  break;
	case 16: func = pcm_write_16bit; break;
	case 24: func = pcm_write_24bit; break;
	case 32: func = pcm_write_32bit; break;
	default: rb_raise(rb_eWaveSemanticError, 
		"unrecognized (or unsupported) bits per sample: %d (for wave format type: %d)", 
		bits_per_sample, wave_format_type);
		break;
	}
	
	bytes_per_sec = samples_per_sec * bits_per_sample / 8 * channels;
	block_size = bits_per_sample / 8 * channels;
	
	data_chunk_ID[0] = 'd';
	data_chunk_ID[1] = 'a';
	data_chunk_ID[2] = 't';
	data_chunk_ID[3] = 'a';
	data_chunk_size = length * bits_per_sample / 8 * channels;
	
	riff_chunk_size = 36 + data_chunk_size;
	if (riff_chunk_size % 2 == 1)  riff_chunk_size++;
	
	io_buf = rb_str_new(0, 0);
	rb_str_buf_cat(io_buf, riff_chunk_ID, 4);
	rb_str_cat_uintle(io_buf, riff_chunk_size, 4);
	rb_str_buf_cat(io_buf, file_format_type, 4);
	io_writepartial(io, io_buf);
	
	io_buf = rb_str_new(0, 0);
	rb_str_buf_cat(io_buf, fmt_chunk_ID, 4);
	rb_str_cat_uintle(io_buf, fmt_chunk_size, 4);
	rb_str_cat_uintle(io_buf, wave_format_type, 2);
	rb_str_cat_uintle(io_buf, channels, 2);
	rb_str_cat_uintle(io_buf, samples_per_sec, 4);
	rb_str_cat_uintle(io_buf, bytes_per_sec, 4);
	rb_str_cat_uintle(io_buf, block_size, 2);
	rb_str_cat_uintle(io_buf, bits_per_sample, 2);
	io_writepartial(io, io_buf);
	
	io_buf = rb_str_new(0, 0);
	rb_str_buf_cat(io_buf, data_chunk_ID, 4);
	if (data_chunk_size % 2 == 1)
		rb_str_cat_uintle(io_buf, data_chunk_size + 1, 4);
	else
		rb_str_cat_uintle(io_buf, data_chunk_size, 4);

	io_writepartial(io, io_buf);
	
	io_buf = rb_str_new(0, 0);
	buffer_size = BUFFER_SIZE / block_size * block_size;
	idx = 0;
	for (long data_offset = 0; data_offset < data_chunk_size; data_offset += buffer_size)
	{
		if ((1. * data_offset + buffer_size) > data_chunk_size)
			buffer_size = data_chunk_size - data_offset;
		
		if (data_offset == 0)
			io_buf = rb_str_buf_z_new(buffer_size);
		else if (RSTRING_LEN(io_buf) != buffer_size)
			rb_str_buf_z_resize(io_buf, buffer_size);
		
		unsigned char *buf_ptr = (unsigned char *)RSTRING_PTR(io_buf);
		for ( ; ; )
		{
			for (long i = 0; i < channels; i++)
			{
				double *s_ptr = mat[i];
				func(buf_ptr+(i*block_size/channels), s_ptr+idx);
			}
			idx += samples_per_block;
			buf_ptr += block_size;
			
			if (buf_ptr == (unsigned char *)RSTRING_END(io_buf))
			{
				io_writepartial(io, io_buf);
				break;
			}
		}
		
	}
	if (data_chunk_size % 2 == 1)
		io_writepartial(io, rb_str_buf_z_new(1));
		
	rb_str_resize(io_buf, 0);
	rb_io_close(io);
	
	return Qtrue; // TODO: must be return a wrote byte-size
}


static VALUE
test_wave_write_linear_pcm(VALUE unused_obj, VALUE fname, VALUE pcm_ary, VALUE bits)
{
	if (TYPE(pcm_ary) != T_ARRAY)
		rb_raise(rb_eTypeError, "not an Array");
	return wave_write_linear_pcm(pcm_ary, NUM2INT(bits), StringValuePtr(fname));
}

void
InitVM_RIFF(void)
{
	rb_define_const(rb_cWaveRIFF, "SupportedVersion", rb_str_new_cstr(SupportedVersion));
	rb_define_singleton_method(rb_cWaveRIFF, "write_linear_pcm", test_wave_write_linear_pcm, 3);
	rb_define_singleton_method(rb_cWaveRIFF, "read_linear_pcm", test_wave_read_linear_pcm, 1);
}

