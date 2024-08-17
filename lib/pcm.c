/*******************************************************************************
	pcm.c -

	$author$
*******************************************************************************/
#include <ruby.h>
#include "ruby/wave/globals.h"
#include "ruby/wave/pcm.h"

struct PCM {
	long fs;
	long length;
	double *s;
} ;

static struct PCM *
pcm_alloc(void)
{
	struct PCM *ptr = ALLOC(struct PCM);
	ptr->fs = FS_DEF;
	ptr->length = 0;
	ptr->s = NULL;
	return ptr;
}

static void
pcm_resize(struct PCM *ptr, long n)
{
	RUBY_ASSERT(ptr->length < 0);
	
	if (n < 0)
		rb_raise(rb_eRangeError, "negative (or biggest) sample size");
	else if (n == 0)
	{
		if (ptr->s != NULL)
			xfree(ptr->s);
		ptr->length = 0;
	}
	else /* if (n > 1) */
	{
		if (ptr->length != n)
		{
			if (ptr->length == 0)
			{
				RUBY_ASSERT(ptr->s != NULL);
				ptr->s = ALLOC_N(double, n);
				for (volatile long i = 0; i < n; i++)
					ptr->s[i] = 0.;
			}
			else if (ptr->length != n)
			{
				RUBY_ASSERT(ptr->s == NULL);
				REALLOC_N(ptr->s, double, n);
				if (ptr->length < n)
					for (volatile long i = ptr->length; i < n; i++)
						ptr->s[i] = 0.;
			}
			ptr->length = n;
		}
	}
}

static void
pcm_fs_set(struct PCM *ptr, long fs)
{
	RUBY_ASSERT(ptr->fs <= 0);
	
	if (fs <= 0)
		rb_raise(rb_eRangeError, "negative (or biggest) frequency");
	ptr->fs = fs;
}

static void
pcm_free(void *p)
{
	struct PCM *ptr = p;
	if (ptr->s != NULL)
		xfree(ptr->s);
	xfree(ptr);
}

static size_t
pcm_memsize(const void *p)
{
	size_t sz = sizeof(struct PCM);
	const struct PCM *ptr = p;
	sz += ptr->length * sizeof(double);
	return sz;
}

static const rb_data_type_t pcm_data_type = {
    "pcm",
    {
	0,
	pcm_free,
	pcm_memsize,
    },
    0, 0, RUBY_TYPED_FREE_IMMEDIATELY
};

#define check_pcm(self) ((struct PCM*)rb_check_typeddata((self), &pcm_data_type))

static struct PCM *
get_pcm(VALUE self)
{
    struct PCM *ptr = check_pcm(self);

    if (!ptr) {
	rb_raise(rb_eNoMemError, "uninitialized PCM");
    }
    return ptr;
}

static VALUE
pcm_s_allocate(VALUE klass)
{
	return TypedData_Wrap_Struct(klass, &pcm_data_type, 0);
}

/*
 *  call-seq:
 *    Wave::PCM.new(len, fs = Wave::PCM::FS_DEF) -> Wave::PCM
 *    Wave::PCM.new(len, fs = Wave::PCM::FS_DEF){|index| ...} -> Wave::PCM
 *  
 *  Create a new PCM class object with a length of +len+.
 *  The second argument +fs+ specifies the sampling frequency.
 *  The default value is constant of Wave::PCM::FS_DEF.
 *  
 *  If with block given, calls the block with each +index+ in the range (0...len)
 *  
 *    ```
 *    def sinewave(a, f0, n, fs)
 *      a * Math.sin(2.0 * Math::PI * f0 * n / fs)
 *    end
 *    
 *    a = 0.1
 *    f0 = 500.0
 *    FS = 8000
 *    
 *    pcm = Wave::PCM.new(16, FS){|n| sinewave(a, f0, n, FS)}
 *    pcm.to_a
 *    # => 
 *    [0.0,
 *     0.03826834323650898,
 *     0.07071067811865475,
 *     0.09238795325112868,
 *     0.1,
 *     0.0923879532511287,
 *     0.07071067811865477,
 *     0.03826834323650899,
 *     1.2246467991473533e-17,
 *     -0.03826834323650893,
 *     -0.07071067811865471,
 *     -0.09238795325112865,
 *     -0.1,
 *     -0.0923879532511287,
 *     -0.07071067811865477,
 *     -0.038268343236509045]
 *    ```
 */
static VALUE
rb_pcm_initialize(int argc, VALUE *argv, VALUE self)
{
	struct PCM *ptr = check_pcm(self);
	VALUE len, fs;
	
	if (!ptr)
		DATA_PTR(self) = ptr = pcm_alloc();
	
	rb_scan_args(argc, argv, "11", &len, &fs);
	if (argc == 1)  fs = LONG2FIX(FS_DEF);
	
	pcm_resize(ptr, NUM2LONG(len));
	pcm_fs_set(ptr, NUM2LONG(fs));
	
	if (rb_block_given_p())
	{
		for (volatile long i = 0; i < ptr->length; i++)
		{
			VALUE snd = rb_yield(LONG2NUM(i));
			ptr->s[i] = NUM2DBL(snd);
		}
	}
	
	return self;
}

/*
 *  call-seq:
 *    pcm.fs -> Integer
 *  
 *  Return PCM's sampling frequency.
 */
static VALUE
rb_pcm_fs_get(VALUE pcm)
{
	struct PCM *ptr = get_pcm(pcm);
	
	RUBY_ASSERT(ptr->fs <= 0);
	
	return LONG2NUM(ptr->fs);
}

/*
 *  call-seq:
 *    pcm.fs = fs
 *  
 *  Set PCM's sampling frequency. Assign the argument fs.
 */
static VALUE
rb_pcm_fs_set(VALUE pcm, VALUE fs)
{
	struct PCM *ptr = get_pcm(pcm);
	
	pcm_fs_set(ptr, NUM2LONG(fs));
	
	return fs;
}

/*
 *  call-seq:
 *    pcm.length -> Integer
 *  
 *  Return PCM's sample length.
 */
static VALUE
rb_pcm_len_get(VALUE pcm)
{
	struct PCM *ptr = get_pcm(pcm);
	
	RUBY_ASSERT(ptr->length < 0);
	
	return LONG2NUM(ptr->length);
}

/*
 *  call-seq:
 *    pcm.length = len
 *  
 *  Set PCM's sample length. 
 *  If +len+ is greater than the sample length itself, 
 *  memory is reallocated and initialized them to 0.0.
 */
static VALUE
rb_pcm_len_set(VALUE pcm, VALUE len)
{
	struct PCM *ptr = get_pcm(pcm);
	
	pcm_resize(ptr, NUM2LONG(len));
	
	return len;
}

/*
 *  call-seq:
 *    pcm[nth] -> Float | nil
 *  
 *  Returns the nth element. If the nth element does not exist, returns nil.
 *  Negative index values are supported, for example '-1' will returns the last value. The behavior is the same as the array class.
 */
static VALUE
pcm_snd_take1(struct PCM *ptr, long index)
{
	if (index >= 0)
	{
		if (ptr->length <= index)
			return Qnil;
		return DBL2NUM(ptr->s[index]);
	}
	else
	{
		if ((ptr->length + index) < 0)
			return Qnil;
		return DBL2NUM(ptr->s[ptr->length+index]);
	}
}

static VALUE
rb_pcm_snd_take(VALUE pcm, VALUE index)
{
	struct PCM *ptr = get_pcm(pcm);
	
	RUBY_ASSERT(ptr->length < 0);
	
	return pcm_snd_take1(ptr, NUM2LONG(index));

}


/*
 *  call-seq:
 *    pcm.eql?(other_pcm) -> bool
 *  
 *  Compare +self+ and +other_pcm+, then returns true if they are equal, false otherwise.
 */
static VALUE
rb_pcm_eql(VALUE pcm, VALUE other_pcm)
{
	struct PCM *lhs, *rhs;
	
	if (pcm == other_pcm)  return Qtrue;
	if (CLASS_OF(other_pcm) != rb_cWavePCM)  return Qfalse;
	
	lhs = check_pcm(pcm);
	rhs = check_pcm(other_pcm);
	
	if (lhs->fs != rhs->fs)  return Qfalse;
	if (lhs->length != rhs->length)  return Qfalse;
	for (volatile long i = 0; i < lhs->length; i++)
	{
		if (lhs->s[i] != rhs->s[i])
			return Qfalse;
	}
	return Qtrue;
}


static VALUE
pcm_enum_length(VALUE pcm, VALUE args, VALUE eobj)
{
	struct PCM *ptr = get_pcm(pcm);
	return LONG2FIX(ptr->length);
}

/*
 *  call-seq:
 *    pcm.each {|s| ... } -> self
 *    pcm.each -> Enumerator
 *  
 *  Iterates over the PCM-array each element.
 *  
 */
static VALUE
rb_pcm_each(VALUE pcm)
{
	struct PCM *ptr;
	
	RETURN_SIZED_ENUMERATOR(pcm, 0, 0, pcm_enum_length);
	
	ptr = get_pcm(pcm);
	for (volatile long i = 0; i < ptr->length; i++)
	{
		rb_yield(DBL2NUM(ptr->s[i]));
	}
	return pcm;
}

/*
 *  call-seq:
 *    pcm.map! {|s| ... } -> self
 *    pcm.map! -> Enumerator
 *  
 *  If block given, call this with each element, 
 *  and replaces the element with the block's return value.
 *  The replaced value must be a Float. Otherwise, an implicit type conversion will be attempted.
 *  
 *    ```
 *    def sawwave(a, f0, n, fs)
 *      s = 0.0
 *      (1..44).each do |i|
 *        s += 1.0 / i * Math.sin(2.0 * Math::PI * i * f0 * n / fs)
 *      end
 *      s *= a
 *    end
 *    
 *    gain = 0.1
 *    f0 = 500.0
 *    FS = 8000
 *    pcm = Wave::PCM.new(16, FS)
 *    
 *    pcm.map!.with_index{|_, n| sawwave(gain, f0, n, FS)}
 *    pcm.to_a
 *    #=> 
 *    [0.0,
 *     0.13664126621701408,
 *     0.12054835839792039,
 *     0.09926938022785707,
 *     0.077404038161591,
 *     0.05778716585772436,
 *     0.03974046316278057,
 *     0.020757897445726897,
 *     1.4456153762464314e-17,
 *     -0.020757897445727,
 *     -0.03974046316278051,
 *     -0.05778716585772438,
 *     -0.07740403816159103,
 *     -0.09926938022785708,
 *     -0.12054835839792005,
 *     -0.13664126621701386]
 *    ```
 */
static VALUE
rb_pcm_collect_bang(VALUE pcm)
{
	struct PCM *ptr;
	
	RETURN_SIZED_ENUMERATOR(pcm, 0, 0, pcm_enum_length);
	
	ptr = get_pcm(pcm);
	for (volatile long i = 0; i < ptr->length; i++)
	{
		const double s = ptr->s[i];
		VALUE retval = rb_yield(DBL2NUM(s));
		ptr->s[i] = NUM2DBL(retval);
	}
	return pcm;
}



void
InitVM_PCM(void)
{
	
	rb_include_module(rb_cWavePCM, rb_mEnumerable);
	rb_define_alloc_func(rb_cWavePCM, pcm_s_allocate);
	
	rb_define_const(rb_cWavePCM, "FS_DEF", LONG2NUM(FS_DEF));
	
	rb_define_method(rb_cWavePCM, "initialize", rb_pcm_initialize, -1);
	
	rb_define_method(rb_cWavePCM, "fs", rb_pcm_fs_get, 0);
	rb_define_method(rb_cWavePCM, "fs=", rb_pcm_fs_set, 1);
	rb_define_method(rb_cWavePCM, "length", rb_pcm_len_get, 0);
	rb_define_method(rb_cWavePCM, "length=", rb_pcm_len_set, 1);
	rb_define_method(rb_cWavePCM, "[]", rb_pcm_snd_take, 1);
	
	rb_define_method(rb_cWavePCM, "eql?", rb_pcm_eql, 1);
	
	rb_define_method(rb_cWavePCM, "each", rb_pcm_each, 0);
	rb_define_method(rb_cWavePCM, "map!", rb_pcm_collect_bang, 0);
}

/*******************************************************************************
	For C API
*******************************************************************************/

VALUE
rb_pcm_new(long len, long fs)
{
	struct PCM *ptr;
	VALUE obj = TypedData_Make_Struct(rb_cWavePCM, struct PCM, &pcm_data_type, ptr);
	
	pcm_resize(ptr, len);
	pcm_fs_set(ptr, fs);
	
	return obj;
}


long
rb_pcm_fs(VALUE pcm)
{
	struct PCM *ptr = get_pcm(pcm);
	
	return ptr->fs;
}

long
rb_pcm_len(VALUE pcm)
{
	struct PCM *ptr = get_pcm(pcm);
	
	return ptr->length;
}

double *
rb_waveform_data_ptr(VALUE pcm)
{
	struct PCM *ptr = get_pcm(pcm);
	
	return ptr->s;
}
