#ifndef RB_WAVE_PCM_H_INCLUDED
#define RB_WAVE_PCM_H_INCLUDED
/**
 * @file
 * @author     $Author$
 */
#include <ruby/internal/value.h> // VALUE


/**
 * Default  sampling  frequency.   We  do  not  believe  that  the  default  for
 * the PCM class  is  44kHz:  The mainstream  has  now  changed.  Yes, don't  be
 * behind the times.
 */
#define FS_DEF  48000


/**
 *  Create a new Wave::PCM object in C level. 
 * 
 * @param[in]  len             Length of the waveform data.
 * @param[in]  fs              Sampling frequency.
 * @exception  rb_eRangeError  Parameter out of range.
 * @return     Wave::PCM object.
 */
VALUE rb_pcm_new(long len, long fs);

/* Variant for creation in 44kHz */
#define rb_pcm_44k_new(len)  rb_pcm_new(len, 44100)

/* Variant for creation in 48kHz */
#define rb_pcm_48k_new(len)  rb_pcm_new(len, 48000)

/**
 * Queries sampling frequency number of the waveform data.
 * 
 * @param[in]  pcm  Wave:PCM in question.
 * @return     Its number of sampling frequency.
 * @pre        `pcm` must be an instance of Wave::PCM.
 */
long rb_pcm_fs(VALUE pcm);

/**
 * Queries length of the waveform data.
 * 
 * @param[in]  pcm  Wave:PCM in question.
 * @return     Its number of elements.
 * @pre        `pcm` must be an instance of Wave::PCM.
 */
long rb_pcm_len(VALUE pcm);
#define RPCM_LEN  rb_pcm_len                            /* alias rb_pcm_len() */

/**
 * Pointer to  PCM class  waveform data.  Returns  the beginning  of  the array.
 * The implementation  is  double  type.   It  is  `NULL`  when  the  length  of
 * arrays is 0.
 * 
 * It does require a technique.
 * 
 * ```CXX
 * VALUE
 * pcm_print(VALUE pcm)
 * {
 *   double *s = WaveformDataPtr(pcm);
 *   for (long i = 0; i < RPCM_LEN(pcm); i++)
 *   {
 *     volatile const double value = s[i];
 *     rb_p(DBL2NUM(value)); // `p` stdout on Ruby level
 *   }
 *   return pcm;
 * }
 * ```
 * 
 * @param[in]  obj            Wave::PCM object.
 * @exception  rb_eTypeError  `obj` is not a Wave::PCM object.
 * @note       Do not free() to this pointer. It will core-dump.
 */
double *rb_waveform_data_ptr(VALUE obj);
#define WaveformDataPtr  rb_waveform_data_ptr


#endif /* RB_WAVE_PCM_H_INCLUDED */
