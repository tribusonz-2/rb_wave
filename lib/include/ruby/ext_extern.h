#ifndef RUBY_EXT_EXTERN_H_INCLUDED
#define RUBY_EXT_EXTERN_H_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#ifdef     USE_GLOBAL_VARIABLE
# define   RUBY_EXT_EXTERN
# define   GLOBAL_VAL(v) = (v)
#else
# define   RUBY_EXT_EXTERN    extern
# define   GLOBAL_VAL(v)
#endif

#if defined(__cplusplus)
}
#endif

#endif /* RUBY_EXT_EXTERN_H_INCLUDED */
