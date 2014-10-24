#ifndef MC04TYPES_H
#define MC04TYPES_H

/* MISRA C 2004 compliant numeric typedef */

typedef char char_t;

/* Probably already defined by inttypes.h or types.h*/
# ifndef __int8_t_defined
#  define __int8_t_defined
typedef signed char int8_t;
typedef signed short int16_t;
typedef signed int int32_t;
#endif

#ifndef _SYS_TYPES_H
#ifndef _INT64_T
#define _INT64_T
typedef signed long int64_t;
#endif
#endif

typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int uint32_t;
#ifndef _UINT64_T
#define _UINT64_T
typedef unsigned long uint64_t;
#endif
typedef float float32_t;
typedef double float64_t;
typedef long double float128_t;

#endif
