#ifndef CPP_STATIC_CHECK_H
#define CPP_STATIC_CHECK_H

#include <boost/static_assert.hpp>
#include <boost/integer_traits.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

#if defined(__GXX_EXPERIMENTAL_CXX0X) || __cplusplus >= 201103L
#define __HAVE_C11_STD 
#endif

#if defined(__MAKECINT__) || defined(G__DICTIONARY) 

#include <iostream>
#include <cassert>
#include <cmath>

#define STATIC_CHECK(expr, msg) \
           if (!(expr) ) std::cerr << "ERROR:   "  << #msg << std::endl; \
           assert(expr);

#else 

namespace CPP
{

  // template<typename TYPE>
  // TYPE round(TYPE r, unsigned decimal = 0) 
  // {
  //   unsigned long s = (unsigned long) pow(10., (int)decimal);
  //   TYPE sr = s*r;
  //   return (sr > 0.0) ? floor(sr + 0.5) : ceil(sr - 0.5);
  // }

  namespace errors
  {
    char buffer[1024];
  }

#ifndef USE_OLD_SC


     template<bool> struct CompileTimeChecker 
     {
        CompileTimeChecker(void *) {}
     };
     template<> struct CompileTimeChecker<false> {};

}  // end namespace CPP

#define STATIC_CHECK(expr, msg) \
   { class ERROR_##msg {}; \
   ERROR_##msg e; \
   (void) (CPP::CompileTimeChecker<(expr) != 0> (&e)); }


#else 

////////////////////////////////////////////////////////////////////////////////
// Helper structure for the STATIC_CHECK macro
////////////////////////////////////////////////////////////////////////////////

    template<int> struct CompileTimeError;
    template<> struct CompileTimeError<true> {};

}  // end namespace CPP

////////////////////////////////////////////////////////////////////////////////
// macro STATIC_CHECK
// Invocation: STATIC_CHECK(expr, id)
// where:
// expr is a compile-time integral or pointer expression
// id is a C++ identifier that does not need to be defined
// If expr is zero, id will appear in a compile-time error message.
////////////////////////////////////////////////////////////////////////////////

#define STATIC_CHECK(expr, msg) \
    { CPP::CompileTimeError<((expr) != 0)> ERROR_##msg; (void)ERROR_##msg; } 

#endif

#endif 

#endif // CPP_STATIC_CHECK_H

////////////////////////////////////////////////////////////////////////////////
// macro FOR_EACH( F, ... )
//
// evaluate preprocessor macro F for each argument in list
////////////////////////////////////////////////////////////////////////////////

#define CONCATENATE(arg1, arg2)   arg1##arg2 

#define FOR_EACH_1(what, x, ...) what(1,x)
#define FOR_EACH_2(what, x, ...) \
  what(2,x),			 \
  FOR_EACH_1(what,  __VA_ARGS__)
#define FOR_EACH_3(what, x, ...) \
  what(3,x),			 \
  FOR_EACH_2(what, __VA_ARGS__)
#define FOR_EACH_4(what, x, ...) \
  what(4,x),			 \
  FOR_EACH_3(what,  __VA_ARGS__)
#define FOR_EACH_5(what, x, ...) \
  what(5,x),			 \
 FOR_EACH_4(what,  __VA_ARGS__)
#define FOR_EACH_6(what, x, ...) \
  what(6,x),			 \
  FOR_EACH_5(what,  __VA_ARGS__)
#define FOR_EACH_7(what, x, ...) \
  what(7,x),			 \
  FOR_EACH_6(what,  __VA_ARGS__)
#define FOR_EACH_8(what, x, ...) \
  what(8,x),			 \
  FOR_EACH_7(what,  __VA_ARGS__)
#define FOR_EACH_9(what, x, ...) \
  what(8,x),			 \
  FOR_EACH_8(what,  __VA_ARGS__)

#define FOR_EACH_NARG(...) FOR_EACH_NARG_(__VA_ARGS__, FOR_EACH_RSEQ_N())
#define FOR_EACH_NARG_(...) FOR_EACH_ARG_N(__VA_ARGS__) 
#define FOR_EACH_ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, N, ...) N 
#define FOR_EACH_RSEQ_N() 9, 8, 7, 6, 5, 4, 3, 2, 1, 0

#define FOR_EACH_(N, what, x, ...) CONCATENATE(FOR_EACH_, N)(what, x, __VA_ARGS__)
#define FOR_EACH(what, x, ...) FOR_EACH_(FOR_EACH_NARG(x, __VA_ARGS__), what, x, __VA_ARGS__)

