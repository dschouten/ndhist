// @author doug schouten <doug dot schouten at triumf dot ca>
// @file tuple.h

#ifndef TUPLE_H
#define TUPLE_H

#include <vector>
#include <map>
#include <stdexcept>
#include <cassert>
#include <cstdarg>
#include <cstdio>
#include <iostream>

#include <cppmacros.h>

#include <boost/array.hpp>

/*
 * tuple<N,type> 
 *
 * Class to represent a tuple (of fixed length N) of objects
 * of a uniform type
 *
 * Currently, maximum dimension is set to TUPLE_MAX_NDIMENSIONS = 10 
 * unless using the C11 standard
 *
 */

#ifndef __HAVE_C11_STD
#define TUPLE_MAX_NDIMENSIONS 10
#else
#define TUPLE_MAX_NDIMENSIONS boost::integer_traits<int>::const_max
#endif

namespace hepstd {
  
  template <unsigned int NDIM, typename DATA = float>
  class tuple 
  {
  public:
    typedef DATA data_t;
    static const unsigned int ndimension = NDIM;
    
  public:
#ifdef __HAVE_C11_STD
    
    // @todo
    
#else
#define TPL_CARG(z, n, empty) const DATA& X##n,
#define TPL_CTOR(z, n, empty)			\
    tuple( BOOST_PP_REPEAT( n, TPL_CARG, 0 )	\
	   const DATA& X##n ); 
    BOOST_PP_REPEAT( TUPLE_MAX_NDIMENSIONS, TPL_CTOR, 0 );
#undef TPL_CTOR
#undef TPL_CARG  
    void assign( data_t*[] ) const; 
#endif
    tuple( const data_t*[] );
    tuple( const data_t* );
    tuple( const std::vector<data_t>& );
    
    virtual ~tuple( ) { }
    
  public:
    
    //
    // API
    //
    
    const data_t& operator[]( unsigned int ) const;
    const data_t& at( unsigned int ) const;
    
    typedef typename boost::array< data_t, NDIM >::const_iterator iterator_t;
    typedef typename boost::array< data_t, NDIM >::const_reverse_iterator reverse_iterator_t;
    
    iterator_t begin() { return _arry.begin(); }
    iterator_t end() { return _arry.end(); }
    
    reverse_iterator_t rbegin() { return _arry.rbegin(); }
    reverse_iterator_t rend() { return _arry.rend(); }
    
    
  protected:
    boost::array< data_t, NDIM > _arry;
    
  private:
    void _init( const data_t*, unsigned int );
  };
  
  // --------------===============--------------===============--------------===============
  
#ifdef __HAVE_C11_STD
  
  // @todo
  
#else
#define TPL_CARG(z, n, delim) const DATA& X##n,
#define TPL_CTOR(n)				\
  tuple<NDIM,DATA>::tuple(			\
			  BOOST_PP_REPEAT( n, TPL_CARG, 0 )	\
			  const DATA& X##n			\
						) 
#define TPL_CSET(z, k, empty) if( ndimension >= k+1 ) arr[k] = X##k;
#define TPL_CTOR_FUN(z, n, empty)		     \
  template<unsigned int NDIM, typename DATA>	     \
  TPL_CTOR(n)					     \
  {							  \
    STATIC_CHECK((NDIM>0&&NDIM<=TUPLE_MAX_NDIMENSIONS),   \
		 Tuple_Dimension_Invalid);		  \
    data_t arr[n];					  \
    data_t arr_dummy[ndimension];			  \
    STATIC_CHECK( (sizeof(arr)<sizeof(arr_dummy)),	  \
		  Dimension_Items_Mismatch );		  \
    BOOST_PP_REPEAT(n, TPL_CSET, 0);			  \
    if( ndimension >= n ) arr[n] = X##n;		  \
    _init( arr, n+1 );					  \
  }
  BOOST_PP_REPEAT( TUPLE_MAX_NDIMENSIONS, TPL_CTOR_FUN, 0 );
  
#undef TPL_CTOR_FUN
#undef TPL_CSET
#undef TPL_CTOR
#undef TPL_CARG
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA>
  void tuple<NDIM,DATA>::assign( tuple<NDIM,DATA>::data_t* x[] ) const
  {
    for( unsigned int iarg = 0; iarg < ndimension; ++iarg )
    {
      (*(x[iarg])) = _arry[iarg];
    }
  }
#endif
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA>
  tuple<NDIM,DATA>::tuple( const tuple<NDIM,DATA>::data_t* p_arr[] )
  {
    STATIC_CHECK((NDIM>0&&NDIM<=TUPLE_MAX_NDIMENSIONS),Tuple_Dimension_Invalid); 
    data_t v_arr[ndimension];
    for( unsigned int i = 0; i < ndimension; ++i )
    {
      v_arr[i] = *(p_arr[i]);
    }
    _init( v_arr, ndimension );
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA>
  tuple<NDIM,DATA>::tuple( const tuple<NDIM,DATA>::data_t* arr )
  {
    STATIC_CHECK((NDIM>0&&NDIM<=TUPLE_MAX_NDIMENSIONS),Tuple_Dimension_Invalid); 
    _init( arr, ndimension );
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA>
  tuple<NDIM,DATA>::tuple( const std::vector<tuple<NDIM,DATA>::data_t>& vec )
  {
    STATIC_CHECK((NDIM>0&&NDIM<=TUPLE_MAX_NDIMENSIONS),Tuple_Dimension_Invalid); 
    data_t* arr = new data_t[vec.size()];
    std::copy( vec.begin(), vec.end(), arr );
    _init( arr, ndimension );
    delete[] arr;
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA>
  void tuple<NDIM,DATA>::_init( const tuple<NDIM,DATA>::data_t* arr, unsigned int n )
  {
    if( n != ndimension )
    {
      sprintf( ::CPP::errors::buffer, "invalid initialization of tuple in %s @ L%d", __FILE__, __LINE__ );
      throw std::runtime_error( std::string(::CPP::errors::buffer).c_str() );
    }
    for( unsigned int i = 0; i < n; ++i )
    {
      _arry[i] = arr[i];
    }
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA>
  const typename tuple<NDIM,DATA>::data_t& tuple<NDIM,DATA>::operator[]( unsigned int i ) const
  {
    if( i >= ndimension )
    {
      sprintf( ::CPP::errors::buffer, "out of tuple bounds in %s @ L%d", __FILE__, __LINE__ );
      throw std::runtime_error( std::string(::CPP::errors::buffer).c_str() );
    }
    return _arry[i];
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA>
  const typename tuple<NDIM,DATA>::data_t& tuple<NDIM,DATA>::at( unsigned int i ) const
  {
    return (*this)[i];
  }
  
#ifndef __HAVE_C11_STD
  
  //
  // TUPLE_ASSIGN( mytuple, type, X0, X1, ..., XN )
  //
  // where N is the dimension of mytuple, will declare variables of type 
  // and fill them with the 0,1,...Nth corresponding entry in mytuple
  //
  
#define __TUPLE_REF( n, var ) &( var )
#define TUPLE_ASSIGN( tpl, typ, ... )					\
  typ __VA_ARGS__ ;							\
  do {									\
    auto typ *_arr_##tpl_A_##__LINE__[TUPLE_MAX_NDIMENSIONS] = {0x0};	\
    auto typ *_arr_##tpl_B_##__LINE__[tpl.ndimension] = { 0x0 };	\
    auto typ *_arr_##tpl_##__LINE__[] = { FOR_EACH( __TUPLE_REF, __VA_ARGS__ ) }; \
    STATIC_CHECK( sizeof(_arr_##tpl_A_##__LINE__) >=			\
		  sizeof(_arr_##tpl_##__LINE__), Tuple_Exceeds_Max_Dimension ); \
    STATIC_CHECK( sizeof(_arr_##tpl_B_##__LINE__) ==			\
		  sizeof(_arr_##tpl_##__LINE__), Tuple_Dimension_Items_Mismatch ); \
    tpl.assign( _arr_##tpl_##__LINE__ );				\
  } while(0)
#endif 
  
  typedef tuple<1,double> double_single_t;
  typedef tuple<2,double> double_couple_t;
  typedef tuple<3,double> double_triple_t;
  typedef tuple<4,double> double_quadruple_t;
  typedef tuple<5,double> double_quintuple_t;
  typedef tuple<6,double> double_sextuple_t;
  typedef tuple<7,double> double_septuple_t;
  typedef tuple<8,double> double_octuple_t;
  typedef tuple<9,double> double_nonuple_t;
  typedef tuple<10,double> double_decuple_t;

  typedef tuple<1,float> float_single_t;
  typedef tuple<2,float> float_couple_t;
  typedef tuple<3,float> float_triple_t;
  typedef tuple<4,float> float_quadruple_t;
  typedef tuple<5,float> float_quintuple_t;
  typedef tuple<6,float> float_sextuple_t;
  typedef tuple<7,float> float_septuple_t;
  typedef tuple<8,float> float_octuple_t;
  typedef tuple<9,float> float_nonuple_t;
  typedef tuple<10,float> float_decuple_t;
  
  typedef tuple<1,int> int_single_t;
  typedef tuple<2,int> int_couple_t;
  typedef tuple<3,int> int_triple_t;
  typedef tuple<4,int> int_quadruple_t;
  typedef tuple<5,int> int_quintuple_t;
  typedef tuple<6,int> int_sextuple_t;
  typedef tuple<7,int> int_septuple_t;
  typedef tuple<8,int> int_octuple_t;
  typedef tuple<9,int> int_nonuple_t;
  typedef tuple<10,int> int_decuple_t;
  
  typedef tuple<1,unsigned int> uint_single_t;
  typedef tuple<2,unsigned int> uint_couple_t;
  typedef tuple<3,unsigned int> uint_triple_t;
  typedef tuple<4,unsigned int> uint_quadruple_t;
  typedef tuple<5,unsigned int> uint_quintuple_t;
  typedef tuple<6,unsigned int> uint_sextuple_t;
  typedef tuple<7,unsigned int> uint_septuple_t;
  typedef tuple<8,unsigned int> uint_octuple_t;
  typedef tuple<9,unsigned int> uint_nonuple_t;
  typedef tuple<10,unsigned int> uint_decuple_t;
  
}

#endif
