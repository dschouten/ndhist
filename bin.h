// @author doug schouten <doug dot schouten at triumf dot ca>
// @file bin.h

#ifndef BIN_H
#define BIN_H

#include <cstdarg>
#include <vector>
#include <map>
#include <stdexcept>
#include <assert.h>
#include <limits>

#include <tuple.h>

namespace hepstd {
  
  template <unsigned int NDIM, typename DATA>
  class bin
  {
  public:
    typedef unsigned int            size_t;
    typedef DATA                    data_t;
    typedef std::pair<DATA,DATA>    interval_t;
    typedef std::vector<interval_t> boundary_t;
    typedef tuple<NDIM,DATA>        tuple_t;
    
    static const unsigned int ndimension = NDIM;
    
  public:
    bin( ) : _edges() { }
    
    bin( DATA xlo, DATA xup, ... );
    bin( const std::vector<DATA>& lower, const std::vector<DATA>& upper );
    bin( const boundary_t& );
    
    bin( const tuple_t& );
    
    bin( const bin& B ) : _edges( NDIM ), _bits( B._bits ) { std::copy( B.begin(), B.end(), _edges.begin() ); }
    bin& operator=( const bin& B );
    virtual ~bin() { }
    
    template<typename P>
    bool operator<( const bin<NDIM,P>& ) const;
    
    template<typename P>
    bool contains( const tuple<NDIM,P>& ) const;
    
    template<typename P>
    bool contains( const bin<NDIM,P>& ) const;
    
    const boundary_t& boundary( ) const;
    
    tuple_t center( ) const; 
    
    interval_t operator[]( unsigned int ) const;
    
    void set_bit( int, bool );
    bool get_bit( int ) const;
    
  public:
    typedef typename boundary_t::const_iterator iterator_t;
    
    iterator_t begin() const { return _edges.begin(); }
    iterator_t end() const { return _edges.end(); }
    
  protected:
    boundary_t _edges;
    
  private:
    long _bits;
  };
  
  // --------------===============--------------===============--------------===============
  template <unsigned int NDIM, typename DATA>
  bin<NDIM,DATA>::bin( const bin<NDIM,DATA>::tuple_t& t )
  {
    DATA epsilon = 0;
    for( unsigned int i = 0; i < ndimension; ++i )
    {
      epsilon = std::numeric_limits<DATA>::epsilon() * t[i];
      _edges.push_back( interval_t( t[i], t[i] + epsilon ) );
    }
  }
  
  // --------------===============--------------===============--------------===============
  template <unsigned int NDIM, typename DATA>
  bin<NDIM,DATA>::bin( DATA xlo, DATA xup, ... )
  {
    va_list args;
    va_start( args, xup );
    _edges.push_back( interval_t( xlo, xup ) );
    for( unsigned int iarg = NDIM-1; iarg > 0; --iarg )
    {
      DATA lo = va_arg( args, DATA );
      DATA up = va_arg( args, DATA );
      _edges.push_back( interval_t( lo, up ) );
    }
    va_end( args );
    assert( _edges.size() == NDIM );
  }
  
  // --------------===============--------------===============--------------===============
  template <unsigned int NDIM, typename DATA>
  bin<NDIM,DATA>::bin( const std::vector<DATA>& lower, const std::vector<DATA>& upper ) : _edges( NDIM )
  {
    if( lower.size() != NDIM || upper.size() != NDIM )
    {
      sprintf( ::CPP::errors::buffer, "invalid initialization of bin %d / %d / %d in %s @ L%d", 
	       (int)lower.size(), (int)upper.size(), NDIM, __FILE__, __LINE__ );
      throw std::runtime_error( ::CPP::errors::buffer );
    }
    for( unsigned int i = 0; i < NDIM; ++i )
    {
      _edges[i] = interval_t( lower[i], upper[i] );
    }
  }
  
  // --------------===============--------------===============--------------===============
  template <unsigned int NDIM, typename DATA>
  bin<NDIM,DATA>::bin( const boundary_t& edges ) : _edges( NDIM )
  {
    if( edges.size() != NDIM )
    {
      sprintf( ::CPP::errors::buffer, "invalid initialization of bin : (%d != %d) in %s @ L%d", 
	       (int)edges.size(), NDIM, __FILE__, __LINE__ );
      throw std::runtime_error( ::CPP::errors::buffer );
    }
    std::copy( edges.begin(), edges.end(), _edges.begin() );
  }
  
  // --------------===============--------------===============--------------===============
  
  template <unsigned int NDIM, typename DATA>
  bin<NDIM,DATA>& bin<NDIM,DATA>::operator=( const bin<NDIM,DATA>& B ) 
  {
    _edges = B._edges;
    _bits  = B._bits;
    return (*this); 
  }
  
  // --------------===============--------------===============--------------===============
  template <unsigned int NDIM, typename DATA>
  template<typename P>
  bool bin<NDIM,DATA>::operator<( const bin<NDIM,P>& B ) const
  {
    return ( _edges < B.boundary() ); // both std::vector and std::pair < operators follow strict weak ordering 
  }
  
  // --------------===============--------------===============--------------===============
  template <unsigned int NDIM, typename DATA>
  template<typename P>
  bool bin<NDIM,DATA>::contains( const bin<NDIM,P>& B ) const
  {
    typename boundary_t::const_iterator itr = B.boundary().begin();
    typename boundary_t::const_iterator itr_end = B.boundary().end();
    unsigned idim = 0;
    for( ; itr != itr_end; ++itr )
    {
      if( (DATA)(itr->first) < _edges[idim].first || (DATA)(itr->second) > _edges[idim].second )
      {
	return false;
      }
      idim += 1;
    }
    return true;
  }
  
  // --------------===============--------------===============--------------===============
  template <unsigned int NDIM, typename DATA>
  template<typename P>
  bool bin<NDIM,DATA>::contains( const tuple<NDIM,P>& T ) const
  {
    return contains( bin<NDIM,P>( T ) );
  }
  
  // --------------===============--------------===============--------------===============
  template <unsigned int NDIM, typename DATA>
  const std::vector< std::pair<DATA,DATA> >& bin<NDIM,DATA>::boundary( ) const
  {
    return _edges;
  }
  
  // --------------===============--------------===============--------------===============
  template <unsigned int NDIM, typename DATA>
  tuple<NDIM,DATA> bin<NDIM,DATA>::center( ) const
  {
    std::vector<data_t> buff;
    for( unsigned idim = 0; idim < NDIM; ++idim )
    {
      buff.push_back( (_edges[idim].first + _edges[idim].second) / 2. );
    }
    return tuple_t( buff );
  }
  
  // --------------===============--------------===============--------------===============
  template <unsigned int NDIM, typename DATA>
  typename bin<NDIM,DATA>::interval_t bin<NDIM,DATA>::operator[]( unsigned int i ) const
  {
    return _edges[i]; // @todo bounds checking ...
  }
  
  // --------------===============--------------===============--------------===============
  template <unsigned int NDIM, typename DATA>
  void bin<NDIM,DATA>::set_bit( int ibit, bool bit ) 
  {
    long mask = ~(0x1<<ibit);
    if( !bit ) _bits = ( _bits & (~(0x1<<ibit)) );
    else       _bits = ( _bits | (0x1<<ibit) );
  }
  
  // --------------===============--------------===============--------------===============
  template <unsigned int NDIM, typename DATA>
  bool bin<NDIM,DATA>::get_bit( int ibit ) const
  {
    return ( ( _bits & (~(0x1<<ibit)) ) != _bits );
  }
  
}

#endif
