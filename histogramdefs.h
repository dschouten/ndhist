
#ifndef HISTOGRAM_DEFS_H
#define HISTOGRAM_DEFS_H

#include <histogram.h>

using ::hepstd::h3dd;
using ::hepstd::h4dd;
using ::hepstd::h5dd;

// --------------===============--------------===============--------------===============--------------===============--------------===============
// --------------===============--------------===============--------------===============--------------===============--------------===============

class hist3d : public h3dd
{
  
  hist3d( ) {}
public:  
  hist3d( const h3dd& H ) : h3dd( H ) { }
  
  bin_t get_bin			( int I, int J, int K ) const			{ return h3dd::get_bin( binselect_t(I,J,K) ); }
  bin_t get_bin			( double X, double Y, double Z ) const		{ return h3dd::get_bin( point_t(X,Y,Z) ); } 
  
  entry_t get_bin_content	( int I, int J, int K ) const			{ return h3dd::get_bin_content( h3dd::get_bin( binselect_t(I,J,K) ) ); }
  entry_t get_bin_error		( int I, int J, int K ) const			{ return h3dd::get_bin_error( h3dd::get_bin( binselect_t(I,J,K) ) ); }  
  
  entry_t get_integral		( double X1, double Y1, double Z1, 
				  double X2, double Y2, double Z2,
				  bool overflow = false ) const			{ return h3dd::get_integral( h3dd::get_bin( point_t(X1,Y1,Z1) ),
													     h3dd::get_bin( point_t(X2,Y2,Z2) ), 
													     overflow ); }
  entry_t get_integral          ( int I1, int J1, int K1, 
				  int I2, int J2, int K2, 
				  bool overflow = false ) const			{ return h3dd::get_integral( h3dd::get_bin( binselect_t(I1,J1,K1) ),
													     h3dd::get_bin( binselect_t(I2,J2,K2) ), 
													     overflow ); }
  
  void set_bin_content		( int I, int K, int J, const entry_t& w )	{ h3dd::set_bin_content( h3dd::get_bin( binselect_t(I,J,K) ), w ); }
  void set_bin_error		( int I, int K, int J, const entry_t& w )	{ h3dd::set_bin_error( h3dd::get_bin( binselect_t(I,J,K) ), w ); }
};

// --------------===============--------------===============--------------===============--------------===============--------------===============
// --------------===============--------------===============--------------===============--------------===============--------------===============

class hist4d : public h4dd
{
  hist4d( ) {}
public:  
  hist4d( const h4dd& H ) : h4dd( H ) { }
  
  bin_t get_bin			( int I, int J, int K, int L ) const			{ return h4dd::get_bin( binselect_t(I,J,K,L) ); }
  bin_t get_bin			( double W, double X, double Y, double Z ) const	{ return h4dd::get_bin( point_t(W,X,Y,Z) ); } 
  
  entry_t get_bin_content	( int I, int J, int K, int L ) const			{ return h4dd::get_bin_content( h4dd::get_bin( binselect_t(I,J,K,L) ) ); }
  entry_t get_bin_error		( int I, int J, int K, int L ) const			{ return h4dd::get_bin_error( h4dd::get_bin( binselect_t(I,J,K,L) ) ); }  
  
  entry_t get_integral		( double W1, double X1, double Y1, double Z1, 
				  double W2, double X2, double Y2, double Z2,
				  bool overflow = false ) const				{ return h4dd::get_integral( h4dd::get_bin( point_t(W1,X1,Y1,Z1) ),
														     h4dd::get_bin( point_t(W2,X2,Y2,Z2) ), 
														     overflow ); }
  
  entry_t get_integral          ( int I1, int J1, int K1, int L1, 
				  int I2, int J2, int K2, int L2,
				  bool overflow = false ) const				{ return h4dd::get_integral( h4dd::get_bin( binselect_t(I1,J1,K1,L1) ),
														     h4dd::get_bin( binselect_t(I2,J2,K2,L2) ), 
														     overflow ); }
  
  void set_bin_content		( int I, int K, int J, int L, const entry_t& w )	{ h4dd::set_bin_content( h4dd::get_bin( binselect_t(I,J,K,L) ), w ); }
  void set_bin_error		( int I, int K, int J, int L, const entry_t& w )	{ h4dd::set_bin_error( h4dd::get_bin( binselect_t(I,J,K,L) ), w ); }
};

// --------------===============--------------===============--------------===============--------------===============--------------===============
// --------------===============--------------===============--------------===============--------------===============--------------===============

class hist5d : public h5dd
{
  hist5d( ) {}
public:  
  hist5d( const h5dd& H ) : h5dd( H ) { }
  
  bin_t get_bin			( int I, int J, int K, int L, int M ) const			{ return h5dd::get_bin( binselect_t(I,J,K,L,M) ); }
  bin_t get_bin			( double V, double W, double X, double Y, double Z ) const	{ return h5dd::get_bin( point_t(V,W,X,Y,Z) ); } 
  
  entry_t get_bin_content	( int I, int J, int K, int L, int M ) const			{ return h5dd::get_bin_content( h5dd::get_bin( binselect_t(I,J,K,L,M) ) ); }
  entry_t get_bin_error		( int I, int J, int K, int L, int M ) const			{ return h5dd::get_bin_error( h5dd::get_bin( binselect_t(I,J,K,L,M) ) ); }  
  
  entry_t get_integral		( double V1, double W1, double X1, double Y1, double Z1, 
				  double V2, double W2, double X2, double Y2, double Z2,
				  bool overflow = false ) const					{ return h5dd::get_integral( h5dd::get_bin( point_t(V1,W1,X1,Y1,Z1) ),
															     h5dd::get_bin( point_t(V2,W2,X2,Y2,Z2) ), 
															     overflow ); }
  entry_t get_integral          ( int I1, int J1, int K1, int L1, int M1,
				  int I2, int J2, int K2, int L2, int M2,
				  bool overflow = false ) const					{ return h5dd::get_integral( h5dd::get_bin( binselect_t(I1,J1,K1,L1,M1) ),
															     h5dd::get_bin( binselect_t(I2,J2,K2,L2,M2) ), 
															     overflow ); }
  
  void set_bin_content		( int I, int K, int J, int L, int M, const entry_t& w )		{ h5dd::set_bin_content( h5dd::get_bin( binselect_t(I,J,K,L,M) ), w ); }
  void set_bin_error		( int I, int K, int J, int L, int M, const entry_t& w )		{ h5dd::set_bin_error( h5dd::get_bin( binselect_t(I,J,K,L,M) ), w ); }
};

#endif
