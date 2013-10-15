
#ifndef HISTOGRAM_DEFS_H
#define HISTOGRAM_DEFS_H

#include <histogram.h>

using ::hepstd::h1dd;
using ::hepstd::h2dd;
using ::hepstd::h3dd;
using ::hepstd::h4dd;
using ::hepstd::h5dd;

using ::hepstd::h1ff;
using ::hepstd::h2ff;
using ::hepstd::h3ff;
using ::hepstd::h4ff;
using ::hepstd::h5ff;

// ==========================================================================================

class hist1d : public h1dd
{
  hist1d() {}
public:  
  hist1d(unsigned nx, data_t* xbins) 
  {
    ::boost::array<unsigned int, 1> nbins_arr = { {nx} };
    ::boost::array<data_t*, 1> bins_arr = { {xbins} };
    h1dd::_set(nbins_arr, bins_arr);
  }
  
  hist1d(unsigned nx, data_t xlow, data_t xhigh)
  {
    ::boost::array<unsigned int, 1> nbins_arr = { {nx} };
    ::boost::array<data_t*, 1> bins_arr;
    bins_arr[0] = new data_t[nx+1];
    for(unsigned int ix=0; ix < nx+1; ++ix)
    {
      bins_arr[0][ix] = xlow + ix*((xhigh-xlow)/nx);
    }
    h1dd::_set(nbins_arr, bins_arr);
    delete[] bins_arr[0];  
  }
  
  hist1d(const h1dd& H) : h1dd(H) { }
  
  bin_t get_bin			(int I) const			{ return h1dd::get_bin(binselect_t(I)); }
  bin_t get_bin			(double X) const		{ return h1dd::get_bin(point_t(X)); } 
  
  entry_t get_bin_content	(int I) const			{ return h1dd::get_bin_content(h1dd::get_bin(binselect_t(I))); }
  entry_t get_bin_error		(int I) const			{ return h1dd::get_bin_error(h1dd::get_bin(binselect_t(I))); }  
  
  entry_t get_integral		(double X1, 
				 double X2,
				 bool overflow = false) const	{ return h1dd::get_integral(h1dd::get_bin(point_t(X1)),
											    h1dd::get_bin(point_t(X2)), 
											    overflow); }
  entry_t get_integral          (int I1,
				 int I2,
				 bool overflow = false) const	{ return h1dd::get_integral(h1dd::get_bin(binselect_t(I1)),
											    h1dd::get_bin(binselect_t(I2)), 
											    overflow); }
  
  void set_bin_content		(int I, const entry_t& w)	{ h1dd::set_bin_content(h1dd::get_bin(binselect_t(I)), w); }
  void set_bin_error		(int I, const entry_t& w)	{ h1dd::set_bin_error(h1dd::get_bin(binselect_t(I)), w); }
};

// ==========================================================================================

class hist2d : public h2dd
{
  hist2d() {}
public:  
  hist2d(unsigned nx, data_t* xbins, unsigned ny, data_t* ybins) 
  {
    ::boost::array<unsigned int,2> nbins_arr = { {nx, ny} };
    ::boost::array<data_t*, 2> bins_arr = { {xbins, ybins} };
    h2dd::_set(nbins_arr, bins_arr);
  }
  
  hist2d(unsigned nx, data_t xlow, data_t xhigh, unsigned ny, data_t ylow, data_t yhigh)
  {
    ::boost::array<unsigned int,2> nbins_arr = { {nx, ny} };
    ::boost::array<data_t*, 2> bins_arr;
    bins_arr[0] = new data_t[nx+1];
    for(unsigned int ix=0; ix < nx+1; ++ix)
    {
      bins_arr[0][ix] = xlow + ix*((xhigh-xlow)/nx);
    }
    bins_arr[1] = new data_t[ny+1];
    for(unsigned int iy=0; iy < ny+1; ++iy)
    {
      bins_arr[0][iy] = ylow + iy*((yhigh-ylow)/ny);
    }
    h2dd::_set(nbins_arr, bins_arr);
    delete[] bins_arr[0];
    delete[] bins_arr[1];    
  }
  
  hist2d(const h2dd& H) : h2dd(H) { }
  
  bin_t get_bin			(int I, int J) const			{ return h2dd::get_bin(binselect_t(I,J)); }
  bin_t get_bin			(double X, double Y) const		{ return h2dd::get_bin(point_t(X,Y)); } 
  
  entry_t get_bin_content	(int I, int J) const			{ return h2dd::get_bin_content(h2dd::get_bin(binselect_t(I,J))); }
  entry_t get_bin_error		(int I, int J) const			{ return h2dd::get_bin_error(h2dd::get_bin(binselect_t(I,J))); }  
  
  entry_t get_integral		(double X1, double Y1, 
				 double X2, double Y2,
				 bool overflow = false) const		{ return h2dd::get_integral(h2dd::get_bin(point_t(X1,Y1)),
												    h2dd::get_bin(point_t(X2,Y2)), 
												    overflow); }
  entry_t get_integral          (int I1, int J1,
				 int I2, int J2,
				 bool overflow = false) const		{ return h2dd::get_integral(h2dd::get_bin(binselect_t(I1,J1)),
												    h2dd::get_bin(binselect_t(I2,J2)), 
												    overflow); }
  
  void set_bin_content		(int I, int J, const entry_t& w)	{ h2dd::set_bin_content(h2dd::get_bin(binselect_t(I,J)), w); }
  void set_bin_error		(int I, int J, const entry_t& w)	{ h2dd::set_bin_error(h2dd::get_bin(binselect_t(I,J)), w); }
};

// ==========================================================================================

class hist3d : public h3dd
{
  
  hist3d() {}
public:  
  hist3d(const h3dd& H) : h3dd(H) { }
  
  bin_t get_bin			(int I, int J, int K) const			{ return h3dd::get_bin(binselect_t(I,J,K)); }
  bin_t get_bin			(double X, double Y, double Z) const		{ return h3dd::get_bin(point_t(X,Y,Z)); } 
  
  entry_t get_bin_content	(int I, int J, int K) const			{ return h3dd::get_bin_content(h3dd::get_bin(binselect_t(I,J,K))); }
  entry_t get_bin_error		(int I, int J, int K) const			{ return h3dd::get_bin_error(h3dd::get_bin(binselect_t(I,J,K))); }  
  
  entry_t get_integral		(double X1, double Y1, double Z1, 
				 double X2, double Y2, double Z2,
				 bool overflow = false) const			{ return h3dd::get_integral(h3dd::get_bin(point_t(X1,Y1,Z1)),
													    h3dd::get_bin(point_t(X2,Y2,Z2)), 
													    overflow); }
  entry_t get_integral          (int I1, int J1, int K1, 
				 int I2, int J2, int K2, 
				 bool overflow = false) const			{ return h3dd::get_integral(h3dd::get_bin(binselect_t(I1,J1,K1)),
													    h3dd::get_bin(binselect_t(I2,J2,K2)), 
													    overflow); }
  
  void set_bin_content		(int I, int K, int J, const entry_t& w)	{ h3dd::set_bin_content(h3dd::get_bin(binselect_t(I,J,K)), w); }
  void set_bin_error		(int I, int K, int J, const entry_t& w)	{ h3dd::set_bin_error(h3dd::get_bin(binselect_t(I,J,K)), w); }
};

// ==========================================================================================

class hist4d : public h4dd
{
  hist4d() {}
public:  
  hist4d(const h4dd& H) : h4dd(H) { }
  
  bin_t get_bin			(int I, int J, int K, int L) const			{ return h4dd::get_bin(binselect_t(I,J,K,L)); }
  bin_t get_bin			(double W, double X, double Y, double Z) const	{ return h4dd::get_bin(point_t(W,X,Y,Z)); } 
  
  entry_t get_bin_content	(int I, int J, int K, int L) const			{ return h4dd::get_bin_content(h4dd::get_bin(binselect_t(I,J,K,L))); }
  entry_t get_bin_error		(int I, int J, int K, int L) const			{ return h4dd::get_bin_error(h4dd::get_bin(binselect_t(I,J,K,L))); }  
  
  entry_t get_integral		(double W1, double X1, double Y1, double Z1, 
				 double W2, double X2, double Y2, double Z2,
				 bool overflow = false) const				{ return h4dd::get_integral(h4dd::get_bin(point_t(W1,X1,Y1,Z1)),
														    h4dd::get_bin(point_t(W2,X2,Y2,Z2)), 
														    overflow); }
  
  entry_t get_integral          (int I1, int J1, int K1, int L1, 
				 int I2, int J2, int K2, int L2,
				 bool overflow = false) const				{ return h4dd::get_integral(h4dd::get_bin(binselect_t(I1,J1,K1,L1)),
														    h4dd::get_bin(binselect_t(I2,J2,K2,L2)), 
														    overflow); }
  
  void set_bin_content		(int I, int K, int J, int L, const entry_t& w)	{ h4dd::set_bin_content(h4dd::get_bin(binselect_t(I,J,K,L)), w); }
  void set_bin_error		(int I, int K, int J, int L, const entry_t& w)	{ h4dd::set_bin_error(h4dd::get_bin(binselect_t(I,J,K,L)), w); }
};

// ==========================================================================================

class hist5d : public h5dd
{
  hist5d() {}
public:  
  hist5d(const h5dd& H) : h5dd(H) { }
  
  bin_t get_bin			(int I, int J, int K, int L, int M) const			{ return h5dd::get_bin(binselect_t(I,J,K,L,M)); }
  bin_t get_bin			(double V, double W, double X, double Y, double Z) const	{ return h5dd::get_bin(point_t(V,W,X,Y,Z)); } 
  
  entry_t get_bin_content	(int I, int J, int K, int L, int M) const			{ return h5dd::get_bin_content(h5dd::get_bin(binselect_t(I,J,K,L,M))); }
  entry_t get_bin_error		(int I, int J, int K, int L, int M) const			{ return h5dd::get_bin_error(h5dd::get_bin(binselect_t(I,J,K,L,M))); }  
  
  entry_t get_integral		(double V1, double W1, double X1, double Y1, double Z1, 
				 double V2, double W2, double X2, double Y2, double Z2,
				 bool overflow = false) const					{ return h5dd::get_integral(h5dd::get_bin(point_t(V1,W1,X1,Y1,Z1)),
															    h5dd::get_bin(point_t(V2,W2,X2,Y2,Z2)), 
															    overflow); }
  entry_t get_integral          (int I1, int J1, int K1, int L1, int M1,
				 int I2, int J2, int K2, int L2, int M2,
				 bool overflow = false) const					{ return h5dd::get_integral(h5dd::get_bin(binselect_t(I1,J1,K1,L1,M1)),
															    h5dd::get_bin(binselect_t(I2,J2,K2,L2,M2)), 
															    overflow); }
  
  void set_bin_content		(int I, int K, int J, int L, int M, const entry_t& w)		{ h5dd::set_bin_content(h5dd::get_bin(binselect_t(I,J,K,L,M)), w); }
  void set_bin_error		(int I, int K, int J, int L, int M, const entry_t& w)		{ h5dd::set_bin_error(h5dd::get_bin(binselect_t(I,J,K,L,M)), w); }
};

// ==========================================================================================

class hist1f : public h1ff
{
  hist1f() {}
public:  
  hist1f(unsigned nx, data_t* xbins) 
  {
    ::boost::array<unsigned int, 1> nbins_arr = { {nx} };
    ::boost::array<data_t*, 1> bins_arr = { {xbins} };
    h1ff::_set(nbins_arr, bins_arr);
  }
  
  hist1f(unsigned nx, data_t xlow, data_t xhigh)
  {
    ::boost::array<unsigned int, 1> nbins_arr = { {nx} };
    ::boost::array<data_t*, 1> bins_arr;
    bins_arr[0] = new data_t[nx+1];
    for(unsigned int ix=0; ix < nx+1; ++ix)
    {
      bins_arr[0][ix] = xlow + ix*((xhigh-xlow)/nx);
    }
    h1ff::_set(nbins_arr, bins_arr);
    delete[] bins_arr[0];  
  }
  
  hist1f(const h1ff& H) : h1ff(H) { }
  
  bin_t get_bin			(int I) const			{ return h1ff::get_bin(binselect_t(I)); }
  bin_t get_bin			(float X) const		{ return h1ff::get_bin(point_t(X)); } 
  
  entry_t get_bin_content	(int I) const			{ return h1ff::get_bin_content(h1ff::get_bin(binselect_t(I))); }
  entry_t get_bin_error		(int I) const			{ return h1ff::get_bin_error(h1ff::get_bin(binselect_t(I))); }  
  
  entry_t get_integral		(float X1, 
				 float X2,
				 bool overflow = false) const	{ return h1ff::get_integral(h1ff::get_bin(point_t(X1)),
											    h1ff::get_bin(point_t(X2)), 
											    overflow); }
  entry_t get_integral          (int I1,
				 int I2,
				 bool overflow = false) const	{ return h1ff::get_integral(h1ff::get_bin(binselect_t(I1)),
											    h1ff::get_bin(binselect_t(I2)), 
											    overflow); }
  
  void set_bin_content		(int I, const entry_t& w)	{ h1ff::set_bin_content(h1ff::get_bin(binselect_t(I)), w); }
  void set_bin_error		(int I, const entry_t& w)	{ h1ff::set_bin_error(h1ff::get_bin(binselect_t(I)), w); }
};

// ==========================================================================================

class hist2f : public h2ff
{
  hist2f() {}
public:  
  hist2f(unsigned nx, data_t* xbins, unsigned ny, data_t* ybins) 
  {
    ::boost::array<unsigned int,2> nbins_arr = { {nx, ny} };
    ::boost::array<data_t*, 2> bins_arr = { {xbins, ybins} };
    h2ff::_set(nbins_arr, bins_arr);
  }
  
  hist2f(unsigned nx, data_t xlow, data_t xhigh, unsigned ny, data_t ylow, data_t yhigh)
  {
    ::boost::array<unsigned int,2> nbins_arr = { {nx, ny} };
    ::boost::array<data_t*, 2> bins_arr;
    bins_arr[0] = new data_t[nx+1];
    for(unsigned int ix=0; ix < nx+1; ++ix)
    {
      bins_arr[0][ix] = xlow + ix*((xhigh-xlow)/nx);
    }
    bins_arr[1] = new data_t[ny+1];
    for(unsigned int iy=0; iy < ny+1; ++iy)
    {
      bins_arr[0][iy] = ylow + iy*((yhigh-ylow)/ny);
    }
    h2ff::_set(nbins_arr, bins_arr);
    delete[] bins_arr[0];
    delete[] bins_arr[1];    
  }
  
  hist2f(const h2ff& H) : h2ff(H) { }
  
  bin_t get_bin			(int I, int J) const			{ return h2ff::get_bin(binselect_t(I,J)); }
  bin_t get_bin			(float X, float Y) const		{ return h2ff::get_bin(point_t(X,Y)); } 
  
  entry_t get_bin_content	(int I, int J) const			{ return h2ff::get_bin_content(h2ff::get_bin(binselect_t(I,J))); }
  entry_t get_bin_error		(int I, int J) const			{ return h2ff::get_bin_error(h2ff::get_bin(binselect_t(I,J))); }  
  
  entry_t get_integral		(float X1, float Y1, 
				 float X2, float Y2,
				 bool overflow = false) const		{ return h2ff::get_integral(h2ff::get_bin(point_t(X1,Y1)),
												    h2ff::get_bin(point_t(X2,Y2)), 
												    overflow); }
  entry_t get_integral          (int I1, int J1,
				 int I2, int J2,
				 bool overflow = false) const		{ return h2ff::get_integral(h2ff::get_bin(binselect_t(I1,J1)),
												    h2ff::get_bin(binselect_t(I2,J2)), 
												    overflow); }
  
  void set_bin_content		(int I, int J, const entry_t& w)	{ h2ff::set_bin_content(h2ff::get_bin(binselect_t(I,J)), w); }
  void set_bin_error		(int I, int J, const entry_t& w)	{ h2ff::set_bin_error(h2ff::get_bin(binselect_t(I,J)), w); }
};

// ==========================================================================================

class hist3f : public h3ff
{
  
  hist3f() {}
public:  
  hist3f(const h3ff& H) : h3ff(H) { }
  
  bin_t get_bin			(int I, int J, int K) const			{ return h3ff::get_bin(binselect_t(I,J,K)); }
  bin_t get_bin			(float X, float Y, float Z) const		{ return h3ff::get_bin(point_t(X,Y,Z)); } 
  
  entry_t get_bin_content	(int I, int J, int K) const			{ return h3ff::get_bin_content(h3ff::get_bin(binselect_t(I,J,K))); }
  entry_t get_bin_error		(int I, int J, int K) const			{ return h3ff::get_bin_error(h3ff::get_bin(binselect_t(I,J,K))); }  
  
  entry_t get_integral		(float X1, float Y1, float Z1, 
				 float X2, float Y2, float Z2,
				 bool overflow = false) const			{ return h3ff::get_integral(h3ff::get_bin(point_t(X1,Y1,Z1)),
													    h3ff::get_bin(point_t(X2,Y2,Z2)), 
													    overflow); }
  entry_t get_integral          (int I1, int J1, int K1, 
				 int I2, int J2, int K2, 
				 bool overflow = false) const			{ return h3ff::get_integral(h3ff::get_bin(binselect_t(I1,J1,K1)),
													    h3ff::get_bin(binselect_t(I2,J2,K2)), 
													    overflow); }
  
  void set_bin_content		(int I, int K, int J, const entry_t& w)	{ h3ff::set_bin_content(h3ff::get_bin(binselect_t(I,J,K)), w); }
  void set_bin_error		(int I, int K, int J, const entry_t& w)	{ h3ff::set_bin_error(h3ff::get_bin(binselect_t(I,J,K)), w); }
};

// ==========================================================================================

class hist4f : public h4ff
{
  hist4f() {}
public:  
  hist4f(const h4ff& H) : h4ff(H) { }
  
  bin_t get_bin			(int I, int J, int K, int L) const			{ return h4ff::get_bin(binselect_t(I,J,K,L)); }
  bin_t get_bin			(float W, float X, float Y, float Z) const	{ return h4ff::get_bin(point_t(W,X,Y,Z)); } 
  
  entry_t get_bin_content	(int I, int J, int K, int L) const			{ return h4ff::get_bin_content(h4ff::get_bin(binselect_t(I,J,K,L))); }
  entry_t get_bin_error		(int I, int J, int K, int L) const			{ return h4ff::get_bin_error(h4ff::get_bin(binselect_t(I,J,K,L))); }  
  
  entry_t get_integral		(float W1, float X1, float Y1, float Z1, 
				 float W2, float X2, float Y2, float Z2,
				 bool overflow = false) const				{ return h4ff::get_integral(h4ff::get_bin(point_t(W1,X1,Y1,Z1)),
														    h4ff::get_bin(point_t(W2,X2,Y2,Z2)), 
														    overflow); }
  
  entry_t get_integral          (int I1, int J1, int K1, int L1, 
				 int I2, int J2, int K2, int L2,
				 bool overflow = false) const				{ return h4ff::get_integral(h4ff::get_bin(binselect_t(I1,J1,K1,L1)),
														    h4ff::get_bin(binselect_t(I2,J2,K2,L2)), 
														    overflow); }
  
  void set_bin_content		(int I, int K, int J, int L, const entry_t& w)	{ h4ff::set_bin_content(h4ff::get_bin(binselect_t(I,J,K,L)), w); }
  void set_bin_error		(int I, int K, int J, int L, const entry_t& w)	{ h4ff::set_bin_error(h4ff::get_bin(binselect_t(I,J,K,L)), w); }
};

// ==========================================================================================

class hist5f : public h5ff
{
  hist5f() {}
public:  
  hist5f(const h5ff& H) : h5ff(H) { }
  
  bin_t get_bin			(int I, int J, int K, int L, int M) const			{ return h5ff::get_bin(binselect_t(I,J,K,L,M)); }
  bin_t get_bin			(float V, float W, float X, float Y, float Z) const	{ return h5ff::get_bin(point_t(V,W,X,Y,Z)); } 
  
  entry_t get_bin_content	(int I, int J, int K, int L, int M) const			{ return h5ff::get_bin_content(h5ff::get_bin(binselect_t(I,J,K,L,M))); }
  entry_t get_bin_error		(int I, int J, int K, int L, int M) const			{ return h5ff::get_bin_error(h5ff::get_bin(binselect_t(I,J,K,L,M))); }  
  
  entry_t get_integral		(float V1, float W1, float X1, float Y1, float Z1, 
				 float V2, float W2, float X2, float Y2, float Z2,
				 bool overflow = false) const					{ return h5ff::get_integral(h5ff::get_bin(point_t(V1,W1,X1,Y1,Z1)),
															    h5ff::get_bin(point_t(V2,W2,X2,Y2,Z2)), 
															    overflow); }
  entry_t get_integral          (int I1, int J1, int K1, int L1, int M1,
				 int I2, int J2, int K2, int L2, int M2,
				 bool overflow = false) const					{ return h5ff::get_integral(h5ff::get_bin(binselect_t(I1,J1,K1,L1,M1)),
															    h5ff::get_bin(binselect_t(I2,J2,K2,L2,M2)), 
															    overflow); }
  
  void set_bin_content		(int I, int K, int J, int L, int M, const entry_t& w)		{ h5ff::set_bin_content(h5ff::get_bin(binselect_t(I,J,K,L,M)), w); }
  void set_bin_error		(int I, int K, int J, int L, int M, const entry_t& w)		{ h5ff::set_bin_error(h5ff::get_bin(binselect_t(I,J,K,L,M)), w); }
};

#endif
