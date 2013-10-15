// @author doug schouten <doug dot schouten at triumf dot ca>
// @file histogram.h

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <iostream>
#include <cmath>
#include <map>
#include <memory>
#include <boost/array.hpp>

/*************************************************************************************
 *
 * histogram<N,D,W> - generic histogram class in N dimensions, in which the data 
 *                    is of type D and the count weights are of type W
 *
 * the histogram bins are allocated dynamically to conserve memory for high-dimension
 * histograms which are sparsely populated (for dense histograms, a simple N-d array
 * is may be more efficient).
 *
 * constructor is: histogram(array<unsigned,N> NBINS, array<D,N> EDGES)
 *
 * where NBINS is an array specifying the # of bins in each dimensions and EDGES
 * is an array of length NBINS + 1 in the range [a,b] defining the bin edges, eg.
 * { a, a + epsilon, ..., b }
 *
 * a histogram factory function is also defined with a simpler argument list for
 * symmetric binning
 *
 * hist_factory<N,D,W>(NBINS_{0}, a_{0}, b_{0}, ..., NBINS_{N}, a_{N}, b_{N})
 *
 *************************************************************************************/

#include <bin.h>
#include <tuple.h>

#include <histogram_base.h>

namespace hepstd
{
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  class histogram : public histogram_base<NDIM, DATA, WEIGHT>
  {
  public:
    typedef histogram_base<NDIM,DATA,WEIGHT> base_t;    
    typedef histogram<NDIM,DATA,WEIGHT> histogram_t;  
    
    typedef typename base_t::entry_t      entry_t;
    typedef typename base_t::data_t       data_t;
    typedef typename base_t::bin_t        bin_t;
    typedef typename base_t::point_t      point_t;
    typedef typename base_t::binselect_t  binselect_t;
    
  public:
    //
    // constructors
    //
    
    histogram(::boost::array<unsigned int,NDIM>, 
	      ::boost::array<data_t*,NDIM>); // << main constructor
    
    histogram(const histogram_t&); // << copy constructor
    histogram& operator=(const histogram_t&); // << assignment operator
    
    virtual ~histogram() { }
    
  protected:
    histogram() : base_t() {}
    
  public:
    
    //
    // API
    //
    
    virtual void fill(const point_t&, const entry_t& w = 1);
    
    virtual bin_t get_bin		(const binselect_t&) const; 
    virtual bin_t get_bin               (const point_t&) const; 
    
    virtual entry_t get_bin_content	(const bin_t&) const;
    virtual entry_t get_bin_error	(const bin_t&) const;
    
    virtual entry_t get_integral        (const bin_t&, const bin_t&, bool include_overflow = false) const; 
    virtual entry_t get_integral        (bool include_overlow = false) const; 
    
    virtual void set_bin_content	(const bin_t&, const entry_t&); 
    virtual void set_bin_error		(const bin_t&, const entry_t&); 
    
  public:
    
    //
    // allow to iterate over filled bins
    //
    
    typedef typename std::map<bin_t,entry_t>::const_iterator const_filled_bin_iterator;
    
    const_filled_bin_iterator begin() const { return  _counts.begin(); }
    const_filled_bin_iterator   end() const { return  _counts.end(); }
    
    const_filled_bin_iterator err_begin() const { return  _weights.begin(); }
    const_filled_bin_iterator   err_end() const { return  _weights.end(); }
    
  public:
    unsigned long get_num_bins	() const { return get_num_bins(false); } // total number of bins
    unsigned long get_num_bins	(bool only_filled) const; // total number of bins
    
    unsigned get_num_bins 	(unsigned idim) const { return _bin_definitions[idim].size(); } // number of bins along dimension 
    int get_ibin                (unsigned idim, data_t) const; // return index of bin along dimension that contains data value
    
    data_t get_bin_low		(unsigned idim, unsigned ibin) const { return ( (idim < NDIM && ibin < _bin_definitions[idim].size()) ? 
										_bin_definitions[idim][ibin].first : 
										pow(-1.0, 0.5) ); } // bin location along specified dimension
    data_t get_bin_high		(unsigned idim, unsigned ibin) const { return ( (idim < NDIM && ibin < _bin_definitions[idim].size()) ? 
										_bin_definitions[idim][ibin].second : 
										pow(-1.0, 0.5) ); }
    data_t get_bin_center	(unsigned idim, unsigned ibin) const { return ( (idim < NDIM && ibin < _bin_definitions[idim].size()) ? 
										get_bin_high(idim,ibin) + get_bin_low(idim,ibin) / 2 :
										pow(-1.0, 0.5) ); }
    
    bool get_use_sum_weights    () const { return _flag_sumwts; }
    void set_use_sum_weights    (bool flag = true) { _flag_sumwts = flag; if(!flag) _weights.clear(); }

    histogram<1,DATA,WEIGHT> unroll(data_t xlow, data_t xhigh, const vector< vector<bin_t> >& mapping) const; // unroll into 1D histogram from xlow to xhigh
    histogram<1,DATA,WEIGHT> unroll(data_t xlow, data_t xhigh, const vector< vector<binselect_t> >& mapping) const; // unroll into 1D histogram from xlow to xhigh
    
    bool add			(const histogram_t*, entry_t w = (entry_t)1.); // add histogram to self with weight factor
    void scale			(entry_t); // scale self
    bool rebin			(const tuple<NDIM,unsigned int>&); // re-define bin edges by some factor along each dimension
    void clear                  () { _counts.clear(); _weights.clear(); _entries.clear(); } // clear all histogrammed data
    
    //
    // define histogram algebra
    //
    
    histogram_t& operator*=(entry_t S)		        { this->scale(S); return (*this); }
    histogram_t& operator+=(const histogram_t& H)	{ this->add(&H); return (*this); }
    
    histogram_t operator+(const histogram_t& H) const	{ histogram_t B = (*this); B.add(&H,+1); return B; }
    histogram_t operator-(const histogram_t& H) const	{ histogram_t B = (*this); B.add(&H,-1); return B; }
    histogram_t operator*(entry_t S) const		{ histogram_t B = (*this); B.scale(S); return B; }
    
  protected:
    
    void _set(::boost::array<unsigned int,NDIM> nbins_arr,
	      ::boost::array<data_t*,NDIM> bins_arr);
    
    typedef std::map<bin_t,entry_t> map_t;
    std::map<bin_t,entry_t> _counts;
    std::map<bin_t,entry_t> _weights;
    std::map<bin_t,long> _entries;
    
    ::boost::array<typename bin_t::boundary_t, NDIM> _bin_definitions;
    
    bool _flag_sumwts;
    
    void _check_valid_bin(const bin_t&) const; // throws std::runtime_error
    bool _check_bin_alignment(const histogram_t&) const; // check if bin definitions for histograms are aligned
    
    // bool _next_bin(const ::boost::array<std::vector<std::pair<int,int> >, NDIM>&, std::vector<int>&) const;
  };
  
  // --------------===============--------------===============--------------===============
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  histogram<NDIM,DATA,WEIGHT>::histogram(::boost::array<unsigned int,NDIM> nbins_arr,
					 ::boost::array<data_t*,NDIM> bins_arr) : base_t(),
										  _flag_sumwts(true)
  { 
    try
    {
      _set(nbins_arr, bins_arr);
    }
    catch(std::runtime_error e)
    {
      std::cerr << "ERROR " << e.what();
    }
  }
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  histogram<NDIM,DATA,WEIGHT>::histogram(const histogram_t& H) : base_t(H),
								 _counts (H._counts),
								 _weights (H._weights),
								 _entries (H._entries),
								 _bin_definitions (H._bin_definitions),
								 _flag_sumwts (H._flag_sumwts) { }
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  histogram<NDIM,DATA,WEIGHT>& histogram<NDIM,DATA,WEIGHT>::operator=(const histogram_t& H)
  {
    base_t::operator=(H);
    _counts = H._counts;
    _weights = H._weights;
    _entries = H._entries;
    _bin_definitions = H._bin_definitions;
    _flag_sumwts = H._flag_sumwts;
    return (*this);
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  void histogram<NDIM,DATA,WEIGHT>::_set(::boost::array<unsigned int, NDIM> nbins_arr,
					 ::boost::array<data_t*, NDIM> bins_arr)
  { 
    for(unsigned int idim = 0; idim < NDIM; ++idim)
    {
      for(unsigned int ibin = 0; ibin < nbins_arr[idim]; ++ibin)
      {
	if(*(bins_arr[idim] + ibin) >= *(bins_arr[idim] + ibin + 1))
	{
	  sprintf(::CPP::errors::buffer, "invalid bin initializer in %s @ L%d", __FILE__, __LINE__);
	  throw std::runtime_error(::CPP::errors::buffer);
	}
	_bin_definitions[idim].push_back(typename bin_t::interval_t(*(bins_arr[idim] + ibin),
								    *(bins_arr[idim] + ibin + 1)));
      }
    }
  }
  
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  unsigned long histogram<NDIM,DATA,WEIGHT>::get_num_bins(bool only_filled) const
  {
    unsigned long ntot = 1;
    if(!only_filled)
    {
      for(unsigned int i=0; i < NDIM; ++i)
	ntot *= _bin_definitions[i].size();
      return ntot;
    }
    else
    {
      return _counts.size();
    }
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  int histogram<NDIM,DATA,WEIGHT>::get_ibin(unsigned idim, data_t x) const
  {
    unsigned nbins = _bin_definitions[idim].size();
    
    if(x < _bin_definitions[idim][0].first)
      return -1;
    
    if(x > _bin_definitions[idim][nbins-1].second)
      return nbins;
    
    for(unsigned ibin = 0; ibin < nbins; ++ibin)
    {
      if(_bin_definitions[idim][ibin].first <= x &&
	 x <= _bin_definitions[idim][ibin].second)
	return ibin;
    }
    
    return -1;
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram<NDIM,DATA,WEIGHT>::bin_t 
  histogram<NDIM,DATA,WEIGHT>::get_bin(const binselect_t& bin_def) const
  {
    bool is_underflow = false;
    bool is_overflow  = false;
    typename bin_t::boundary_t edges;
    for(unsigned idim = 0; idim < NDIM; ++idim)
    {
      int j = bin_def[idim]; 
      if(j >= (int)(_bin_definitions[idim].size())) 
      {
	is_overflow = true;
	edges.push_back(typename bin_t::interval_t((DATA)((_bin_definitions[idim][j-1]).second),
						   std::numeric_limits<DATA>::max()));
      }
      if(j < 0)
      {
	is_underflow = true;
	edges.push_back(typename bin_t::interval_t(-std::numeric_limits<DATA>::max(),
						   (DATA)((_bin_definitions[idim][0]).first)));
      }
      if(j >= 0 && j < (int)(_bin_definitions[idim].size()))
      {
	edges.push_back(typename bin_t::interval_t((DATA)((_bin_definitions[idim][j]).first), 
						   (DATA)((_bin_definitions[idim][j]).second))); 
      }
    }
    bin_t the_bin(edges);
    if(is_underflow) the_bin.set_bit((int)base_t::k_underflow, is_underflow);
    if(is_overflow) the_bin.set_bit((int)base_t::k_overflow , is_overflow ); // flag this bin as being under/over-flow 
    return the_bin;
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram<NDIM,DATA,WEIGHT>::bin_t 
  histogram<NDIM,DATA,WEIGHT>::get_bin(const point_t& pt) const
  {
    ::boost::array<int, NDIM> ibin_def;
    unsigned nbins;
    for(unsigned idim = 0; idim < NDIM; ++idim)
    {
      nbins = _bin_definitions[idim].size();
      if(pt[idim] < _bin_definitions[idim][0].first)
      {
	ibin_def[idim] = base_t::IUNDERFLOW;
	continue;
      }
      if(pt[idim] > _bin_definitions[idim][nbins-1].second)
      {
	ibin_def[idim] = nbins;
	continue;      
      }
      for(unsigned ibin = 0; ibin < nbins; ++ibin)
      {
	if(pt[idim] <= _bin_definitions[idim][ibin].second && pt[idim] >= _bin_definitions[idim][ibin].first)
	{
	  ibin_def[idim] = ibin;
	  break;
	}
      }
    }
    return get_bin(binselect_t(ibin_def.data()));
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram<NDIM,DATA,WEIGHT>::entry_t histogram<NDIM,DATA,WEIGHT>::get_bin_content(const bin_t& b) const
  {
    typename map_t::const_iterator itr = _counts.find(b);
    if(itr != _counts.end()) 
    {
      return itr->second; 
    }
    return entry_t(0);
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram<NDIM,DATA,WEIGHT>::entry_t histogram<NDIM,DATA,WEIGHT>::get_bin_error(const bin_t& b) const
  {
    typename map_t::const_iterator itr = _weights.find(b);
    if(itr == _weights.end())
    {
      itr = _counts.find(b);
    }
    else
    {
      if(_flag_sumwts)
      {
	return std::sqrt(itr->second);
      }
      itr = _counts.find(b);
    }
    if(itr == _counts.end())
      return entry_t(0);
    else
      return std::sqrt(itr->second);
    return entry_t(0);
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram<NDIM,DATA,WEIGHT>::entry_t histogram<NDIM,DATA,WEIGHT>::get_integral(bool include_overflow) const
  {
    std::vector<int> buffer;
    for(unsigned int idim = 0; idim < NDIM; ++idim)
      buffer.push_back(_bin_definitions[idim].size() - (include_overflow ? 1 : 0));
    
    return get_integral(get_bin(binselect_t(std::vector<int>(NDIM, (include_overflow ? -1 : 0)))),
			get_bin(binselect_t(buffer)), include_overflow);
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram<NDIM,DATA,WEIGHT>::entry_t histogram<NDIM,DATA,WEIGHT>::get_integral(const bin_t& a, const bin_t& b, bool include_overflow) const
  {
    //
    // to optimize efficiency
    //    - if # of bins in hypercube is smaller than number of filled bins, 
    //      loop over all bins in hypercube
    //    - otherwise loop over all filled bins and sum weights for those 
    //      that are inside the cube
    //
    
    try
    {
      _check_valid_bin(a);
      _check_valid_bin(b);
    }
    catch(std::runtime_error err)
    {
      std::cerr << "ERROR: " << "in get_integral(), bin definitions not valid - use get_bin()" << std::endl;
      return entry_t(0);
    }
    
    unsigned long ncube;
    for(unsigned int idim = 0; idim < NDIM; ++idim)
    {
      ncube *= (abs(get_ibin(idim, (a[idim].first + a[idim].second)/2.0) - 
		    get_ibin(idim, (b[idim].first + b[idim].second)/2.0)));
    }
    
    if(ncube > get_num_bins(true) || true) // default to brute force ... 
    {
      entry_t sumwt;
      
      std::vector<data_t> l, h;
      for(unsigned int idim = 0; idim < NDIM; ++idim)
      {
	l.push_back(a[idim].first);
	h.push_back(b[idim].second);
      }
      bin_t cube(l, h);
      const_filled_bin_iterator itr   = begin();
      const_filled_bin_iterator itr_e = end();
      for(; itr != itr_e; ++itr)
      {
	if(((*itr).first.get_bit(base_t::k_underflow) ||
	    (*itr).first.get_bit(base_t::k_overflow)) && !include_overflow)
	  continue;
	
	if(cube.contains((*itr).first.center()))
	{
	  sumwt += (*itr).second;
	}
      }
      
      return sumwt;
    }
    else
    {  
      // @todo implement this ...
    }
    return entry_t(0);
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  
  // --------------===============--------------===============--------------===============
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  void histogram<NDIM,DATA,WEIGHT>::fill(const point_t& p, const entry_t& w)
  {
    
    // 
    // find the relevant bin edges's in each dimension
    // then retrieve that bin (if it exists, otherwise create it)
    // and fill the bin with weight 'w'
    //
    
    std::vector<int> bin_def(NDIM, 0);
    
    for(unsigned idim=0; idim < NDIM; ++idim)
    {
      unsigned nbins = _bin_definitions[idim].size();
      if(p[idim] < _bin_definitions[idim][0].first) // underflow 
      {
	bin_def[idim] = base_t::IUNDERFLOW;
      }
      if(p[idim] > _bin_definitions[idim][nbins].second) // overflow 
      {
	bin_def[idim] = nbins;
      }
      for(unsigned ib=0; ib < nbins; ++ib)
      {
	if(_bin_definitions[idim][ib].first <= p[idim] && p[idim] <= _bin_definitions[idim][ib].second)
	{
	  bin_def[idim] = ib; // @todo replace with binary search
	  break;
	}
      }
    }
    
    bin_t b = get_bin(binselect_t(bin_def));
    
    if(_counts.find(b) == _counts.end())
    {
      _entries[b] = 0;
      _counts[b] = 0;
    }
    _counts[b] += w;
    _entries[b] += 1;
    
    if(_flag_sumwts)
    {
      if(_weights.find(b) == _weights.end()) 
      {
	_weights[b] = 0;
      }
      _weights[b] += w*w;
    }
    
    this->inc_nentries();
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  void histogram<NDIM,DATA,WEIGHT>::set_bin_content(const bin_t& b, const entry_t& w) 
  {
    _check_valid_bin(b); // throws runtime_error if bin argument is invalid
    if(_flag_sumwts) 
    {
      _weights[b] = w*w;
    }
    _counts[b] = w;
    _entries[b] = 1;
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  void histogram<NDIM,DATA,WEIGHT>::set_bin_error(const bin_t& b, const entry_t& w)
  {
    _check_valid_bin(b); // throws runtime_error if bin argument is invalid
    _weights[b] = w*w;
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  void histogram<NDIM,DATA,WEIGHT>::_check_valid_bin(const bin_t& b) const
  {
    
    //
    // the overflow and underflow bins cannot be set manually, 
    // and the bin definition *must* be a valid combination of edges
    //
    // throw runtime_error if the bin boundary does not correspond to the binning of this histogram
    //
    
    bool found = false;
    
    for(unsigned idim=0; idim < NDIM; ++idim)
    {
      if(b.get_bit(base_t::k_underflow) || b.get_bit(base_t::k_overflow))
	continue;
      
      for(unsigned ibin=0; ibin < _bin_definitions[idim].size(); ++ibin)
      {
	if(_bin_definitions[idim][ibin].first  == b[idim].first &&
	   _bin_definitions[idim][ibin].second == b[idim].second)
	{
	  found = true;
	  break;
	}
      }
      if(!found)
      {
	sprintf(::CPP::errors::buffer, "invalid bin definition in %s @ L%d", __FILE__, __LINE__);
	throw std::runtime_error(::CPP::errors::buffer);
      }
      found = false;
    }
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  bool histogram<NDIM,DATA,WEIGHT>::_check_bin_alignment(const histogram_t& H) const
  {
    
    // 
    // check whether other histogram has identical bin boundaries
    // as this one (necessary condition for histogram addition, division)
    //
    
    unsigned nbins;
    for(unsigned idim = 0; idim < NDIM; ++idim)
    {
      unsigned nbins = _bin_definitions[idim].size();
      if(nbins != H.get_num_bins((unsigned)idim))
	return false;
      for(unsigned ibin = 0; ibin < nbins; ++ibin)
      {
	if(get_bin_low(idim, ibin) != H.get_bin_low(idim, ibin))
	  return false;
	if(get_bin_high(idim, ibin) != H.get_bin_high(idim, ibin))
	  return false;
      }
    }
    
    return true;
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  bool histogram<NDIM,DATA,WEIGHT>::add(const histogram_t* H, entry_t w)
  {
    if(! _check_bin_alignment(H))
    {
      std::cerr << "ERROR: " << "attempted to add histograms with different bin definitions" << std::endl;
      return false;
    }
    bin_t b;
    entry_t e, y;
    const_filled_bin_iterator itr = H->begin();
    const_filled_bin_iterator itr_end = H->end();
    for(; itr != itr_end; ++itr)
    {
      y = itr->second;
      b = itr->first;
      if(H->get_use_sum_weights() && get_use_sum_weights())
	e = H->get_bin_error(b);
      else
	e = 0.;
      set_bin_content(b, get_bin_content(b) + w*y);
      if(get_use_sum_weights())
	set_bin_error(b, std::sqrt(pow(get_bin_error(b), 2) + w*w*e*e));
    }
    return true;
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  void histogram<NDIM,DATA,WEIGHT>::scale(entry_t w)
  {
    bin_t b;
    entry_t e, y;
    
    for(const_filled_bin_iterator itr = begin();
	itr != end(); ++itr)
    {
      b = itr->first;
      y = itr->second;
      if(get_use_sum_weights())
      {
	e = get_bin_error(b);
	set_bin_error(b, w*e);
      }
      set_bin_content(b, w*y);
    }
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  bool histogram<NDIM,DATA,WEIGHT>::rebin(const tuple<NDIM,unsigned int>& factor_tpl) 
  {
    
    //
    // first check that # of bins in each dimension is some
    // multiple of the rebinning factor
    //
    
    for(unsigned idim = 0; idim < NDIM; ++idim)
    {
      if((_bin_definitions[idim].size() % factor_tpl[idim]) != 0)
      {
	std::cerr << "ERROR: " << "attempted to rebin histogram with non-divisible factor : " << idim << " " << factor_tpl[idim] << std::endl;
	return false;
      }
    }
    
    //
    // create new merged bins and copy the counts & weights data
    //
    
    ::boost::array<typename bin_t::boundary_t, NDIM> bin_defs;
    
    std::map<bin_t,entry_t> counts = _counts;
    std::map<bin_t,entry_t> weights = _weights;	
    
    for(unsigned idim = 0; idim < NDIM; ++idim)
    {
      unsigned jbin = 0;
      data_t a = _bin_definitions[idim][0].first;
      data_t b = _bin_definitions[idim][0].second;
      for(unsigned ibin = 1; ibin < _bin_definitions[idim].size(); ++ibin)
      {
	if((ibin+1) % factor_tpl[idim] == 0)
	{
	  b = _bin_definitions[idim][ibin].second;
	  bin_defs[idim].push_back(typename bin_t::interval_t(a, b));
	  a = b;
	}
      }
    }
    
    //
    // reset the bin definitions and clear the old counts & weights data
    // then set the bin contents of the new bins
    //
    
    _bin_definitions = bin_defs;
    
    clear();
    
    const_filled_bin_iterator itr = counts.begin();
    for(; itr != counts.begin(); ++itr)
    {
      bin_t bin_def = get_bin((*itr).first.center());
      entry_t y = get_bin_content(bin_def) + (*itr).second;
      set_bin_content(bin_def, y);
      if(get_use_sum_weights())
      {
	const_filled_bin_iterator itr_w = weights.find((*itr).first);
	if(itr_w != weights.end())
	{
	  entry_t w = get_bin_error(bin_def) + itr_w->second;
	  set_bin_error(bin_def, w);
	}
      }
    }
    
    return true;
  }

  // --------------===============--------------===============--------------===============

  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  histogram<1,DATA,WEIGHT> histogram<NDIM,DATA,WEIGHT>::unroll(data_t xlow, data_t xhigh, const std::vector< std::vector<bin_t> >& mapping) const
  {
    typedef typename histogram<1,DATA,WEIGHT>::binselect_t binselect1D_t;
    typedef typename histogram<1,DATA,WEIGHT>::bin_t bin1D_t;
    data_t xpt = xlow;    
    ::boost::array<unsigned int, 1> nbins_arr = { { _entries.size() } };
    ::boost::array<data_t*, 1> bins_arr = { { new data_t[_entries.size()+1] } };
    
    for(unsigned int ib=0; ib <= _entries.size(); ++ib)
    {
      bins_arr[0][ib] = xlow + ib*((xhigh - xlow)/_entries.size());
    }
    histogram<1,DATA,WEIGHT> h(nbins_arr, bins_arr);
    const_filled_bin_iterator citr = begin();
    const_filled_bin_iterator eitr = err_begin();
    long nentries = get_nentries();
    while(citr != end())
    {
      // use bin mapping to find output index for each bin
      bool mapped = false;
      for(unsigned int ib=0; ib < mapping.size(); ++ib)
      {
	for(unsigned int jb=0; jb < mapping[ib].size(); ++jb)
	{
	  if(mapping[ib][jb].contains(citr->first))
	  {
	    bin1D_t b = h.get_bin(binselect1D_t(ib));
            h.set_bin_content(b, h.get_bin_content(b) + citr->second);
	    h.set_bin_error(b, std::sqrt(pow(h.get_bin_error(b),2) + eitr->second));
	    mapped = true;
	    break;
	  }
	}
	if(mapped) break;
      }
      if(!mapped)
      {
	// exclude these entries from total
	nentries -= _entries[citr->first];
      }
      ++citr;
      ++eitr;
    }
    h.set_nentries(nentries);
    return h;
  }

  // --------------===============--------------===============--------------===============

  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  histogram<1,DATA,WEIGHT> histogram<NDIM,DATA,WEIGHT>::unroll(data_t xlow, data_t xhigh, const std::vector< std::vector<binselect_t> >& mapping) const
  {
    typedef typename histogram<1,DATA,WEIGHT>::binselect_t binselect1D_t;
    typedef typename histogram<1,DATA,WEIGHT>::bin_t bin1D_t;
    data_t xpt = xlow;    
    ::boost::array<unsigned int, 1> nbins_arr = { { _entries.size() } };
    ::boost::array<data_t*, 1> bins_arr = { { new data_t[_entries.size()+1] } };
    
    for(unsigned int ib=0; ib <= _entries.size(); ++ib)
    {
      bins_arr[0][ib] = xlow + ib*((xhigh - xlow)/_entries.size());
    }
    histogram<1,DATA,WEIGHT> h(nbins_arr, bins_arr);

    long nentries = 0;
    
    for(unsigned int ib=0; ib < mapping.size(); ++ib)
    {
      for(unsigned int jb=0; jb < mapping[ib].size(); ++jb)
      {
	bin_t src_b = get_bin(mapping[ib][jb]);
	bin1D_t dest_b = h.get_bin(binselect1D_t(ib));
	h.set_bin_content(dest_b, h.get_bin_content(dest_b) + _counts[src_b]);
	h.set_bin_error(dest_b, std::sqrt(pow(h.get_bin_error(dest_b),2) + _weights[src_b]));
	nentries += _entries[src_b];
      }
    }
    h.set_nentries(nentries);
    return h;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  //
  // here are some hacks to create factories for histograms of dimension 1,...,10
  //
  ///////////////////////////////////////////////////////////////////////////////////////////
  
  namespace histfactory 
  {
    
#ifndef __HAVE_C11_STD
    
#define ARRY_CARG(z,n,nil) darr_##n,
#define ARRY_NARG(z,n,nil) N##n,
#define HIST_CARG_A(z,n,nil) unsigned N##n, const DATA& dl_##n, const DATA& du_##n,	
#define HIST_CTOR_IMP_A(z,n,nil)					\
    template<unsigned int NDIM,						\
	     typename DATA>						\
    std::auto_ptr<histogram<NDIM,DATA,double> >				\
    create(BOOST_PP_REPEAT(n, HIST_CARG_A, 0)				\
	   unsigned N##n, const DATA& dl_##n, const DATA& du_##n)	\
    {									\
      ::boost::array<unsigned, NDIM> nbins =				\
	{ BOOST_PP_REPEAT(n, ARRY_NARG, 0) N##n };			\
      ::boost::array<DATA*, NDIM> edges = { 0x0 };			\
      for(unsigned int idim = 0; idim < NDIM; ++idim)			\
      {									\
	edges[idim] = new DATA[N##n+1];					\
	for(unsigned int iedge = 0; iedge < N##n; ++iedge)		\
	{								\
	  edges[idim][iedge] = (du_##n - dl_##n)/N##n * iedge + dl_##n;	\
	}								\
	edges[idim][N##n] = du_##n;					\
      }									\
      histogram<NDIM,DATA,double>* p =					\
	new histogram<NDIM,DATA,double>(nbins, edges);			\
      for(unsigned int idim = 0; idim < n; ++idim)			\
      {									\
	delete[] (edges[idim]);						\
      }									\
      return std::auto_ptr<histogram<NDIM,DATA,double> >(p);		\
    }
    
    BOOST_PP_REPEAT(TUPLE_MAX_NDIMENSIONS, HIST_CTOR_IMP_A, 0)
    
#undef ARRY_CARG
#undef ARRY_NARG
#undef HIST_CARG_A
#undef HIST_CTOR_IMP_A
    
#define ARRY_CARG(z,n,nil) darr_##n,
#define ARRY_NARG(z,n,nil) N##n,
#define HIST_CARG_B(z,n,nil) unsigned N##n, const DATA* darr_##n,	
#define HIST_CTOR_IMP_B(z,n,nil)					\
    template<unsigned int NDIM,						\
	     typename DATA>						\
    std::auto_ptr<histogram<NDIM,DATA,double> >				\
    create(BOOST_PP_REPEAT(n, HIST_CARG_B, 0)				\
	   unsigned N##n, const DATA* darr_##n)				\
    {									\
      ::boost::array<unsigned, NDIM> nbins =				\
	{ BOOST_PP_REPEAT(n, ARRY_NARG, 0) N##n };			\
      ::boost::array<DATA*, NDIM> edges    =				\
	  { BOOST_PP_REPEAT(n, ARRY_CARG, 0) darr_#n };			\
      histogram<NDIM,DATA,double>* p =					\
	new histogram<NDIM,DATA,double>(nbins, edges);			\
      return std::auto_ptr<histogram<NDIM,DATA,double> >(p);		\
    }
  }
  
#undef ARRY_CARG
#undef ARRY_NARG
#undef HIST_CARG_B
#undef HIST_CTOR_IMP_B
  
#define DECLARE_HIST_TYPEDEFS(z,n,nil)			\
  typedef histogram<n,float,float> h##n##ff;		\
  typedef std::auto_ptr<h##n##ff> h##n##ff_ptr;		\
  typedef histogram<n,float,float> h##n##f;		\
  typedef std::auto_ptr<h##n##ff> h##n##f_ptr;		\
  typedef histogram<n,float,double> h##n##fd;		\
  typedef std::auto_ptr<h##n##fd> h##n##fd_ptr;		\
  typedef histogram<n,double,float> h##n##df;		\
  typedef std::auto_ptr<h##n##df> h##n##df_ptr;		\
  typedef histogram<n,double,double> h##n##dd;		\
  typedef std::auto_ptr<h##n##dd> h##n##dd_ptr;		\
  typedef histogram<n,double,double> h##n##d;		\
  typedef std::auto_ptr<h##n##dd> h##n##d_ptr;	
  
  BOOST_PP_REPEAT(TUPLE_MAX_NDIMENSIONS, DECLARE_HIST_TYPEDEFS, 0)
  
#undef DECLARE_HIST_TYPEDEFS
  
  // example usage:
  //
  // ::histfactory::h1d_ptr = ::histfactory::create<1,double>(10, 0, 1);
  //
  // ::histfactory::h5d_ptr = ::histfactory::create<5,double>(10, 0, 1, ..., 10, -1, 0);
  // 
  
#else
  
  // @todo
  
#endif
  
}

///////////////////////////////////////////////////////////////////////////////////////////

#endif
