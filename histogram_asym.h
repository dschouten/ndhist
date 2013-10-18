// @author doug schouten <doug dot schouten at triumf dot ca>
// @file histogram.h

#ifndef HISTOGRAM_ASYM_H
#define HISTOGRAM_ASYM_H

#include <iostream>
#include <cmath>
#include <map>
#include <memory>
#include <vector>
#include <algorithm>
#include <boost/array.hpp>

/*************************************************************************************
 *
 * histogram_asym<N,D,W> - this histogram type allows for binned volumes that 
 * are not rectangular. Graphically, histograms generally allow only binning like
 *
 * -----------------------
 * |    |        |   |   |
 * -----------------------
 * |    |        |   |   |
 * -----------------------
 * |    |        |   |   |
 * -----------------------
 * |    |        |   |   |
 * -----------------------
 *
 * in eg., two dimensions, but with this histogram type, one can build histograms like
 *
 * ------------------
 * |    |        |   |
 * -------------------
 *      |    |   |   |
 *      ------------------
 *      |    |   |   |   |
 * -----------------------
 * |    |        |   |   |
 * ------        ---------
 *
 * the draw-back is that bins are defined individually (although bin
 * boundary checking is performed). Furthermore, once a bin is defined
 * it must be stored, so the memory resources needed increases linearly 
 * with the number of bins (in contrast to a histogram where only the edges
 * along each dimension need to be stored)
 *
 * Any histogram_asym can be written out in a 1D unpacked histogram
 * where the ordering along the axis is according to the order the bins were
 * defined
 *
 * histogram_asym<N,D,W>
 * 
 * N is the number of dimensions, in which the data is of type D and 
 * the count weights are of type W
 *
 *************************************************************************************/

#include <bin.h>
#include <tuple.h>

#include <histogram_base.h>
#include <histogram.h>

namespace hepstd
{
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  class histogram_asym : public histogram_base<NDIM, DATA, WEIGHT>
  {
  public:
    typedef histogram_base<NDIM,DATA,WEIGHT> base_t;  
    typedef histogram_asym<NDIM,DATA,WEIGHT> histogram_t;
    typedef histogram<1,DATA,WEIGHT> histogram1D_t;

    typedef typename base_t::entry_t      entry_t;
    typedef typename base_t::data_t       data_t;
    typedef typename base_t::bin_t        bin_t;
    typedef typename base_t::point_t      point_t;
    typedef typename base_t::binselect_t  binselect_t;
    
  public:
    histogram_asym();
    histogram_asym(const std::vector<bin_t>&);
    
    histogram_asym(const histogram_t&); // << copy constructor
    histogram_asym& operator=(const histogram_t&); // << assignment operator
    
    histogram_asym(const base_t&); // << build from filled bins of histogram
    
    virtual ~histogram_asym() { }
    
  public:
    
    //
    // add a new bin definition (return 0 if failed, otherwise returns number of bins after adding)
    //
    
    unsigned int add_bin(const bin_t&);
    
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
    // allow to iterate over bins
    //    
    
    typedef typename std::vector<entry_t>::const_iterator const_bin_iterator;
    
    const_bin_iterator begin() const { return  _counts.begin(); }
    const_bin_iterator   end() const { return  _counts.end(); }
    
    const_bin_iterator err_begin() const { return  _weights.begin(); }
    const_bin_iterator   err_end() const { return  _weights.end(); }
    
  public:
    
    bool get_use_sum_weights     () const { return _flag_sumwts; }
    void set_use_sum_weights     (bool flag = true) { _flag_sumwts = flag; if(!flag) _weights.clear(); } 
     
    // bool add			 (const histogram_t*, entry_t w = (entry_t)1.); // add histogram to self with weight factor
    // void scale			 (entry_t); // scale self
    // bool rebin			 (const tuple<NDIM,unsigned int>&); // re-define bin edges by some factor along each dimension
    // void clear                   (); // clear all histogrammed data
    // void reset                   (); // reset the histogram (remove all bin definitions)

    histogram1D_t unroll (data_t xlow, data_t xhigh, bool uniform=true) const; // return 1D histogram (default with uniform bins)
    
    //
    // define histogram algebra 
    //
    
    // histogram_t& operator*=(entry_t S)		        { this->scale(S); return (*this); }
    // histogram_t& operator+=(const histogram_t& H)	{ this->add(&H); return (*this); }
    
    // histogram_t operator+(const histogram_t& H) const	{ histogram_t B = (*this); B.add(&H,+1); return B; }
    // histogram_t operator-(const histogram_t& H) const	{ histogram_t B = (*this); B.add(&H,-1); return B; }
    // histogram_t operator*(entry_t S) const		{ histogram_t B = (*this); B.scale(S); return B; }
    
  private:
    
    std::map<bin_t,long> _bin_indices;
    std::vector<entry_t> _counts;
    std::vector<entry_t> _weights;
    std::vector<long> _entries;
    std::vector<bin_t> _bins;

    ::boost::array<data_t,NDIM> _xmin;
    ::boost::array<data_t,NDIM> _xmax;
        
    bool _flag_sumwts;
    
    void _check_valid_bin(const bin_t&) const; // throws std::runtime_error
    bool _check_overlap(const bin_t&) const; // check if bin overlaps with any bin in histogram (note complexity is O(nbins))
    bool _check_bin_alignment(const histogram_t&) const; // check if bin definitions for histograms are aligned (note complexity is O(nbins))

    class sort_along_dim
    {
    public:
      sort_along_dim(unsigned int idim) { _idim = idim; }

      bool operator()(const bin_t& x, const bin_t& y) const 
      { return x[_idim].second <= y[_idim].first; }

    private:
      sort_along_dim();
      unsigned int _idim;
    };

  };
  
  
  // --------------===============--------------===============--------------===============

  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  histogram_asym<NDIM,DATA,WEIGHT>::histogram_asym() : base_t()
  {
    std::fill(_xmin.begin(), _xmin.end(), 0);
    std::fill(_xmax.begin(), _xmax.end(), 0);
  }    

  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  histogram_asym<NDIM,DATA,WEIGHT>::histogram_asym(const std::vector<bin_t>& bins) : base_t(),
										     _flag_sumwts(true)
  { 
    std::fill(_xmin.begin(), _xmin.end(), 0);
    std::fill(_xmax.begin(), _xmax.end(), 0);
    for(unsigned int ib = 0; ib < bins.size(); ++ib)
    {
      add_bin(bins[ib]);
    }
  }
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  histogram_asym<NDIM,DATA,WEIGHT>::histogram_asym(const histogram_t& H) : base_t(),
									   _bin_indices (H._bin_indices),
									   _counts (H._counts),
									   _weights (H._weights),
									   _entries (H._entries),
									   _bins (H._bins),
									   _xmin (H._xmin),
									   _xmax (H._xmax),
									   _flag_sumwts (H._flag_sumwts) { }
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  histogram_asym<NDIM,DATA,WEIGHT>& histogram_asym<NDIM,DATA,WEIGHT>::operator=(const histogram_t& H)
  {
    base_t::operator=(H);
    _bin_indices = H._bin_indices;
    _counts = H._counts;
    _weights = H._weights;
    _entries = H._entries;
    _bins = H._bins;
    _xmin = H._xmin;
    _xmax = H._xmax;
    _flag_sumwts = H._flag_sumwts;
    return (*this);
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////

  // --------------===============--------------===============--------------===============

  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  unsigned int histogram_asym<NDIM,DATA,WEIGHT>::add_bin(const bin_t& b)
  {
    if(!_check_overlap(b))
    {
      return 0;
    }
    _bins.push_back(b);
    unsigned int ib = _bins.size()-1;
    _bin_indices[b] = ib;
    _counts[ib] = entry_t(0);
    if(_flag_sumwts)
    {
      _weights[ib] = entry_t(0);
    }
    _entries[ib] = 0;
    for(unsigned int idim=0; idim < NDIM; ++idim)
    {
      if(b[idim] < _xmin[idim]) _xmin[idim] = b[idim];
      if(b[idim] > _xmax[idim]) _xmax[idim] = b[idim];
    }
    return _bins.size();
  }

  // --------------===============--------------===============--------------===============

  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  void histogram_asym<NDIM,DATA,WEIGHT>::fill(const point_t& p, const entry_t& w)
  {
    for(unsigned int ib=0; ib < _bins.size(); ++ib)
    {
      if(_bins[ib].contains(p))
      {
	_entries[ib] += 1;
	if(_flag_sumwts)
	{
	  _weights[ib] += w*w;
	}
	_counts[ib] += w;
      }
    }
  }

  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram_asym<NDIM,DATA,WEIGHT>::bin_t histogram_asym<NDIM,DATA,WEIGHT>::get_bin(const binselect_t& x) const
  {
    // order bins along each dimension, then find the corresponding index

    bool is_underflow = false;
    bool is_overflow  = false;
    bool is_undefined = false;
    typename bin_t::boundary_t edges;
    for(unsigned idim = 0; idim < NDIM; ++idim)
    {
      int j = x[idim];
      if(j >= _bin_indices.size())
      {
	is_overflow = true;
	edges.push_back(typename bin_t::interval_t((DATA)(_xmax[idim]), +std::numeric_limits<DATA>::max()));
      }
      if(j == this->IUNDERFLOW)
      {	
	is_underflow = true;
	edges.push_back(typename bin_t::interval_t(-std::numeric_limits<DATA>::max(), (DATA)(_xmin[idim])));
      }
      if(j == this->IUNDEFINED)
      {
	is_undefined = true;
	edges.push_back(typename bin_t::interval_t(-std::numeric_limits<DATA>::max(), +std::numeric_limits<DATA>::max()));
      }
      if(j >= 0 && j < (int)(_bin_indices.size()))
      {
	std::vector<bin_t> loc(_bins);
	std::sort(loc.begin(), loc.end(), sort_along_dim(idim));
	edges.push_back(loc[j].first);
      }
    }
    bin_t the_bin(edges);
    if(is_undefined) the_bin.set_bit((int)this->k_undefined, is_undefined);
    if(is_underflow) the_bin.set_bit((int)this->k_underflow, is_underflow);
    if(is_overflow) the_bin.set_bit((int)this->k_overflow , is_overflow ); // flag this bin as being under/over-flow 
    return the_bin;    
  }

  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram_asym<NDIM,DATA,WEIGHT>::bin_t histogram_asym<NDIM,DATA,WEIGHT>::get_bin(const point_t& x) const
  {
    ::boost::array<int,NDIM> ibin_def;
    for(unsigned idim = 0; idim < NDIM; ++idim)
    {
      if(x[idim] < _xmin[idim] || x[idim] > _xmax[idim])
      {
	return bin_t(); // @fixme use underflow/overflow
      }
      for(unsigned ibin = 0; ibin < _bins.size(); ++ibin)
      {
	ibin_def[idim] = this->IUNDEFINED;
	if(x[idim] <= _bins[ibin].second and x[idim] >= _bins[ibin].first)
	{
	  ibin_def[idim] = ibin;
	  break;
	}
	if(ibin_def[idim] == this->IUNDEFINED)
	  return bin_t(); // @fixme use undefined
      }
    }
    return get_bin(binselect_t(ibin_def.data()));
  }
  
  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram_asym<NDIM,DATA,WEIGHT>::entry_t histogram_asym<NDIM,DATA,WEIGHT>::get_bin_content(const bin_t& b) const 
  {
    if( _bin_indices.find(b) != _bin_indices.end() )
    {
      return _counts[_bin_indices[b]];
    }
    return entry_t(0);
  }
  
  // --------------===============--------------===============--------------===============
    
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram_asym<NDIM,DATA,WEIGHT>::entry_t histogram_asym<NDIM,DATA,WEIGHT>::get_bin_error(const bin_t& b) const
  {
    if( _bin_indices.find(b) != _bin_indices.end() )
    {
      if(_flag_sumwts) return std::sqrt(_weights[_bin_indices[b]]);
      else return std::sqrt(_counts[_bin_indices[b]]);
    }
    return entry_t(0);
  }
  
  // --------------===============--------------===============--------------===============

  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram_asym<NDIM,DATA,WEIGHT>::entry_t histogram_asym<NDIM,DATA,WEIGHT>::get_integral(const bin_t&, const bin_t&, bool include_overflow) const
  {
    return entry_t(0); // @ fix me
  }

  // --------------===============--------------===============--------------===============

  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram_asym<NDIM,DATA,WEIGHT>::entry_t histogram_asym<NDIM,DATA,WEIGHT>::get_integral(bool include_overflow) const
  {
    return entry_t(0); // @ fix me
  }

  // --------------===============--------------===============--------------===============

  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  void histogram_asym<NDIM,DATA,WEIGHT>::set_bin_content(const bin_t&, const entry_t&) 
  {} // @ fix me

  // --------------===============--------------===============--------------===============
  
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  void histogram_asym<NDIM,DATA,WEIGHT>::set_bin_error(const bin_t&, const entry_t&) 
  {} // @ fix me

  // --------------===============--------------===============--------------===============

  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  typename histogram_asym<NDIM,DATA,WEIGHT>::histogram1D_t histogram_asym<NDIM,DATA,WEIGHT>::unroll(data_t xlow, data_t xhigh, bool uniform) const
  {
    data_t vol = 0; 
    data_t xpt = xlow;    
    ::boost::array<unsigned int, 1> nbins_arr = { { _bins.size() } };
    ::boost::array<data_t*, 1> bins_arr = { { new data_t[_bins.size()+1] } };
    if(!uniform)
    {
      BOOST_FOREACH(bin_t b, _bins) { vol += b.volume(); }
    }
    for(unsigned int ib=0; ib <= _bins.size(); ++ib)
    {
      if(uniform)
      {
	bins_arr[0][ib] = xlow + ib*((xhigh - xlow)/_bins.size());
      }
      else 
      {
	bins_arr[0][ib] = xpt;
	if(ib>0) 
	  bins_arr[0][ib] += (xhigh - xlow)*_bins[ib-1].volume()/vol;
	xpt = bins_arr[0][ib];
      }
    }
    histogram1D_t h(nbins_arr, bins_arr);
    for(unsigned int ib=0; ib < _bins.size(); ++ib)
    {
      h.set_bin_content(h.get_bin(histogram1D_t::binselect_t(ib)), _counts[ib]);
      h.set_bin_error(h.get_bin(histogram1D_t::binselect_t(ib)), std::sqrt(_weights[ib]));
    }
    h.set_entries(this->_nentries);
    return h;
  }

}

#endif
