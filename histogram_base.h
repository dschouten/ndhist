#ifndef HISTOGRAM_BASE
#define HISTOGRAM_BASE

namespace hepstd
{
  template<unsigned int NDIM, typename DATA, typename WEIGHT>
  class histogram_base
  {
  public:
    enum bin_status { k_underflow = 1, 
		      k_overflow  = 2,
                      k_undefined = 3 };
    
    static const int IUNDERFLOW      = -1;
    static const int IUNDEFINED      = -2;
    static const unsigned NDIMENSION = NDIM;
    
    typedef WEIGHT                            entry_t;
    typedef DATA                              data_t;
    typedef bin<NDIM,DATA>		      bin_t;
    typedef tuple<NDIM,DATA>		      point_t;
    typedef tuple<NDIM,int>		      binselect_t;
    
  public:

    histogram_base() : _nentries(0) {}
    histogram_base(const histogram_base& H) { this->_nentries = H.get_nentries(); }

    histogram_base& operator=(const histogram_base& H) { this->_nentries = H.get_nentries(); return (*this); }

    //
    // API
    //

    virtual void fill(const point_t&, const entry_t& w = 1) = 0;
    
    virtual bin_t get_bin		(const binselect_t&) const = 0; // find bin by index
    virtual bin_t get_bin               (const point_t&) const = 0; // find bin containing N-dimensional point
    
    virtual entry_t get_bin_content	(const bin_t&) const = 0;
    virtual entry_t get_bin_error	(const bin_t&) const = 0;
    
    virtual entry_t get_integral        (const bin_t&, const bin_t&, bool include_overflow = false) const = 0; // calculate integral in hypercube define by two corners
    virtual entry_t get_integral        (bool include_overlow = false) const = 0; // calculate total integral  
        
    virtual void set_bin_content	(const bin_t&, const entry_t&) = 0; 
    virtual void set_bin_error		(const bin_t&, const entry_t&) = 0;    

    virtual long get_nentries           () const { return _nentries; }    
    virtual void set_nentries           (long n) { _nentries = n; } 
    virtual void inc_nentries           () { _nentries += 1; }

  protected:
    long _nentries;
  };
}

#endif
