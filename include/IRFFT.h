// -*- Mode: c++; tab-width: 8; c-basic-offset: 2; indent-tabs-mode: t -*-
// NOTE: the first line of this file sets up source code indentation rules
// for Emacs; it is also a hint to anyone modifying this file.

/*
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


// File         : fft.hxx
// Author       : Pavel A. Koshevoy
// Created      : 2005/06/03 10:16
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : Wrapper class and helper functions for working with
//                FFTW3 and ITK images.

#ifndef FFT_HXX_
#define FFT_HXX_

// FFTW includes:
// #include <fftw3.h>

// ITK includes:
#include <itkImage.h>
#include <itkComplexToComplexFFTImageFilter.h>

// system includes:
#include <complex>
#include <assert.h>
#include <stdlib.h>


namespace itk_fft
{
  using TransformDirectionEnum = itk::ComplexToComplexFFTImageFilterEnums::TransformDirection;
  //----------------------------------------------------------------
  // set_num_fftw_threads
  //
  // set's number of threads used by fftw, returns previous value:
  // extern std::size_t set_num_fftw_threads(std::size_t num_threads);
  
  //----------------------------------------------------------------
  // itk_image_t
  // 
  typedef itk::Image<float, 2> itk_image_t;
  
  //----------------------------------------------------------------
  // itk_imageptr_t
  // 
  typedef itk_image_t::Pointer itk_imageptr_t;
  
  //----------------------------------------------------------------
  // fft_complex_t
  // 
  typedef std::complex<float> fft_complex_t;
  
  //----------------------------------------------------------------
  // itk_complex_image_t
  // 
  typedef itk::Image<fft_complex_t, 2> itk_complex_image_t;
  
  //----------------------------------------------------------------
  // itk_complex_imageptr_t
  // 
  typedef itk_complex_image_t::Pointer itk_complex_imageptr_t;

  //----------------------------------------------------------------
  // fft_data_t
  // 
  class fft_data_t
  {
  public:
    fft_data_t(): image_(), nx_(0), ny_(0) {}
    fft_data_t(const unsigned int w, const unsigned int h);
    explicit fft_data_t(const itk_imageptr_t & real);
    fft_data_t(const itk_imageptr_t & real,
	       const itk_imageptr_t & imag);
    fft_data_t(const fft_data_t & data);
    ~fft_data_t() { cleanup(); }
    
    fft_data_t & operator = (const fft_data_t & data);
    
    void cleanup();
    
    void resize(const unsigned int w, const unsigned int h);
    
    void fill(const float real, const float imag = 0.0);
    
    void setup(const itk_imageptr_t & real,
	       const itk_imageptr_t & imag = itk_imageptr_t());
    
    // ITK helpers:
    itk_imageptr_t component(const bool imag = 0) const;
    inline itk_imageptr_t real() const { return component(0); }
    inline itk_imageptr_t imag() const { return component(1); }
    
    // Apply a low-pass filter to this image. This function will
    // zero-out high-frequency components, where the cutoff frequency
    // is specified by radius r in [0, 1]. The sharpness of the
    // cutoff may be controlled by parameter s, where s == 0 results in
    // an ideal low-pass filter, and s == 1 is a low pass filter defined
    // by a scaled and shifted cosine function, 1 at the origin,
    // 0.5 at the cutoff frequency and 0 at twice the cutoff frequency.
    void apply_lp_filter(const double r, const double s = 0);
    
    // accessors:
    inline const itk_complex_image_t * data() const { return image_.get(); }
    inline       itk_complex_image_t * data()       { return image_.get(); }
    
    inline unsigned int nx() const { return nx_; }
    inline unsigned int ny() const { return ny_; }
    
    inline const fft_complex_t & operator() (const unsigned int & x,
					     const unsigned int & y) const
    { return at(x, y); }
    
    inline fft_complex_t & operator() (const unsigned int & x,
				       const unsigned int & y)
    { return at(x, y); }
    
    inline const fft_complex_t & at(const unsigned int & x,
				    const unsigned int & y) const
    { return image_->GetPixel({ x, y }); }
    
    inline fft_complex_t & at(const unsigned int & x,
			      const unsigned int & y)
    { return image_->GetPixel({ x, y }); }
    
  protected:
    itk_complex_imageptr_t image_;
    unsigned int nx_;
    unsigned int ny_;
  };
  
  
  //----------------------------------------------------------------
  // fft
  // 
  extern bool
  fft(itk_image_t::ConstPointer & in, fft_data_t & out);
  
  //----------------------------------------------------------------
  // fft
  // 
  extern bool
  fft(const fft_data_t & in, fft_data_t & out, TransformDirectionEnum sign = TransformDirectionEnum::FORWARD);
  
  //----------------------------------------------------------------
  // ifft
  // 
  extern bool
  ifft(const fft_data_t & in, fft_data_t & out);
  
  //----------------------------------------------------------------
  // ifft
  // 
  inline fft_data_t
  ifft(const fft_data_t & in)
  {
    fft_data_t out;
#ifndef NDEBUG // get around an annoying compiler warning:
    bool ok =
#endif
      ifft(in, out);
    assert(ok);
    return out;
  }
  
  
  // //----------------------------------------------------------------
  // // fn_fft_c_t
  // // 
  // typedef fft_complex_t(*fn_fft_c_t)(const fft_complex_t & in);
  
  // //----------------------------------------------------------------
  // // fn_fft_cc_t
  // // 
  // typedef fft_complex_t(*fn_fft_cc_t)(const fft_complex_t & a,
	// 			      const fft_complex_t & b);
  
  // //----------------------------------------------------------------
  // // elem_by_elem
  // // 
  // extern void
  // elem_by_elem(fn_fft_c_t f,
	//        const fft_data_t & in,
	//        fft_data_t & out);
  
  // //----------------------------------------------------------------
  // // elem_by_elem
  // // 
  // extern void
  // elem_by_elem(fft_complex_t(*f)(const float & a,
	// 			 const fft_complex_t & b),
	//        const float & a,
	//        const fft_data_t & b,
	//        fft_data_t & c);
  
  // //----------------------------------------------------------------
  // // elem_by_elem
  // // 
  // extern void
  // elem_by_elem(fft_complex_t(*f)(const fft_complex_t & a,
	// 			 const float & b),
	//        const fft_data_t & a,
	//        const float & b,
	//        fft_data_t & c);
  
  // //----------------------------------------------------------------
  // // elem_by_elem
  // // 
  // extern void
  // elem_by_elem(fft_complex_t(*f)(const fft_complex_t & a,
	// 			 const fft_complex_t & b),
	//        const fft_complex_t & a,
	//        const fft_data_t & b,
	//        fft_data_t & c);
  
  // //----------------------------------------------------------------
  // // elem_by_elem
  // // 
  // extern void
  // elem_by_elem(fft_complex_t(*f)(const fft_complex_t & a,
	// 			 const fft_complex_t & b),
	//        const fft_data_t & a,
	//        const fft_complex_t & b,
	//        fft_data_t & c);
  
  // //----------------------------------------------------------------
  // // elem_by_elem
  // // 
  // extern void
  // elem_by_elem(fn_fft_cc_t f,
	//        const fft_data_t & a,
	//        const fft_data_t & b,
	//        fft_data_t & c);
  
  
  // //----------------------------------------------------------------
  // // conj
  // // 
  // // element-by-element complex conjugate:
  // // 
  // inline void
  // conj(const fft_data_t & in, fft_data_t & out)
  // { elem_by_elem(std::conj, in, out); }
  
  // //----------------------------------------------------------------
  // // conj
  // // 
  // inline fft_data_t
  // conj(const fft_data_t & in)
  // {
  //   fft_data_t out;
  //   conj(in, out);
  //   return out;
  // }
  
  
  // //----------------------------------------------------------------
  // // sqrt
  // // 
  // // element-by-element square root:
  // // 
  // inline void
  // sqrt(const fft_data_t & in, fft_data_t & out)
  // { elem_by_elem(std::sqrt, in, out); }
  
  // //----------------------------------------------------------------
  // // sqrt
  // // 
  // inline fft_data_t
  // sqrt(const fft_data_t & in)
  // {
  //   fft_data_t out;
  //   sqrt(in, out);
  //   return out;
  // }
  
  
  // //----------------------------------------------------------------
  // // _mul
  // //
  // template <class a_t, class b_t>
  // inline fft_complex_t
  // _mul(const a_t & a, const b_t & b)
  // { return a * b; }
  
  // //----------------------------------------------------------------
  // // mul
  // // 
  // // element-by-element multiplication, c = a * b:
  // // 
  // template <class a_t, class b_t>
  // inline void
  // mul(const a_t & a, const b_t & b, fft_data_t & c)
  // { elem_by_elem(_mul, a, b, c); }
  
  // //----------------------------------------------------------------
  // // mul
  // // 
  // template <class a_t, class b_t>
  // inline fft_data_t
  // mul(const a_t & a, const b_t & b)
  // {
  //   fft_data_t c;
  //   mul<a_t, b_t>(a, b, c);
  //   return c;
  // }
  
  
  //----------------------------------------------------------------
  // _div
  // 
  template <class a_t, class b_t>
  inline fft_complex_t
  _div(const a_t & a, const b_t & b)
  { return a / b; }
  
  // //----------------------------------------------------------------
  // // div
  // // 
  // // element-by-element division, c = a / b:
  // // 
  // template <class a_t, class b_t>
  // inline void
  // div(const a_t & a, const b_t & b, fft_data_t & c)
  // { elem_by_elem(_div, a, b, c); }
  
  // //----------------------------------------------------------------
  // // div
  // // 
  // template <class a_t, class b_t>
  // inline fft_data_t
  // div(const a_t & a, const b_t & b)
  // {
  //   fft_data_t c;
  //   div<a_t, b_t>(a, b, c);
  //   return c;
  // }
  
  
  //----------------------------------------------------------------
  // _add
  // 
  template <class a_t, class b_t>
  inline fft_complex_t
  _add(const a_t & a, const b_t & b)
  { return a + b; }
  
  // //----------------------------------------------------------------
  // // add
  // // 
  // // element-by-element addition, c = a + b:
  // // 
  // template <class a_t, class b_t>
  // inline void
  // add(const a_t & a, const b_t & b, fft_data_t & c)
  // { elem_by_elem(_add, a, b, c); }
  
  // //----------------------------------------------------------------
  // // add
  // // 
  // template <class a_t, class b_t>
  // inline fft_data_t
  // add(const a_t & a, const b_t & b)
  // {
  //   fft_data_t c;
  //   add<a_t, b_t>(a, b, c);
  //   return c;
  // }
  
  
  // //----------------------------------------------------------------
  // // _sub
  // // 
  // template <class a_t, class b_t>
  // inline fft_complex_t
  // _sub(const a_t & a, const b_t & b)
  // { return a - b; }
  
  // //----------------------------------------------------------------
  // // sub
  // // 
  // // element-by-element subtraction, c = a - b:
  // // 
  // template <class a_t, class b_t>
  // inline void
  // sub(const a_t & a, const b_t & b, fft_data_t & c)
  // { elem_by_elem(_sub, a, b, c); }
  
  // //----------------------------------------------------------------
  // // sub
  // // 
  // template <class a_t, class b_t>
  // inline fft_data_t
  // sub(const a_t & a, const b_t & b)
  // {
  //   fft_data_t c;
  //   sub<a_t, b_t>(a, b, c);
  //   return c;
  // }
}


#endif // FFT_HXX_
