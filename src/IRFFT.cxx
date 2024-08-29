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


// File         : fft.cxx
// Author       : Pavel A. Koshevoy
// Created      : 2005/06/03 10:16
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : Wrapper class and helper functions for working with
//                FFTW3 and ITK images.

// system includes:
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <string.h>
#include <math.h>

// ITK includes:
#include <itkImage.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkForwardFFTImageFilter.h>
#include <itkInverseFFTImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkComposeImageFilter.h>
#include <itkComplexToRealImageFilter.h>
#include <itkComplexToImaginaryImageFilter.h>
#include <itkComplexToComplexFFTImageFilter.h>

// local includes:
#include "IRFFT.h"
#include "IRMutexInterface.h"
#include "itkIRUtils.h"


namespace itk_fft
{
	// static the_mutex_interface_t * mutex = nullptr; 
	// //----------------------------------------------------------------
	// // fftw_mutex
	// // 
	// static the_mutex_interface_t *
	// 	fftw_mutex()
	// {
	// 	if(mutex == nullptr)
	// 		mutex = the_mutex_interface_t::create();
	// 	return mutex;
	// }

	//----------------------------------------------------------------
	// NUM_FFTW_THREADS
	// 
	static std::size_t NUM_FFTW_THREADS = 1;

	//----------------------------------------------------------------
	// set_num_fftw_threads
	// 
	// set's number of threads used by fftw, returns previous value:
	// std::size_t set_num_fftw_threads(std::size_t num_threads)
	// {
	// 	the_lock_t<the_mutex_interface_t> lock(fftw_mutex());
	// 	std::size_t prev = NUM_FFTW_THREADS;
	// 	if (num_threads > 0)
	// 	{

	// 		NUM_FFTW_THREADS = num_threads;
	// 	}

	// 	return prev;
	// }


	// //----------------------------------------------------------------
	// // fft_cache_t
	// // 
	// class fft_cache_t
	// {
	// public:
	// 	fft_cache_t():
	// 	  fwd_(nullptr),
	// 		  inv_(nullptr),
	// 		  w_(0),
	// 		  h_(0),
	// 		  h_complex_(0),
	// 		  h_padded_(0),
	// 		  buffer_(nullptr)
	// 	  {}

	// 	  ~fft_cache_t()
	// 	  {
	// 		  clear();
	// 	  }

	// 	  void clear()
	// 	  {
	// 		  if (buffer_)
	// 		  {
	// 			  // fftw is not thread safe:
	// 			//   the_lock_t<the_mutex_interface_t> lock(fftw_mutex());

	// 			  fftwf_destroy_plan(fwd_);
	// 			  fwd_ = nullptr;

	// 			  fftwf_destroy_plan(inv_);
	// 			  inv_ = nullptr;

	// 			  fftwf_free(buffer_);
	// 			  buffer_ = nullptr;
	// 		  }

	// 		  w_ = 0;
	// 		  h_ = 0;
	// 		  h_complex_ = 0;
	// 		  h_padded_ = 0;
	// 	  }

	// 	  void update(const unsigned int & w, const unsigned int & h)
	// 	  {
	// 		  if (w_ == w && h_ == h)
	// 		  {
	// 			  return;
	// 		  }

	// 		  // fftw is not thread safe:
	// 		//   the_lock_t<the_mutex_interface_t> lock(fftw_mutex());

	// 		//   if (buffer_)
	// 		//   {
	// 		// 	  fftwf_destroy_plan(fwd_);
	// 		// 	  fftwf_destroy_plan(inv_);
	// 		// 	  fftwf_free(buffer_);
	// 		//   }

	// 		  w_ = w;
	// 		  h_ = h;

	// 		  h_complex_ = h_ / 2 + 1;
	// 		  h_padded_ = h_complex_ * 2;

	// 		  buffer_ = (float *)(fftwf_malloc(w_ * h_padded_ * sizeof(float)));
	// 		  fftwf_plan_with_nthreads(NUM_FFTW_THREADS);
	// 		  fwd_ = fftwf_plan_dft_r2c_2d(w_,
	// 			  h_,
	// 			  buffer_,
	// 			  (fftwf_complex *)buffer_,
	// 			  FFTW_DESTROY_INPUT | FFTW_ESTIMATE);

	// 		  fft_data_t dummy_in(4, 4);
	// 		  fft_data_t dummy_out(4, 4);
	// 		  fftwf_plan_with_nthreads(NUM_FFTW_THREADS);
	// 		  inv_ = fftwf_plan_dft_2d(w_,
	// 			  h_,
	// 			  (fftwf_complex *)(dummy_in.data()),
	// 			  (fftwf_complex *)(dummy_out.data()),
	// 			  FFTW_BACKWARD,
	// 			  FFTW_ESTIMATE);
	// 	  }

	// 	  fftwf_plan fwd_; // analysis, forward fft
	// 	  fftwf_plan inv_; // synthesis, inverse fft
	// 	  unsigned int w_;
	// 	  unsigned int h_;
	// 	  unsigned int h_complex_;
	// 	  unsigned int h_padded_;
	// 	  float * buffer_;
	// };

	//----------------------------------------------------------------
	// tss
	// 
	// static boost::thread_specific_ptr<fft_cache_t> tss;
	using ForwardFFTFilterType = itk::ForwardFFTImageFilter<itk_image_t, itk_complex_image_t>;
	static thread_local ForwardFFTFilterType::Pointer forwardFFTFilter = ForwardFFTFilterType::New();
	using InverseFFTFilterType = itk::ComplexToComplexFFTImageFilter<itk_complex_image_t, itk_complex_image_t>;
	static thread_local InverseFFTFilterType::Pointer inverseFFTFilter = InverseFFTFilterType::New();

	//----------------------------------------------------------------
	// fft_data_t::fft_data_t
	// 
	fft_data_t::fft_data_t(const itk_image_t::Pointer & real):
		image_(),
		nx_(0),
		ny_(0)
	{
		setup(real);
	}

	//----------------------------------------------------------------
	// fft_data_t::fft_data_t
	// 
	fft_data_t::fft_data_t(const unsigned int w, const unsigned int h):
		image_(),
		nx_(0),
		ny_(0)
	{
		resize(w, h);
	}

	//----------------------------------------------------------------
	// fft_data_t::fft_data_t
	// 
	fft_data_t::fft_data_t(const itk_image_t::Pointer & real,
		const itk_image_t::Pointer & imag):
		image_(),
		nx_(0),
		ny_(0)
	{
		setup(real, imag);
	}

	//----------------------------------------------------------------
	// fft_data_t::fft_data_t
	// 
	fft_data_t::fft_data_t(const fft_data_t & data):
		image_(),
		nx_(0),
		ny_(0)
	{
		*this = data;
	}

	//----------------------------------------------------------------
	// fft_data_t::operator =
	// 
	fft_data_t &
		fft_data_t::operator = (const fft_data_t & data)
	{
		assert(this != &data);

		if (data.image_ != nullptr)
		{
			using DuplicatorType = itk::ImageDuplicator<itk_complex_image_t>;
			DuplicatorType::Pointer duplicator = DuplicatorType::New();
			duplicator->SetInputImage(data.image_);
			duplicator->Update();
			image_ = duplicator->GetOutput();
			resize(data.nx_, data.ny_);
		}
		else
		{
			cleanup();
		}

		return *this;
	}

	//----------------------------------------------------------------
	// fft_data_t::cleanup
	// 
	void
		fft_data_t::cleanup()
	{
		// if (data_ != nullptr) fftwf_free((fft_complex_t *)(data_));
		// data_ = nullptr;
		nx_ = 0;
		ny_ = 0;
	}

	//----------------------------------------------------------------
	// fft_data_t::resize
	// 
	void
		fft_data_t::resize(const unsigned int w, const unsigned int h)
	{
		const unsigned int new_sz = w  * h;
		const unsigned int old_sz = nx_ * ny_;
		if (old_sz == new_sz) return;

		// if (data_ != nullptr) fftwf_free((fft_complex_t *)(data_));
		// data_ = (fft_complex_t *)(fftwf_malloc(new_sz * sizeof(fft_complex_t)));
		// assert(data_ != nullptr);
		image_->SetRegions({ w, h });
		image_->Allocate();

		nx_ = w;
		ny_ = h;
	}

	//----------------------------------------------------------------
	// fft_data_t::fill
	// 
	void
		fft_data_t::fill(const float real, const float imag)
	{
		const fft_complex_t value(real, imag);
		this->image_->FillBuffer(value);
	}

	//----------------------------------------------------------------
	// fft_data_t::setup
	// 
	void
		fft_data_t::setup(const itk_image_t::Pointer & real,
		const itk_image_t::Pointer & imag)
	{
		using ComposeFilterType = itk::ComposeImageFilter<itk_image_t, itk_complex_image_t>;
		ComposeFilterType::Pointer composeFilter = ComposeFilterType::New();
		composeFilter->SetInput1(real);
		composeFilter->SetInput2(imag);
		composeFilter->Update();
		image_ = composeFilter->GetOutput();
	}

	//----------------------------------------------------------------
	// fft_data_t::component
	// 
	itk_image_t::Pointer
		fft_data_t::component(const bool imag) const
	{
		if (imag)
		{
			using ComplexToImaginaryFilterType = itk::ComplexToImaginaryImageFilter<itk_complex_image_t, itk_image_t>;
			ComplexToImaginaryFilterType::Pointer complexToImaginaryFilter = ComplexToImaginaryFilterType::New();
			complexToImaginaryFilter->SetInput(image_);
			complexToImaginaryFilter->Update();
			return complexToImaginaryFilter->GetOutput();
		}
		using ComplexToRealFilterType = itk::ComplexToRealImageFilter<itk_complex_image_t, itk_image_t>;
		ComplexToRealFilterType::Pointer complexToRealFilter = ComplexToRealFilterType::New();
		complexToRealFilter->SetInput(image_);
		complexToRealFilter->Update();
		return complexToRealFilter->GetOutput();
	}

	static std::vector<double> syLookup;
	static std::vector<double> sySquaredLookup;

	//----------------------------------------------------------------
	// fft_data_t::apply_lp_filter
	// 
	void
		fft_data_t::apply_lp_filter(const double r, const double s)
	{
		if (r > ::sqrt(2.0)) return;

		const unsigned int hx = nx_ / 2;
		const unsigned int hy = ny_ / 2;

		const double r0 = (r - s);
		const double r1 = (r + s);
		const double dr = r1 - r0;

		const double r0sqr = r0 * r0;
		const double r1sqr = r1 * r1;

		//The use of this static cache needs to be wrapped in a lock, copied, or some other protection if multithreaded calls with different image sizes are used.
		if((unsigned int)syLookup.size() != ny_)
		{
			syLookup.clear(); 
			sySquaredLookup.clear(); 

			syLookup.reserve(ny_);
			sySquaredLookup.reserve(ny_);
			for (unsigned int y = 0; y < ny_; y++)
			{
				double sy = 2.0 * (double((y + hy) % ny_) - double(hy)) / double(ny_);
				double y2 = sy * sy;
				syLookup.push_back(sy);
				sySquaredLookup.push_back(sy); 
			}
		}

		fft_complex_t * data_ = this->image_->GetBufferPointer();
		unsigned int i = 0; 
		for (unsigned int x = 0; x < nx_; x++)
		{
			double sx = 2.0 * (double((x + hx) % nx_) - double(hx)) / double(nx_);
			double x2 = sx * sx;

			for (unsigned int y = 0; y < ny_; y++)
			{
				double sy = syLookup[y];
				double y2 = sySquaredLookup[y];

				double d2 = x2 + y2;
				double v;

				if(d2 < r0sqr)
					v = 1.0;
				else if(d2 > r1sqr)
					v = 0.0; 
				else
				{
					double d = ::sqrt(x2 + y2);
					v = (1.0 + cos(M_PI * (d - r0) / dr)) / 2.0;
				}

				data_[i++] *= v;
			}
		}
	}


	//----------------------------------------------------------------
	// fft
	// 
	bool
		fft(itk_image_t::ConstPointer & in, fft_data_t & out)
	{
		forwardFFTFilter->SetInput(in);
		forwardFFTFilter->GraftOutput(out.data());
		forwardFFTFilter->UpdateLargestPossibleRegion();

		return true;
	}


	//----------------------------------------------------------------
	// ifft
	// 
	bool
		ifft(const fft_data_t & in, fft_data_t & out)
	{
		out.resize(in.nx(), in.ny());
		inverseFFTFilter->SetInput(in.data());
		inverseFFTFilter->GraftOutput(out.data());
		inverseFFTFilter->UpdateLargestPossibleRegion();

		return true;
	}

	// //----------------------------------------------------------------
	// // elem_by_elem
	// // 
	// void
	// 	elem_by_elem(fn_fft_c_t f,
	// 	const fft_data_t & in,
	// 	fft_data_t & out)
	// {
	// 	const unsigned nx = in.nx();
	// 	const unsigned ny = in.ny();

	// 	out.resize(nx, ny);
	// 	fft_complex_t * dout = out.data();
	// 	const fft_complex_t * din = in.data();

	// 	const unsigned int size = nx * ny;
	// 	for (unsigned int i = 0; i < size; i++)
	// 	{
	// 		dout[i] = f(din[i]);
	// 	}
	// }


	// //----------------------------------------------------------------
	// // elem_by_elem
	// //
	// void
	// 	elem_by_elem(fft_complex_t(*f)(const double & a,
	// 	const fft_complex_t & b),
	// 	const double & a,
	// 	const fft_data_t & b,
	// 	fft_data_t & c)
	// {
	// 	assert(&b != &c);

	// 	const unsigned nx = b.nx();
	// 	const unsigned ny = b.ny();
	// 	c.resize(nx, ny);

	// 	fft_complex_t * dc = c.data();
	// 	const fft_complex_t * db = b.data();

	// 	const unsigned int size = nx * ny;
	// 	for (unsigned int i = 0; i < size; i++)
	// 	{
	// 		dc[i] = f(a, db[i]);
	// 	}
	// }


	// //----------------------------------------------------------------
	// // elem_by_elem
	// // 
	// void
	// 	elem_by_elem(fft_complex_t(*f)(const fft_complex_t & a,
	// 	const double & b),
	// 	const fft_data_t & a,
	// 	const double & b,
	// 	fft_data_t & c)
	// {
	// 	assert(&a != &c);

	// 	const unsigned nx = a.nx();
	// 	const unsigned ny = a.ny();
	// 	c.resize(nx, ny);

	// 	fft_complex_t * dc = c.data();
	// 	const fft_complex_t * da = a.data();

	// 	const unsigned int size = nx * ny;
	// 	for (unsigned int i = 0; i < size; i++)
	// 	{
	// 		dc[i] = f(da[i], b);
	// 	}
	// }


	// //----------------------------------------------------------------
	// // elem_by_elem
	// //
	// void
	// 	elem_by_elem(fft_complex_t(*f)(const fft_complex_t & a,
	// 	const fft_complex_t & b),
	// 	const fft_complex_t & a,
	// 	const fft_data_t & b,
	// 	fft_data_t & c)
	// {
	// 	assert(&b != &c);

	// 	const unsigned nx = b.nx();
	// 	const unsigned ny = b.ny();
	// 	c.resize(nx, ny);

	// 	fft_complex_t * dc = c.data();
	// 	const fft_complex_t * db = b.data();

	// 	const unsigned int size = nx * ny;
	// 	for (unsigned int i = 0; i < size; i++)
	// 	{
	// 		dc[i] = f(a, db[i]);
	// 	}
	// }


	// //----------------------------------------------------------------
	// // elem_by_elem
	// // 
	// void
	// 	elem_by_elem(fft_complex_t(*f)(const fft_complex_t & a,
	// 	const fft_complex_t & b),
	// 	const fft_data_t & a,
	// 	const fft_complex_t & b,
	// 	fft_data_t & c)
	// {
	// 	assert(&a != &c);

	// 	const unsigned nx = a.nx();
	// 	const unsigned ny = a.ny();
	// 	c.resize(nx, ny);

	// 	fft_complex_t * dc = c.data();
	// 	const fft_complex_t * da = a.data();

	// 	const unsigned int size = nx * ny;
	// 	for (unsigned int i = 0; i < size; i++)
	// 	{
	// 		dc[i] = f(da[i], b);
	// 	}
	// }


	// //----------------------------------------------------------------
	// // elem_by_elem
	// //
	// void
	// 	elem_by_elem(fn_fft_cc_t f,
	// 	const fft_data_t & a,
	// 	const fft_data_t & b,
	// 	fft_data_t & c)
	// {
	// 	assert(&a != &c);
	// 	assert(&b != &c);

	// 	const unsigned nx = a.nx();
	// 	const unsigned ny = a.ny();
	// 	assert(nx == b.nx() && ny == b.ny());

	// 	c.resize(nx, ny);
	// 	fft_complex_t * dc = c.data();

	// 	const fft_complex_t * da = a.data();
	// 	const fft_complex_t * db = b.data();

	// 	const unsigned int size = nx * ny;
	// 	for (unsigned int i = 0; i < size; i++)
	// 	{
	// 		dc[i] = f(da[i], db[i]);
	// 	}
	// }

} // namespace itk_fft
