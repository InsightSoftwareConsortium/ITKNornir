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


// File         : itkNumericTransform.h
// Author       : Pavel A. Koshevoy
// Created      : 2005/06/03 10:16
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : A numerical inverse transform class, based on the
//                Newton-Raphson method for nonlinear systems of equations.

#ifndef __itkNumericInverse_h
#define __itkNumericInverse_h

// system includes:
#include <iostream>
#include <assert.h>
#include <vector>

// ITK includes:
#include <itkTransform.h>
#include <itkMacro.h>

// VXL includes:
#include <vnl/algo/vnl_svd.h>

namespace help
{
  //----------------------------------------------------------------
  // nonlinear_system_evaluator_t
  // 
  template <class ScalarType>
  class nonlinear_system_evaluator_t
  {
  public:
    virtual ~nonlinear_system_evaluator_t() {}
    
    // evaluate the nonlinear system of equations
    // and its Jacobian (dF/dx) at a given point:
    virtual void
    operator() (const std::vector<ScalarType> & x,
		std::vector<ScalarType> & F,
		std::vector<std::vector<ScalarType> > & J) const = 0;
  };
  
  //----------------------------------------------------------------
  // NewtonRaphson
  // 
  template <class ScalarType>
  bool
  NewtonRaphson(const nonlinear_system_evaluator_t<ScalarType> & usrfun,
		std::vector<ScalarType> & x, // estimated root point
		const unsigned int & ntrial, // number of iterations
		const ScalarType & tolx, // convergence tolerance in x
		const ScalarType & tolf) // convergence tolerance in F
  {
    // shortcut:
    const unsigned int n = x.size();
    
    std::vector<ScalarType> F(n);
    std::vector<std::vector<ScalarType> > J(n);
    for (unsigned int i = 0; i < n; i++) J[i].resize(n);
    
    vnl_matrix<ScalarType> A(n, n);
    vnl_vector<ScalarType> b(n);
    
    for (unsigned int k = 0; k < ntrial; k++)
    {
      // evaluate the function at the current position:
      usrfun(x, F, J);
      
      // check for function convergence:
      ScalarType errf = ScalarType(0);
      for (unsigned int i = 0; i < n; i++)
      {
	errf += fabs(F[i]);
      }
      if (errf <= tolf) break;
      
      // setup the left hand side of the linear system:
      for (unsigned int i = 0; i < n; i++)
      {
	for (unsigned int j = 0; j < n; j++)
	{
	  A[i][j] = J[i][j];
	}
      }
      
      // setup the right hand side of the linear system:
      for (unsigned int i = 0; i < n; i++)
      {
	b[i] = -F[i];
      }
      
      vnl_svd<double> svd(A);
      vnl_vector<double> dx = svd.solve(b);
      
      // check for root convergence:
      ScalarType errx = ScalarType(0);
      for (unsigned int i = 0; i < n; i++)
      {
	errx += fabs(dx[i]);
	x[i] += dx[i];
      }
      if (errx <= tolx) break;
    }
    
    return true;
  }
  
} // namespace help

namespace itk
{
  //----------------------------------------------------------------
  // NumericInverse
  //
  template <class TTransform>
  class NumericInverse :
    public help::nonlinear_system_evaluator_t<typename TTransform::ScalarType>
  {
  public:
    typedef TTransform TransformType;
    typedef typename TTransform::ScalarType ScalarType;
    
    NumericInverse(const TransformType & transform):
      transform_(transform)
    {}
    
    // virtual:
    void operator() (const std::vector<ScalarType> & x,
		     std::vector<ScalarType> & F,
		     std::vector<std::vector<ScalarType> > & J) const
    {
      transform_.eval(x, F, J);
      
      const unsigned int n = x.size();
      for (unsigned int i = 0; i < n; i++)
      {
	F[i] -= y_[i];
      }
    }
    
    // If y = Transform(x), then x = BackTransform(y).
    // given y, find x:
    bool transform(const std::vector<ScalarType> & y,
		   std::vector<ScalarType> & x,
		   bool x_is_initialized = false) const
    {
      y_ = y;
      if (!x_is_initialized) x  = y;
      return NewtonRaphson(*this, x, 50, 1e-12, 1e-12);
    }
    
  private:
    // the transform whose inverse we are trying to evaluate:
    const TransformType & transform_;
    
    // the point for which we are tryying to find the inverse mapping:
    mutable std::vector<ScalarType> y_;
  };
  
} // namespace itk

#endif // __itkNumericInverse_h
