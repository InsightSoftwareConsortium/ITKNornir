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


// File         : itkInverseTransform.h
// Author       : Pavel A. Koshevoy
// Created      : 2005/06/03 10:16
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : A inverse transform class.

#ifndef __itkInverseTransform_h
#define __itkInverseTransform_h

// system includes:
#include <assert.h>

// ITK includes:
#include <itkTransform.h>
#include <itkMacro.h>


//----------------------------------------------------------------
// itk::InverseTransform
// 
namespace itk
{
  //----------------------------------------------------------------
  // InverseTransform
  //
  template <class ForwardTransform>
  class InverseTransform :
    public Transform< typename ForwardTransform::ScalarType,
		      ForwardTransform::OutputSpaceDimension,
		      ForwardTransform::InputSpaceDimension >
  {
  public:
    /** Standard class typedefs. */
    typedef InverseTransform Self;
    
    typedef Transform< typename ForwardTransform::ScalarType,
		       ForwardTransform::OutputSpaceDimension,
		       ForwardTransform::InputSpaceDimension > Superclass;
    
    typedef SmartPointer< Self >	Pointer;
    typedef SmartPointer< const Self >	ConstPointer;
    
    /** Base inverse transform type. */
    typedef typename Superclass::InverseTransformType InverseTransformType;
    typedef SmartPointer< InverseTransformType > InverseTransformPointer;
    
    /** New method for creating an object using a factory. */
    itkNewMacro(Self);
    
    /** Run-time type information (and related methods). */
    itkTypeMacro( InverseTransform, Transform );
    
    /** Standard scalar type for this class. */
    typedef typename Superclass::ScalarType ScalarType;
    
    /** Type of the input parameters. */
    typedef typename Superclass::ParametersType ParametersType;
    
    /** Type of the Jacobian matrix. */
    typedef typename Superclass::JacobianType JacobianType;
    
    /** Standard coordinate point type for this class. */
    typedef typename Superclass::InputPointType InputPointType;
    typedef typename Superclass::OutputPointType OutputPointType;
    
    /** Dimension of the domain space. */
    itkStaticConstMacro(InputSpaceDimension,
			unsigned int,
			ForwardTransform::OutputSpaceDimension);
    itkStaticConstMacro(OutputSpaceDimension,
			unsigned int,
			ForwardTransform::InputSpaceDimension);
    
    /** Set the forward transform pointer. */
    void SetForwardTransform(const ForwardTransform * forward)
    { forward_ = forward; }
    
    /**  Method to transform a point. */
    virtual OutputPointType TransformPoint(const InputPointType & y) const
    {
      assert(forward_ != nullptr);
      return forward_->BackTransformPoint(y);
    }
    
    virtual const JacobianType & GetJacobian(const InputPointType &) const
    {
      itkExceptionMacro(<< "GetJacobian is not implemented "
			"for InverseTransform");
      return this->m_Jacobian;
    };
    
    virtual unsigned int GetNumberOfParameters() const 
    { return 0; }
    
    virtual InverseTransformPointer GetInverse() const
    { return const_cast<ForwardTransform *>(forward_); }
    
  protected:
    InverseTransform(): Superclass(0, 0) {}
    
  private:
    // disable default copy constructor and assignment operator:
    InverseTransform(const Self & other);
    const Self & operator = (const Self & t);
    
    // the transform whose inverse we are trying to evaluate:
    const ForwardTransform * forward_;
  };
  
} // namespace itk

#endif // __itkInverseTransform_h
