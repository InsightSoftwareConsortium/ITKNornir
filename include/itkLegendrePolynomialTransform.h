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


// File         : itkLegendrePolynomialTransform.h
// Author       : Pavel A. Koshevoy
// Created      : 2006/02/22 15:55
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : A bivariate centered/normalized Legendre polynomial
//                transform and helper functions.

#ifndef __itkLegendrePolynomialTransform_h
#define __itkLegendrePolynomialTransform_h

// system includes:
#include <iostream>
#include <assert.h>

// ITK includes:
#include <itkTransform.h>
#include <itkMacro.h>
#include <itkImage.h>

// local includes:
#include "itkInverseTransform.h"

// VXL includes:
#include <vnl/algo/vnl_svd.h>


//----------------------------------------------------------------
// itk::LegendrePolynomialTransform
//
// Let
//   A = (u - uc) / Xmax
//   B = (v - vc) / Ymax
//
// where uc, vc correspond to the center of the image expressed in
// the coordinate system of the mosaic.
//
// The transform is defined as
//   x(u, v) = Xmax * Sa
//   y(u, v) = Ymax * Sb
//
// where
//   Sa = sum(i in [0, N], sum(j in [0, i], a_jk * Pj(A) * Qk(B)));
//   Sb = sum(i in [0, N], sum(j in [0, i], b_jk * Pj(A) * Qk(B)));
//
// where k = i - j and (Pj, Qk) are Legendre polynomials
// of degree (j, k) respectively.
//
namespace itk
{
template <class TScalar = itk::SpacePrecisionType, unsigned int N = 2>
class LegendrePolynomialTransform : public Transform<TScalar, 2, 2>
{
public:
  // standard typedefs:
  typedef LegendrePolynomialTransform Self;
  typedef SmartPointer<Self>          Pointer;
  typedef SmartPointer<const Self>    ConstPointer;

  typedef Transform<TScalar, 2, 2> Superclass;

  // Base inverse transform type:
  typedef Superclass                         InverseTransformType;
  typedef SmartPointer<InverseTransformType> InverseTransformPointer;

  // static constant for the degree of the polynomial:
  itkStaticConstMacro(Degree, unsigned int, N);

  // static constant for the number of a_jk (or b_jk) coefficients:
  itkStaticConstMacro(CoefficientsPerDimension, unsigned int, ((N + 1) * (N + 2)) / 2);

  // static constant for the length of the parameter vector:
  itkStaticConstMacro(ParameterVectorLength, unsigned int, (N + 1) * (N + 2));

  // RTTI:
  itkOverrideGetNameOfClassMacro(LegendrePolynomialTransform);

  // macro for instantiation through the object factory:
  itkNewMacro(Self);

  /** Standard scalar type for this class. */
  typedef typename Superclass::ScalarType ScalarType;

  /** Dimension of the domain space. */
  itkStaticConstMacro(InputSpaceDimension, unsigned int, 2);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 2);

  // shortcuts:
  using FixedParametersType = typename Superclass::FixedParametersType;
  typedef typename Superclass::ParametersType ParametersType;
  typedef typename Superclass::JacobianType   JacobianType;

  typedef typename Superclass::InputPointType  InputPointType;
  typedef typename Superclass::OutputPointType OutputPointType;

  // virtual:
  OutputPointType
  TransformPoint(const InputPointType & x) const override;

  // Inverse transformations:
  // If y = Transform(x), then x = BackTransform(y);
  // if no mapping from y to x exists, then an exception is thrown.
  InputPointType
  BackTransformPoint(const OutputPointType & y) const;

  using InputDiffusionTensor3DType = typename Superclass::InputDiffusionTensor3DType;
  using OutputDiffusionTensor3DType = typename Superclass::OutputDiffusionTensor3DType;
  OutputDiffusionTensor3DType
  TransformDiffusionTensor3D(const InputDiffusionTensor3DType & inputTensor, const InputPointType & point) const override
  {
    itkExceptionMacro("TransformDiffusionTensor3D( const InputDiffusionTensor3DType & ) is "
                      "unimplemented for "
                      << this->GetNameOfClass());
  }

  using InputVectorPixelType = typename Superclass::InputVectorPixelType;
  using OutputVectorPixelType = typename Superclass::OutputVectorPixelType;
  OutputVectorPixelType
  TransformDiffusionTensor3D(const InputVectorPixelType & inputTensor, const InputPointType & point) const override
  {
    itkExceptionMacro("TransformDiffusionTensor3D( const InputVectorPixelType & ) is "
                      "unimplemented for "
                      << this->GetNameOfClass());
  }

  // virtual:
  void
  SetFixedParameters(const FixedParametersType & params) override
  {
    this->m_FixedParameters = params;
  }

  // virtual:
  const FixedParametersType &
  GetFixedParameters() const override
  {
    return this->m_FixedParameters;
  }

  // virtual:
  void
  SetParameters(const ParametersType & params) override
  {
    this->m_Parameters = params;
  }

  // virtual:
  const ParametersType &
  GetParameters() const override
  {
    return this->m_Parameters;
  }

  // virtual:
  IdentifierType
  GetNumberOfParameters() const override
  {
    return ParameterVectorLength;
  }

  void
  ComputeJacobianWithRespectToParameters(const InputPointType & point, JacobianType & jacobian) const override;

  // virtual: return an inverse of this transform.
  InverseTransformPointer
  GetInverse() const
  {
    typedef InverseTransform<Self>     InvTransformType;
    typename InvTransformType::Pointer inv = InvTransformType::New();
    inv->SetForwardTransform(this);
    return inv.GetPointer();
  }

  // setup the fixed transform parameters:
  void
  setup( // image bounding box expressed in the physical space:
    const double x_min,
    const double x_max,
    const double y_min,
    const double y_max,

    // normalization parameters:
    const double Xmax = 0.0,
    const double Ymax = 0.0)
  {
    double & uc_ = this->m_FixedParameters[0];
    double & vc_ = this->m_FixedParameters[1];
    double & xmax_ = this->m_FixedParameters[2];
    double & ymax_ = this->m_FixedParameters[3];
    double & a00_ = this->m_Parameters[index_a(0, 0)];
    double & b00_ = this->m_Parameters[index_b(0, 0)];

    // center of the image:
    double xc = (x_min + x_max) / 2.0;
    double yc = (y_min + y_max) / 2.0;
    uc_ = xc;
    vc_ = yc;

    // setup the normalization parameters:
    if (Xmax != 0.0 && Ymax != 0.0)
    {
      xmax_ = Xmax;
      ymax_ = Ymax;
    }
    else
    {
      const double w = x_max - x_min;
      const double h = y_max - y_min;

      // -1 : 1
      xmax_ = w / 2.0;
      ymax_ = h / 2.0;
    }

    // setup a00, b00 (local translation parameters):
    a00_ = xc / xmax_;
    b00_ = yc / ymax_;
  }

  // setup the translation parameters:
  void
  setup_translation( // translation is expressed in the physical space:
    const double tx_Xmax = 0.0,
    const double ty_Ymax = 0.0)
  {
    // incorporate translation into the (uc, vc) fixed parameters:
    double & uc_ = this->m_FixedParameters[0];
    double & vc_ = this->m_FixedParameters[1];

    // FIXME: the signs might be wrong here (20051101):
    uc_ -= tx_Xmax;
    vc_ -= ty_Ymax;

    // FIXME: untested:
    /*
    double & a00_ = this->m_Parameters[index_a(0, 0)];
    double & b00_ = this->m_Parameters[index_b(0, 0)];
    a00_ -= tx_Xmax / GetXmax();
    b00_ -= ty_Ymax / GetYmax();
    */
  }

  // helper required for numeric inverse transform calculation;
  // evaluate F = T(x), J = dT/dx (another Jacobian):
  void
  eval(const std::vector<ScalarType> & x, std::vector<ScalarType> & F, std::vector<std::vector<ScalarType>> & J) const;

  // setup a linear system to solve for the parameters of this
  // transform such that it maps points uv to xy:
  void
  setup_linear_system(const unsigned int                   start_with_degree,
                      const unsigned int                   degrees_covered,
                      const std::vector<InputPointType> &  uv,
                      const std::vector<OutputPointType> & xy,
                      vnl_matrix<double> &                 M,
                      vnl_vector<double> &                 bx,
                      vnl_vector<double> &                 by) const;

  // find the polynomial coefficients such that this transform
  // would map uv to xy:
  void
  solve_for_parameters(const unsigned int                   start_with_degree,
                       const unsigned int                   degrees_covered,
                       const std::vector<InputPointType> &  uv, // mosaic
                       const std::vector<OutputPointType> & xy, // tile
                       ParametersType &                     params) const;

  inline void
  solve_for_parameters(const unsigned int                   start_with_degree,
                       const unsigned int                   degrees_covered,
                       const std::vector<InputPointType> &  uv,
                       const std::vector<OutputPointType> & xy)
  {
    ParametersType params = GetParameters();
    solve_for_parameters(start_with_degree, degrees_covered, uv, xy, params);
    SetParameters(params);
  }

  // calculate the number of coefficients of a given
  // degree range (per dimension):
  inline static unsigned int
  count_coefficients(const unsigned int start_with_degree, const unsigned int degrees_covered)
  {
    return index_a(0, start_with_degree + degrees_covered) - index_a(0, start_with_degree);
  }

  // accessors to the warp origin expressed in the mosaic coordinate system:
  inline const double &
  GetUc() const
  {
    return this->m_FixedParameters[0];
  }

  inline const double &
  GetVc() const
  {
    return this->m_FixedParameters[1];
  }

  // accessors to the normalization parameters Xmax, Ymax:
  inline const double &
  GetXmax() const
  {
    return this->m_FixedParameters[2];
  }

  inline const double &
  GetYmax() const
  {
    return this->m_FixedParameters[3];
  }

  // generate a mask of shared parameters:
  static void
  setup_shared_params_mask(bool shared, std::vector<bool> & mask)
  {
    mask.assign(ParameterVectorLength, shared);
    mask[index_a(0, 0)] = false;
    mask[index_b(0, 0)] = false;
  }

  // Convert the j, k indeces associated with a(j, k) coefficient
  // into an index that can be used with the parameters array:
  inline static unsigned int
  index_a(const unsigned int & j, const unsigned int & k)
  {
    return j + ((j + k) * (j + k + 1)) / 2;
  }

  inline static unsigned int
  index_b(const unsigned int & j, const unsigned int & k)
  {
    return CoefficientsPerDimension + index_a(j, k);
  }

  // virtual: Generate a platform independent name:
  std::string
  GetTransformTypeAsString() const override
  {
    std::string        base = Superclass::GetTransformTypeAsString();
    std::ostringstream name;
    name << base << '_' << N;
    return name.str();
  }

protected:
  LegendrePolynomialTransform();

  // virtual:
  void
  PrintSelf(std::ostream & s, Indent indent) const override;

private:
  // disable default copy constructor and assignment operator:
  LegendrePolynomialTransform(const Self & other);
  const Self &
  operator=(const Self & t);

  // polynomial coefficient accessors:
  inline const double &
  a(const unsigned int & j, const unsigned int & k) const
  {
    return this->m_Parameters[index_a(j, k)];
  }

  inline const double &
  b(const unsigned int & j, const unsigned int & k) const
  {
    return this->m_Parameters[index_b(j, k)];
  }

}; // class LegendrePolynomialTransform

} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkLegendrePolynomialTransform.hxx"
#endif


//----------------------------------------------------------------
// setup_transform
//
template <class transform_t>
typename transform_t::Pointer
setup_transform(const itk::Point<double, 2> & bbox_min, const itk::Point<double, 2> & bbox_max)
{
  double w = bbox_max[0] - bbox_min[0];
  double h = bbox_max[1] - bbox_min[1];
  double Umax = w / 2.0;
  double Vmax = h / 2.0;

  typename transform_t::Pointer t = transform_t::New();
  t->setup(bbox_min[0], bbox_max[0], bbox_min[1], bbox_max[1], Umax, Vmax);

  return t;
}


//----------------------------------------------------------------
// setup_transform
//
template <class transform_t, class image_t>
typename transform_t::Pointer
setup_transform(const image_t * image)
{
  typedef typename image_t::IndexType index_t;

  index_t i00;
  i00[0] = 0;
  i00[1] = 0;

  typename image_t::PointType origin;
  image->TransformIndexToPhysicalPoint(i00, origin);

  index_t i11;
  i11[0] = 1;
  i11[1] = 1;

  typename image_t::PointType spacing;
  image->TransformIndexToPhysicalPoint(i11, spacing);
  spacing[0] -= origin[0];
  spacing[1] -= origin[1];

  typename image_t::SizeType sz = image->GetLargestPossibleRegion().GetSize();

  typename image_t::PointType bbox_min = origin;
  typename image_t::PointType bbox_max;
  bbox_max[0] = bbox_min[0] + spacing[0] * double(sz[0]);
  bbox_max[1] = bbox_min[1] + spacing[1] * double(sz[1]);

  return setup_transform<transform_t>(bbox_min, bbox_max);
}


#endif // __itkLegendrePolynomialTransform_h
