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


// File         : itkRBFTransform.cxx
// Author       : Pavel A. Koshevoy
// Created      : 2007/01/23 10:15
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : A Thin Plate Spline transform.

// local includes:
#include "itkRBFTransform.h"

// ITK includes:
#include <vnl/algo/vnl_svd.h>

// system includes:
#include <iostream>

// namespace access:
using std::cerr;
using std::endl;


namespace itk
{
#if 0
}
#endif

//----------------------------------------------------------------
// RBFTransform::RBFTransform
//
RBFTransform::RBFTransform()
{
  this->m_FixedParameters.SetSize(4);
  this->m_Parameters.SetSize(6);

  double & Xmax = this->m_FixedParameters[0];
  double & Ymax = this->m_FixedParameters[1];
  Xmax = 1;
  Ymax = 1;

  double & uc = this->m_FixedParameters[2];
  double & vc = this->m_FixedParameters[3];
  uc = 0;
  vc = 0;

  double * ab_vec = &(this->m_Parameters[index_a(0)]);
  ab_vec[0] = 0;
  ab_vec[1] = 0;
  ab_vec[2] = 1;
  ab_vec[3] = 0;
  ab_vec[4] = 0;
  ab_vec[5] = 1;
}

//----------------------------------------------------------------
// RBFTransform::TransformPoint
//
RBFTransform::OutputPointType
RBFTransform::TransformPoint(const InputPointType & uv) const
{
  OutputPointType xy;

  const double & Xmax = GetXmax();
  const double & Ymax = GetYmax();

  const double & uc = GetUc();
  const double & vc = GetVc();

  const double * uv_vec = this->uv(0);
  const double * fg_vec = &(this->f(0));
  const double * ab_vec = &(this->a(0));

  const double & u = uv[0];
  const double & v = uv[1];

  // calculate the summation terms:
  double F = 0;
  double G = 0;
  {
    unsigned int num_pts = num_points();
    for (unsigned int i = 0; i < num_pts; i++)
    {
      unsigned int offset = i * 2;

      // calculate the kernel term:
      double Q = kernel(&(uv[0]), uv_vec + offset, Xmax, Ymax);
      F += Q * fg_vec[offset];
      G += Q * fg_vec[offset + 1];
    }
  }

  double A = (u - uc) / Xmax;
  double B = (v - vc) / Ymax;
  xy[0] = Xmax * (ab_vec[0] + ab_vec[2] * A + ab_vec[4] * B + F);
  xy[1] = Ymax * (ab_vec[1] + ab_vec[3] * A + ab_vec[5] * B + G);

  return xy;
}

//----------------------------------------------------------------
// RBFTransform::GetInverse
//
RBFTransform::InverseTransformPointer
RBFTransform::GetInverse() const
{
  unsigned int                 num_pts = num_points();
  std::vector<OutputPointType> xy_vec(num_pts);
  std::vector<InputPointType>  uv_vec(num_pts);

  for (unsigned int i = 0; i < num_pts; i++)
  {
    const double * uv = this->uv(i);
    uv_vec[i][0] = uv[0];
    uv_vec[i][1] = uv[1];

    // find where a given mosaic space control point maps to in the tile:
    xy_vec[i] = TransformPoint(uv_vec[i]);
  }

  // NOTE: this is not an exact inverse, some form of an iterative
  // inverse mapping calculation will be required (like find_inverse):
  OutputPointType tile_min;
  OutputPointType tile_max;
  tile_min[0] = 0;
  tile_min[1] = 0;
  tile_max[0] = GetXmax() * 2;
  tile_max[1] = GetYmax() * 2;

  RBFTransform::Pointer inverse = RBFTransform::New();
  inverse->setup(tile_min, tile_max, num_pts, num_pts ? &(xy_vec[0]) : nullptr, num_pts ? &(uv_vec[0]) : nullptr);

  return inverse.GetPointer();
}

//----------------------------------------------------------------
// RBFTransform::setup
//
void
RBFTransform::setup(                // image bounding box expressed in the image space,
                                    // defines transform normalization parameters:
  const OutputPointType & tile_min, // tile space
  const OutputPointType & tile_max, // tile space

  // landmark correspondences:
  const unsigned int      num_pts, // number of pairs
  const InputPointType *  uv,      // mosaic space
  const OutputPointType * xy)      // tile space
{
  this->m_FixedParameters.SetSize(4 + num_pts * 2);
  this->m_Parameters.SetSize(6 + num_pts * 2);

  // setup the normalization parameters:
  double & Xmax = this->m_FixedParameters[0];
  double & Ymax = this->m_FixedParameters[1];
  Xmax = (tile_max[0] - tile_min[0]) / 2.0;
  Ymax = (tile_max[1] - tile_min[1]) / 2.0;

  // store the control points and calculate the center of mass:
  double * uv_vec = &(this->m_FixedParameters[4]);
  double   u_sum = 0;
  double   v_sum = 0;
  for (unsigned int i = 0; i < num_pts; i++)
  {
    unsigned int offset = i * 2;
    uv_vec[offset] = uv[i][0];
    uv_vec[offset + 1] = uv[i][1];

    u_sum += uv[i][0];
    v_sum += uv[i][1];
  }

  // setup the center of rotation:
  double & uc = this->m_FixedParameters[2];
  double & vc = this->m_FixedParameters[3];
  uc = num_pts == 0 ? 0 : u_sum / double(num_pts);
  vc = num_pts == 0 ? 0 : v_sum / double(num_pts);

  if (num_pts == 0)
  {
    double * ab_vec = &(this->m_Parameters[index_a(0)]);
    ab_vec[0] = 0;
    ab_vec[1] = 0;
    ab_vec[2] = 1;
    ab_vec[3] = 0;
    ab_vec[4] = 0;
    ab_vec[5] = 1;
  }
  else
  {
    // setup the linear system to solve for transform parameters:
    vnl_matrix<double> M(num_pts + 3, num_pts + 3, 0);
    vnl_vector<double> bx(num_pts + 3, 0);
    vnl_vector<double> by(num_pts + 3, 0);

    for (unsigned int r = 0; r < num_pts; r++)
    {
      bx[r] = xy[r][0] / Xmax;
      by[r] = xy[r][1] / Ymax;

      // kernel components:
      for (unsigned int c = 0; c < num_pts; c++)
      {
        M(r, c) = kernel(&(uv[r][0]), &(uv[c][0]), Xmax, Ymax);
      }

      // polynomial components:
      M(r, num_pts) = 1;
      M(r, num_pts + 1) = (uv[r][0] - uc) / Xmax;
      M(r, num_pts + 2) = (uv[r][1] - vc) / Ymax;

      M(num_pts, r) = 1;
      M(num_pts + 1, r) = (uv[r][0] - uc) / Xmax;
      M(num_pts + 2, r) = (uv[r][1] - vc) / Ymax;
    }

    // FIXME:
#if 0
    cerr << "M: " << M << endl
	 << "bx: " << bx << endl
	 << "by: " << by << endl;
#endif

    // use SVD to solve the linear system:
    vnl_svd<double>    svd(M);
    vnl_vector<double> ax = svd.solve(bx);
    vnl_vector<double> ay = svd.solve(by);

    // store the polynomial coefficients into the parameter vector:
    double * ab_vec = &(this->m_Parameters[index_a(0)]);
    ab_vec[0] = ax[num_pts];
    ab_vec[1] = ay[num_pts];
    ab_vec[2] = ax[num_pts + 1];
    ab_vec[3] = ay[num_pts + 1];
    ab_vec[4] = ax[num_pts + 2];
    ab_vec[5] = ay[num_pts + 2];

    // store the Radial Basis Function coefficients and mosaic space
    // landmark coordinates into the transform parameter vector:
    double * fg_vec = &(this->m_Parameters[index_f(0)]);
    for (unsigned int i = 0; i < num_pts; i++)
    {
      unsigned int offset = i * 2;
      fg_vec[offset] = ax[i];
      fg_vec[offset + 1] = ay[i];
    }
  }

  // FIXME: PrintSelf(cerr, itk::Indent(0));
}

//----------------------------------------------------------------
// RBFTransform::GetJacobian
//
void
RBFTransform::ComputeJacobianWithRespectToParameters(const InputPointType & point,
                                                     JacobianType &         jacobian) const
{
  // shortcuts:
  const double & Xmax = GetXmax();
  const double & Ymax = GetYmax();
  const double & uc = GetUc();
  const double & vc = GetVc();

  const double & u = point[0];
  const double & v = point[1];

  const double A = (u - uc) / Xmax;
  const double B = (v - vc) / Ymax;

  jacobian(0, index_a(0)) = Xmax;
  jacobian(0, index_a(1)) = Xmax * A;
  jacobian(0, index_a(2)) = Xmax * B;

  jacobian(1, index_b(0)) = Ymax;
  jacobian(1, index_b(1)) = Ymax * A;
  jacobian(1, index_b(2)) = Ymax * B;

  const unsigned int num_pts = num_points();
  const double *     uv_vec = &(this->m_FixedParameters[index_uv(0)]);

  for (unsigned int i = 0; i < num_pts; i++)
  {
    unsigned int offset = i * 2;
    double       U = (u - uv_vec[offset]) / Xmax;
    double       V = (v - uv_vec[offset + 1]) / Ymax;

    double R = U * U + V * V;
    double lR = R == 0 ? 0 : log(R);
    double Q = R * lR;

    jacobian(0, index_f(i)) = Xmax * Q;
    jacobian(1, index_g(i)) = Ymax * Q;
  }
}

#if 0
//----------------------------------------------------------------
// RBFTransform::eval
// 
void
RBFTransform::eval(const std::vector<double> & point,
		   std::vector<double> & function,
		   std::vector<std::vector<double> > & jacobian) const
{
  // shortcuts:
  const double & Xmax = GetXmax();
  const double & Ymax = GetYmax();
  const double & uc = GetUc();
  const double & vc = GetVc();
  
  const double & u = point[0];
  const double & v = point[1];
  
  const double A = (u - uc) / Xmax;
  const double B = (v - vc) / Ymax;
  
  const unsigned int num_pts = num_points();
  
  double F = 0;
  double G = 0;
  double dF_du = 0;
  double dF_dv = 0;
  double dG_du = 0;
  double dG_dv = 0;
  
  for (unsigned int i = 0; i < num_pts; i++)
  {
    const double * fg = &(this->m_Parameters[index_f(i)]);
    const double & fi = fg[0];
    const double & gi = fg[1];
    
    const double * uv = &(this->m_FixedParameters[index_uv(i)]);
    const double & ui = uv[0];
    const double & vi = uv[1];
    
    double U = (u - ui) / Xmax;
    double V = (v - vi) / Ymax;
    double R = U * U + V * V;
    double r = sqrt(R);
    double lR = R == 0 ? 0 : log(R);
    double Q = R * lR;
    
    F += fi * Q;
    G += gi * Q;
    
    dF_du += fi * (u - ui) * (lnR + 1);
    dF_dv += fi * (v - vi) * (lnR + 1);
    dG_du += gi * (u - ui) * (lnR + 1);
    dG_dv += gi * (v - vi) * (lnR + 1);
  }
  
  dF_du *= 2 / (Xmax * Xmax);
  dF_dv *= 2 / (Ymax * Ymax);
  dG_du *= 2 / (Xmax * Xmax);
  dG_dv *= 2 / (Ymax * Ymax);
  
  function[0] = Xmax * (a(0) + a(1) * A + a(2) * B + F);
  function[1] = Ymax * (b(0) + b(1) * A + b(2) * B + G);
  
  // dx/du:
  jacobian[0][0] = Xmax * (a(1) + dF_du);
  
  // dx/dv:
  jacobian[0][1] = Xmax * (a(2) + dF_dv);
  
  // dy/du:
  jacobian[1][0] = Ymax * (b(1) + dG_du);
  
  // dy/dv:
  jacobian[1][1] = Ymax * (b(2) + dG_dv);
}
#endif

//----------------------------------------------------------------
// RBFTransform::PrintSelf
//
void
RBFTransform::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  unsigned int num_pts = num_points();

  os << indent << "Xmax = " << GetXmax() << endl
     << indent << "Ymax = " << GetYmax() << endl
     << indent << "uc = " << GetUc() << endl
     << indent << "vc = " << GetVc() << endl;

  os << indent << "a0 = " << a(0) << endl << indent << "a1 = " << a(1) << endl << indent << "a2 = " << a(2) << endl;

  for (unsigned int i = 0; i < num_pts; i++)
  {
    os << indent << 'f' << i << " = " << f(i) << endl;
  }

  os << indent << "b0 = " << b(0) << endl << indent << "b1 = " << b(1) << endl << indent << "b2 = " << b(2) << endl;

  for (unsigned int i = 0; i < num_pts; i++)
  {
    os << indent << 'g' << i << " = " << g(i) << endl;
  }
}


#if 0
{
#endif
} // namespace itk
