// -*- Mode: c++; tab-width: 8; c-basic-offset: 2; indent-tabs-mode: nil -*-
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


// File         : the_grid_transform.hxx
// Author       : Pavel A. Koshevoy
// Created      : Thu Nov 30 10:45:47 MST 2006
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : A discontinuous transform -- a uniform grid of vertices is
//                mapped to an image. At each vertex, in addition to image
//                space coordinates, a second set of coordinates is stored.
//                This is similar to texture mapped OpenGL triangle meshes,
//                where the texture coordinates correspond to the image space
//                vertex coordinates.

#ifndef THE_GRID_TRANSFORM_HXX_
#define THE_GRID_TRANSFORM_HXX_

// local includes:
#include "itkIRCommon.h"

// system includes:
#include <vector>
#include <list>


//----------------------------------------------------------------
// vertex_t
//
class vertex_t
{
public:
  // normalized tile space coordinates, typically [0, 1] x [0, 1]:
  pnt2d_t uv_;

  // physical space coordinates:
  pnt2d_t xy_;
};

//----------------------------------------------------------------
// triangle_t
//
class triangle_t
{
public:
  triangle_t();

  // check whether a given xy-point falls within this triangle,
  // return corresponding uv-point (not triangle barycentric coordinates):
  bool
  xy_intersect(const vertex_t * v_arr, const pnt2d_t & xy, pnt2d_t & uv) const;

  // check whether a given uv-point falls within this triangle,
  // return corresponding xy-point  (not triangle barycentric coordinates):
  bool
  uv_intersect(const vertex_t * v_arr, const pnt2d_t & uv, pnt2d_t & xy) const;

  // triangle vertex indices, counterclockwise winding:
  unsigned int vertex_[3];

  // precomputed fast barycentric coordinate calculation coefficients:

  // for intersection calculation in xy-space:
  double xy_pwb[3];
  double xy_pwc[3];

  // for intersection calculation in uv-space:
  double uv_pwb[3];
  double uv_pwc[3];
};

//----------------------------------------------------------------
// the_acceleration_grid_t
//
// The bounding grid triangle/point intersection acceleration
// structure used to speed up grid transform and mesh transform:
//
class the_acceleration_grid_t
{
public:
  the_acceleration_grid_t();

  // find the grid cell containing a given xy-point:
  unsigned int
  xy_cell(const pnt2d_t & xy) const;

  // find the triangle containing a given xy-point,
  // and calculate corresponding uv-point:
  unsigned int
  xy_triangle(const pnt2d_t & xy, pnt2d_t & uv) const;

  // find the grid cell containing a given uv-point:
  unsigned int
  uv_cell(const pnt2d_t & uv) const;

  // find the triangle containing a given uv-point,
  // and calculate corresponding xy-point:
  unsigned int
  uv_triangle(const pnt2d_t & uv, pnt2d_t & xy) const;

  // update the vertex xy coordinates and rebuild the grid
  void
  update(const vec2d_t * xy_shift);
  void
  shift(const vec2d_t & xy_shift);

  // resize the grid:
  void
  resize(unsigned int rows, unsigned int cols);

  // rebuild the acceleration grid:
  void
  rebuild();

private:
  // helper used to rebuild the grid:
  void
  update_grid(unsigned int t_idx);

public:
  // the acceleration structure:
  std::vector<std::vector<unsigned int>> xy_;
  std::vector<std::vector<unsigned int>> uv_;
  unsigned int                           rows_;
  unsigned int                           cols_;

  // the grid bounding box (in xy-space):
  pnt2d_t xy_min_;
  vec2d_t xy_ext_;

  // the triangle mesh:
  std::vector<vertex_t>   mesh_;
  std::vector<triangle_t> tri_;
};


//----------------------------------------------------------------
// the_base_triangle_transform_t
//
class the_base_triangle_transform_t
{
public:
  the_base_triangle_transform_t() {}

  // transform the point:
  bool
  transform(const pnt2d_t & xy, pnt2d_t & uv) const;

  // inverse transform the point:
  bool
  transform_inv(const pnt2d_t & uv, pnt2d_t & xy) const;

  // calculate the derivatives of the transforms with respect to
  // transform parameters:
  bool
  jacobian(const pnt2d_t & xy, unsigned int * idx, double * jac) const;

public:
  // tile bounding box:
  pnt2d_t tile_min_;
  vec2d_t tile_ext_;

  // the acceleration grid (stores triangle vertices, and triangles):
  the_acceleration_grid_t grid_;
};


//----------------------------------------------------------------
// the_grid_transform_t
//
class the_grid_transform_t : public the_base_triangle_transform_t
{
public:
  the_grid_transform_t();

  // check to see whether the transform has already been setup:
  bool
  is_ready() const;

  // vertex accessors:
  inline const vertex_t &
  vertex(size_t row, size_t col) const
  {
    return grid_.mesh_[row * (cols_ + 1) + col];
  }

  inline vertex_t &
  vertex(size_t row, size_t col)
  {
    return grid_.mesh_[row * (cols_ + 1) + col];
  }

  // inverse transform the point:
  bool
  transform_inv(const pnt2d_t & uv, pnt2d_t & xy) const;

  // setup the transform:
  void
  setup(unsigned int                 rows,
        unsigned int                 cols,
        const pnt2d_t &              tile_min,
        const pnt2d_t &              tile_max,
        const std::vector<pnt2d_t> & xy);

private:
  // helper used to setup the triangle mesh:
  void
  setup_mesh();

public:
  // number of rows and columns of quads in the mesh
  // (each quad is made up of 2 triangles):
  size_t rows_;
  size_t cols_;
};


//----------------------------------------------------------------
// the_mesh_transform_t
//
class the_mesh_transform_t : public the_base_triangle_transform_t
{
public:
  // check to see whether the transform has already been setup:
  bool
  is_ready() const;

  // setup the transform:
  bool
  setup(const pnt2d_t &              tile_min,
        const pnt2d_t &              tile_max,
        const std::vector<pnt2d_t> & uv,
        const std::vector<pnt2d_t> & xy,
        unsigned int                 accel_grid_rows = 16,
        unsigned int                 accel_grid_cols = 16);

  // insert a point into the mesh, and re-triangulate
  // using Delaunay triangulation:
  bool
  insert_point(const pnt2d_t & uv, const pnt2d_t & xy, const bool delay_setup = false);

  // insert a point into the mesh (xy-point is extrapolated),
  // and re-triangulate using Delaunay triangulation:
  bool
  insert_point(const pnt2d_t & uv);

private:
  // helper used to setup the triangle mesh:
  bool
  setup_mesh();
};


#endif // THE_GRID_TRANSFORM_HXX_
