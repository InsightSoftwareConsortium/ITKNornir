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


// File         : mosaic_refinement_common.hxx
// Author       : Pavel A. Koshevoy, Joel Spaltenstein
// Created      : Mon Mar 26 10:34:39 MDT 2007
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : Helper functions for automatic mosaic refinement.

#ifndef MOSAIC_REFINEMENT_COMMON_HXX_
#define MOSAIC_REFINEMENT_COMMON_HXX_

// the includes:
#include <itkIRCommon.h>
#include <IRGridCommon.h>
#include <itkImageMosaicVarianceMetric.h>
#include <itkRegularStepGradientDescentOptimizer2.h>
#include <IRAABBox.h>
#include <IRThreadPool.h>

#include "itkImageDuplicator.h"

// system includes:
#include <math.h>
#include <vector>
#include <sstream>

#ifndef WIN32
#  include <unistd.h>
#endif


//----------------------------------------------------------------
// regularize_displacements
//
extern void
regularize_displacements( // computed displacement vectors of the moving image
                          // grid transform control points, in mosaic space:
  std::vector<vec2d_t> & xy_shift,
  std::vector<double> &  mass,

  image_t::Pointer & dx,
  image_t::Pointer & dy,
  image_t::Pointer & db,

  // median filter radius:
  const unsigned int & median_radius);


//----------------------------------------------------------------
// refine_mosaic
//
template <class TImage, class TMask>
void
refine_mosaic(the_log_t & log,

              // all the tiles in the mosaic, and their masks:
              array2d(typename TImage::Pointer) & pyramid,
              std::vector<typename TMask::ConstPointer> & mask,
              std::vector<base_transform_t::Pointer> &    transform,

              unsigned int iterations_per_level)
{
  typedef itk::LinearInterpolateImageFunction<TImage, double> interpolator_t;

  unsigned int pyramid_levels = pyramid.size();
  if (pyramid_levels == 0)
    return;

  std::vector<typename TImage::Pointer> & image = pyramid[pyramid_levels - 1];
  unsigned int                            num_images = image.size();
  if (num_images == 0)
    return;

  log << "iterations per level: " << iterations_per_level << endl
      << "transform type: " << transform[0]->GetTransformTypeAsString() << endl;

  // try global refinement of the mosaic:
  typedef itk::ImageMosaicVarianceMetric<TImage, interpolator_t> mosaic_metric_t;

  typename mosaic_metric_t::Pointer mosaic_metric = mosaic_metric_t::New();
  mosaic_metric->image_.resize(num_images);
  mosaic_metric->mask_.resize(num_images);
  mosaic_metric->transform_.resize(num_images);
  for (unsigned int i = 0; i < num_images; i++)
  {
    mosaic_metric->image_[i] = image[i];
    mosaic_metric->mask_[i] = mask[i];
    mosaic_metric->transform_[i] = transform[i];
  }

  // FIXME: ITK doesn't have an API for this:
  std::vector<bool> param_shared(transform[0]->GetNumberOfParameters(), false);
  std::vector<bool> param_active(transform[0]->GetNumberOfParameters(), true);

  // setup the shared parameters mask:
  mosaic_metric->setup_param_map(param_shared, param_active);
  mosaic_metric->Initialize();

  // setup the optimizer scales:
  typename mosaic_metric_t::params_t parameter_scales = mosaic_metric->GetTransformParameters();
  parameter_scales.Fill(1.0);

  for (unsigned int level = 0; level < pyramid_levels && iterations_per_level > 0; level++)
  {
    for (unsigned int i = 0; i < num_images; i++)
    {
      mosaic_metric->image_[i] = pyramid[level][i];
    }

    typename mosaic_metric_t::measure_t metric_before =
      mosaic_metric->GetValue(mosaic_metric->GetTransformParameters());

    // run several iterations of the optimizer:
    for (unsigned int k = 0; k < 3; k++)
    {
      typename mosaic_metric_t::params_t params_before = mosaic_metric->GetTransformParameters();

      typename mosaic_metric_t::measure_t metric_after =
        std::numeric_limits<typename mosaic_metric_t::measure_t>::max();

      // use global refinement:
      optimizer_t::Pointer                                optimizer = optimizer_t::New();
      typename optimizer_observer_t<optimizer_t>::Pointer observer = optimizer_observer_t<optimizer_t>::New();
      observer->log_ = &log;
      optimizer->AddObserver(itk::IterationEvent(), observer);
      optimizer->SetLog(&log);
      optimizer->SetMinimize(true);
      optimizer->SetNumberOfIterations(iterations_per_level);
      optimizer->SetMinimumStepLength(1e-12); // min_step
      optimizer->SetMaximumStepLength(1e-5);  // max_step
      optimizer->SetGradientMagnitudeTolerance(1e-6);
      optimizer->SetRelaxationFactor(5e-1);
      optimizer->SetCostFunction(mosaic_metric);
      optimizer->SetInitialPosition(params_before);
      optimizer->SetScales(parameter_scales);
      optimizer->SetPickUpPaceSteps(5);
      optimizer->SetBackTracking(true);

      // refine the mosaic:
      try
      {
        log << "\n" << level << '.' << k << ": refining distortion transforms" << endl;

        optimizer->StartOptimization();
      }
      catch (itk::ExceptionObject & exception)
      {
        // oops:
        log << "optimizer threw an exception:" << endl << exception.what() << endl;
      }

      mosaic_metric->SetTransformParameters(optimizer->GetBestParams());
      metric_after = optimizer->GetBestValue();

      typename mosaic_metric_t::params_t params_after = mosaic_metric->GetTransformParameters();

      log << "before: METRIC = " << metric_before << ", PARAMS = " << params_before << endl
          << "after:  METRIC = " << metric_after << ", PARAMS = " << params_after << endl;

      // quantify the improvement:
      double improvement = 1.0 - metric_after / metric_before;
      bool   failed_to_improve = (metric_after - metric_before) >= 0.0;
      bool   negligible_improvement = !failed_to_improve && (improvement < 1e-3);

      if (!failed_to_improve)
      {
        log << "IMPROVEMENT: " << setw(3) << int(100.0 * improvement) << "%" << endl;
      }

      if (failed_to_improve)
      {
        log << "NOTE: minimization failed, ignoring registration results..." << endl;

        // previous transform was better:
        mosaic_metric->SetTransformParameters(params_before);
        break;
      }
      else if (negligible_improvement)
      {
        log << "NOTE: improvement is negligible..." << endl;
        break;
      }

      // avoid recalculating the same metric:
      metric_before = metric_after;
    }
  }
}


//----------------------------------------------------------------
// calc_displacements
//
// This is the original version used for single-threaded execution
//
template <typename TImage, typename TMask>
void
calc_displacements(the_log_t & log,

                   // computed displacement vectors of the moving image
                   // grid transform control points, in mosaic space:
                   std::vector<vec2d_t> & xy_shift,
                   std::vector<double> &  mass,

                   // a flag indicating whether the tiles
                   // have already been transformed into mosaic space:
                   bool tiles_already_warped,

                   // fixed:
                   const TImage *           tile_0,
                   const TMask *            mask_0,
                   const base_transform_t * forward_0,

                   // moving:
                   const TImage *             tile_1,
                   const TMask *              mask_1,
                   const itk::GridTransform * forward_1,

                   // neighborhood size:
                   const unsigned int & neighborhood,

                   // minimum acceptable neighborhood overlap ratio:
                   const double & min_overlap,

                   // median filter radius:
                   const unsigned int & median_radius)
{
  // shortcuts:
  const the_grid_transform_t & gt = forward_1->transform_;
  unsigned int                 mesh_cols = gt.cols_ + 1;
  unsigned int                 mesh_rows = gt.rows_ + 1;
  unsigned                     mesh_size = gt.grid_.mesh_.size();
  xy_shift.assign(mesh_size, vec2d(0, 0));

  // make sure both tiles have the same pixel spacing:
  typename TImage::SpacingType sp = tile_1->GetSpacing();
  if (sp != tile_0->GetSpacing())
    return;

  // setup the local neighborhood:
  typename TImage::SizeType sz;
  sz[0] = neighborhood;
  sz[1] = neighborhood;

  typename TImage::Pointer img[] = { make_image<TImage>(sp, sz), make_image<TImage>(sp, sz) };

  typename TMask::Pointer msk[] = { make_image<TMask>(sp, sz), make_image<TMask>(sp, sz) };

  // for each interpolation point, do a local neighborhood fft matching,
  // and use the resulting displacement vector to adjust the mesh:
  image_t::Pointer dx = make_image<image_t>(mesh_cols, mesh_rows, 1.0, 0.0);

  image_t::Pointer dy = make_image<image_t>(mesh_cols, mesh_rows, 1.0, 0.0);

  image_t::Pointer db = make_image<image_t>(mesh_cols, mesh_rows, 1.0, 0.0);

  typename TImage::Pointer img_large;
  typename TMask::Pointer  msk_large;

  if (!tiles_already_warped)
  {
    typename TImage::SizeType sz_large(sz);
    sz_large[0] *= 2;
    sz_large[1] *= 2;
    img_large = make_image<TImage>(sp, sz_large);
    msk_large = make_image<TMask>(sp, sz_large);
  }

  log << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
#pragma omp parallel for
  for (int i = 0; i < (int)mesh_size; i++)
  {
    // shortcut:
    const vertex_t & vertex = gt.grid_.mesh_[i];

    // find the mosaic space coordinates of this vertex:
    pnt2d_t center;
    gt.transform_inv(vertex.uv_, center);

    // extract a neighborhood of the vertex from both tiles:
    image_t::IndexType index;
    index[0] = i % mesh_cols;
    index[1] = i / mesh_cols;
    dx->SetPixel(index, 0);
    dy->SetPixel(index, 0);
    db->SetPixel(index, 0);

    // feed the two neighborhoods into the FFT translation estimator:
    vec2d_t shift(vec2d(0, 0));
    bool    ok = tiles_already_warped ?

                                      refine_one_point_fft(log,
                                                        shift,

                                                        tile_0,
                                                        mask_0,

                                                        tile_1,
                                                        mask_1,

                                                        center,
                                                        min_overlap,

                                                        img[0].GetPointer(),
                                                        msk[0].GetPointer(),

                                                        img[1].GetPointer(),
                                                        msk[1].GetPointer())
                                      :

                                      refine_one_point_fft(log,
                                                        shift,

                                                        tile_0,
                                                        mask_0,

                                                        tile_1,
                                                        mask_1,

                                                        forward_0,
                                                        forward_1,

                                                        center,
                                                        min_overlap,

                                                        sz,
                                                        sp,

                                                        img_large.GetPointer(),
                                                        msk_large.GetPointer(),

                                                        img[0].GetPointer(),
                                                        msk[0].GetPointer(),

                                                        img[1].GetPointer(),
                                                        msk[1].GetPointer());
    if (!ok)
    {
      continue;
    }

    log << i << ". shift: " << shift << endl;
    dx->SetPixel(index, shift[0]);
    dy->SetPixel(index, shift[1]);
    db->SetPixel(index, 1);
  }

  // regularize the displacement vectors here:
  regularize_displacements(xy_shift, mass, dx, dy, db, median_radius);
}


//----------------------------------------------------------------
// calc_displacements_t
//
//
template <typename TImage, typename TMask>
class calc_displacements_t : public the_transaction_t
{
public:
  calc_displacements_t( // a flag indicating whether the tiles
                        // have already been transformed into mosaic space:
    bool tiles_already_warped,

    // fixed:
    const TImage *           tile_0,
    const TMask *            mask_0,
    const base_transform_t * forward_0,

    // moving:
    const TImage *             tile_1,
    const TMask *              mask_1,
    const itk::GridTransform * forward_1,

    // neighborhood size:
    const unsigned int & neighborhood_size,

    // minimum acceptable neighborhood overlap ratio:
    const double & min_overlap,

    // mesh node displacements:
    image_t::Pointer dx,
    image_t::Pointer dy,

    // mesh node displacement weights:
    image_t::Pointer db,

    // mesh node index:
    const std::list<image_t::IndexType> & index,

    // mesh node mosaic space coordinates:
    const std::list<pnt2d_t> & center)
    : tiles_already_warped_(tiles_already_warped)
    ,

    tile_0_(tile_0)
    , mask_0_(mask_0)
    , forward_0_(forward_0)
    ,

    tile_1_(tile_1)
    , mask_1_(mask_1)
    , forward_1_(forward_1)
    ,

    min_overlap_(min_overlap)
    ,

    dx_(dx)
    , dy_(dy)
    , db_(db)
    ,

    index_(index.begin(), index.end())
    , center_(center.begin(), center.end())
  {
    // make sure both tiles have the same pixel spacing:
    sp_ = tile_1->GetSpacing();
    bool ok = (sp_ == tile_0->GetSpacing());
    if (!ok)
    {
      assert(false);
      return;
    }

    // setup the local neighborhood:
    sz_[0] = neighborhood_size;
    sz_[1] = neighborhood_size;
  }

  // virtual:
  void
  execute(the_thread_interface_t * thread)
  {
    WRAP(the_terminator_t terminator("calc_displacements_t"));

    typename TImage::Pointer img[] = { make_image<TImage>(sp_, sz_), make_image<TImage>(sp_, sz_) };

    typename TMask::Pointer msk[] = { make_image<TMask>(sp_, sz_), make_image<TMask>(sp_, sz_) };

    typename TImage::Pointer img_large;
    typename TMask::Pointer  msk_large;

    if (!tiles_already_warped_)
    {
      typename TImage::SizeType sz_large(sz_);
      sz_large[0] *= 2;
      sz_large[1] *= 2;
      img_large = make_image<TImage>(sp_, sz_large);
      msk_large = make_image<TMask>(sp_, sz_large);
    }

    std::size_t num_nodes = center_.size();
    for (std::size_t i = 0; i < num_nodes; i++)
    {
      // shortcuts:
      const image_t::IndexType & index = index_[i];
      const pnt2d_t &            center = center_[i];

      // feed the two neighborhoods into the FFT translation estimator:
      vec2d_t shift(vec2d(0, 0));

      bool ok = tiles_already_warped_ ?
                                      // if the tiles are already warped we don't
                                      // need the transforms
                  refine_one_point_fft(*null_log(),
                                       shift,

                                       tile_0_,
                                       mask_0_,

                                       tile_1_,
                                       mask_1_,

                                       center,
                                       min_overlap_,

                                       img[0].GetPointer(),
                                       msk[0].GetPointer(),

                                       img[1].GetPointer(),
                                       msk[1].GetPointer())
                                      :

                                      refine_one_point_fft(*null_log(),
                                                           shift,

                                                           tile_0_,
                                                           mask_0_,

                                                           tile_1_,
                                                           mask_1_,

                                                           forward_0_,
                                                           forward_1_,

                                                           center,
                                                           min_overlap_,

                                                           sz_,
                                                           sp_,

                                                           img_large.GetPointer(),
                                                           msk_large.GetPointer(),

                                                           img[0].GetPointer(),
                                                           msk[0].GetPointer(),

                                                           img[1].GetPointer(),
                                                           msk[1].GetPointer());
      if (ok)
      {
        dx_->SetPixel(index, shift[0]);
        dy_->SetPixel(index, shift[1]);
        db_->SetPixel(index, 1);
      }
    }
  }

  // a flag indicating whether the tiles
  // have already been transformed into mosaic space:
  bool tiles_already_warped_;

  // fixed:
  const TImage *           tile_0_;
  const TMask *            mask_0_;
  const base_transform_t * forward_0_;

  // moving:
  const TImage *             tile_1_;
  const TMask *              mask_1_;
  const itk::GridTransform * forward_1_;

  // moving tile spacing (must be the same as fixed tile spacing):
  typename TImage::SpacingType sp_;
  typename TImage::SizeType    sz_;

  // minimum acceptable neighborhood overlap ratio:
  const double min_overlap_;

  // mesh node displacements:
  image_t::Pointer dx_;
  image_t::Pointer dy_;

  // mesh node displacement weights:
  image_t::Pointer db_;

  // mesh node index:
  std::vector<image_t::IndexType> index_;

  // mesh node mosaic space coordinates:
  std::vector<pnt2d_t> center_;
};


//----------------------------------------------------------------
// calc_displacements_mt
//
// calculate transform mesh displacement vectors multi-threaded
//
template <typename TImage, typename TMask>
void
calc_displacements_mt(unsigned int num_threads,

                      the_log_t & log,

                      // computed displacement vectors of the moving image
                      // grid transform control points, in mosaic space:
                      std::vector<vec2d_t> & xy_shift,
                      std::vector<double> &  mass,

                      // a flag indicating whether the tiles
                      // have already been transformed into mosaic space:
                      bool tiles_already_warped,

                      // fixed:
                      typename TImage::ConstPointer           tile_0,
                      typename TMask::ConstPointer            mask_0,
                      typename base_transform_t::ConstPointer forward_0,

                      // moving:
                      typename TImage::ConstPointer               tile_1,
                      typename TMask::ConstPointer                mask_1,
                      typename itk::GridTransform::ConstPointer & forward_1,

                      // neighborhood size:
                      const unsigned int & neighborhood_size,

                      // minimum acceptable neighborhood overlap ratio:
                      const double & min_overlap,

                      // median filter radius:
                      const unsigned int & median_radius)
{
  // make sure both tiles have the same pixel spacing:
  if (tile_1->GetSpacing() != tile_0->GetSpacing())
    return;

  // shortcuts:
  const the_grid_transform_t & gt = forward_1->transform_;
  unsigned int                 mesh_cols = gt.cols_ + 1;
  unsigned int                 mesh_rows = gt.rows_ + 1;
  unsigned                     mesh_size = gt.grid_.mesh_.size();
  xy_shift.assign(mesh_size, vec2d(0, 0));

  // for each interpolation point, do a local neighborhood fft matching,
  // and use the resulting displacement vector to adjust the mesh:
  image_t::Pointer dx = make_image<image_t>(mesh_cols, mesh_rows, 1.0, 0.0);

  image_t::Pointer dy = make_image<image_t>(mesh_cols, mesh_rows, 1.0, 0.0);

  image_t::Pointer db = make_image<image_t>(mesh_cols, mesh_rows, 1.0, 0.0);

  the_thread_pool_t thread_pool(num_threads);
  thread_pool.set_idle_sleep_duration(50); // 50 usec

  // split nodes between threads:
  std::vector<std::list<image_t::IndexType>> node_index_list(num_threads);
  std::vector<std::list<pnt2d_t>>            node_center_list(num_threads);

#pragma omp parallel for
  for (int i = 0; i < mesh_size; i++)
  {
    // shortcuts:
    const vertex_t &   vertex = gt.grid_.mesh_[i];
    const unsigned int which_thread = i % num_threads;

    // find the mosaic space coordinates of this vertex:
    pnt2d_t center;
    gt.transform_inv(vertex.uv_, center);
    node_center_list[which_thread].push_back(center);

    // extract a neighborhood of the vertex from both tiles:
    image_t::IndexType index;
    index[0] = i % mesh_cols;
    index[1] = i / mesh_cols;
    dx->SetPixel(index, 0);
    dy->SetPixel(index, 0);
    db->SetPixel(index, 0);

    node_index_list[which_thread].push_back(index);
  }

  // setup a transaction for each thread:
  for (unsigned int i = 0; i < num_threads; i++)
  {
    calc_displacements_t<TImage, TMask> * t = new calc_displacements_t<TImage, TMask>(tiles_already_warped,

                                                                                      tile_0,
                                                                                      mask_0,
                                                                                      forward_0,

                                                                                      tile_1,
                                                                                      mask_1,
                                                                                      forward_1,

                                                                                      neighborhood_size,
                                                                                      min_overlap,

                                                                                      dx,
                                                                                      dy,
                                                                                      db,

                                                                                      node_index_list[i],
                                                                                      node_center_list[i]);

    thread_pool.push_back(t);
  }

  // run the transactions:
  thread_pool.pre_distribute_work();
  suspend_itk_multithreading_t suspend_itk_mt;
  thread_pool.start();
  thread_pool.wait();

  // regularize the displacement vectors here:
  regularize_displacements(xy_shift, mass, dx, dy, db, median_radius);
}


//----------------------------------------------------------------
// intermediate_result_t
//
class intermediate_result_t
{
public:
  intermediate_result_t() {}

  intermediate_result_t(const intermediate_result_t & r) { *this = r; }

  intermediate_result_t(unsigned int num_neighbors, unsigned int mesh_rows, unsigned int mesh_cols)
    : dx_(num_neighbors)
    , dy_(num_neighbors)
    , db_(num_neighbors)
  {
    for (unsigned int i = 0; i < num_neighbors; i++)
    {
      dx_[i] = make_image<image_t>(mesh_cols, mesh_rows, 1.0, 0.0);

      dy_[i] = make_image<image_t>(mesh_cols, mesh_rows, 1.0, 0.0);

      db_[i] = make_image<image_t>(mesh_cols, mesh_rows, 1.0, 0.0);
    }
  }

  intermediate_result_t &
  operator=(const intermediate_result_t & r)
  {
    if (this != &r)
    {
      dx_.resize(r.dx_.size());
      for (unsigned int i = 0; i < dx_.size(); i++)
      {
        dx_[i] = cast<image_t, image_t>(r.dx_[i]);
      }

      dy_.resize(r.dy_.size());
      for (unsigned int i = 0; i < dy_.size(); i++)
      {
        dy_[i] = cast<image_t, image_t>(r.dy_[i]);
      }

      db_.resize(r.db_.size());
      for (unsigned int i = 0; i < db_.size(); i++)
      {
        db_[i] = cast<image_t, image_t>(r.db_[i]);
      }
    }

    return *this;
  }

  std::vector<image_t::Pointer> dx_;
  std::vector<image_t::Pointer> dy_;
  std::vector<image_t::Pointer> db_;
};

//----------------------------------------------------------------
// calc_intermediate_results_t
//
template <typename img_ptr_t, typename msk_ptr_t>
class calc_intermediate_results_t : public the_transaction_t
{
public:
  typedef typename img_ptr_t::ObjectType::Pointer::ObjectType TImage;
  typedef typename msk_ptr_t::ObjectType::Pointer::ObjectType TMask;

  calc_intermediate_results_t(the_log_t &                                            log,
                              unsigned int                                           thread_offset,
                              unsigned int                                           thread_stride,
                              std::vector<itk::GridTransform::Pointer> &             transform,
                              const std::vector<img_ptr_t> &                         warped_tile,
                              const std::vector<msk_ptr_t> &                         warped_mask,
                              const std::vector<the_dynamic_array_t<unsigned int>> & neighbors,
                              const bool &                                           tiles_already_warped,
                              const unsigned int &                                   neighborhood_size,
                              const double &                                         minimum_overlap,
                              const bool &                                           keep_first_tile_fixed,
                              std::vector<intermediate_result_t> &                   results)
    :

    log_(log)
    , thread_offset_(thread_offset)
    , thread_stride_(thread_stride)
    , transform_(transform)
    , warped_tile_(warped_tile)
    , warped_mask_(warped_mask)
    , neighbors_(neighbors)
    , tiles_already_warped_(tiles_already_warped)
    , neighborhood_(neighborhood_size)
    , minimum_overlap_(minimum_overlap)
    , keep_first_tile_fixed_(keep_first_tile_fixed)
    , results_(results)
  {}

  // virtual:
  void
  execute(the_thread_interface_t * thread)
  {
    WRAP(the_terminator_t terminator("calc_intermediate_results_t::execute"));

    unsigned int num_tiles = warped_tile_.size();
    unsigned int start = keep_first_tile_fixed_ ? 1 : 0;

    // setup the local neighborhood image buffers:
    typename TImage::SpacingType sp = warped_tile_[start]->GetSpacing();
    typename TImage::SizeType    sz;
    sz[0] = neighborhood_;
    sz[1] = neighborhood_;

    typename TImage::Pointer img[] = { make_image<TImage>(sp, sz), make_image<TImage>(sp, sz) };

    typename TMask::Pointer msk[] = { make_image<TMask>(sp, sz), make_image<TMask>(sp, sz) };

    typename TImage::Pointer img_large;
    typename TMask::Pointer  msk_large;

    if (!tiles_already_warped_)
    {
      typename TImage::SizeType sz_large(sz);
      sz_large[0] *= 2;
      sz_large[1] *= 2;
      img_large = make_image<TImage>(sp, sz_large);
      msk_large = make_image<TMask>(sp, sz_large);
    }

    for (unsigned int tile_index = start; tile_index < num_tiles; tile_index++)
    {
      // shortcuts:
      const unsigned int         num_neighbors = neighbors_[tile_index].size();
      const TImage *             tile = warped_tile_[tile_index];
      const TMask *              mask = warped_mask_[tile_index];
      const itk::GridTransform * transform = transform_[tile_index];

      const the_grid_transform_t & gt = transform->transform_;
      const unsigned int           mesh_cols = gt.cols_ + 1;
      const unsigned int           mesh_size = gt.grid_.mesh_.size();

      // shortcuts to intermediate mesh refinement results for this tile:
      intermediate_result_t & results = results_[tile_index];

      for (unsigned int neighbor = 0; neighbor < num_neighbors; neighbor++)
      {
        WRAP(terminator.terminate_on_request());

        // shortcuts:
        image_t::Pointer & dx = results.dx_[neighbor];
        image_t::Pointer & dy = results.dy_[neighbor];
        image_t::Pointer & db = results.db_[neighbor];

        const unsigned int neighbor_index = neighbors_[tile_index][neighbor];
        log_ << thread_offset_ << " thread, matching " << tile_index << ":" << neighbor_index << endl;

        // shortcuts:
        const TImage *           neighbor_tile = warped_tile_[neighbor_index];
        const TMask *            neighbor_mask = warped_mask_[neighbor_index];
        const base_transform_t * neighbor_xform = transform_[neighbor_index];

        for (unsigned int mesh_index = thread_offset_; mesh_index < mesh_size; mesh_index += thread_stride_)
        {
          // shortcut:
          const vertex_t & vertex = gt.grid_.mesh_[mesh_index];

          // find the mosaic space coordinates of this vertex:
          pnt2d_t center;
          gt.transform_inv(vertex.uv_, center);

          // extract a neighborhood of the vertex from both tiles:
          image_t::IndexType index;
          index[0] = mesh_index % mesh_cols;
          index[1] = mesh_index / mesh_cols;
          dx->SetPixel(index, 0);
          dy->SetPixel(index, 0);
          db->SetPixel(index, 0);

          // feed the two neighborhoods into the FFT translation estimator:
          vec2d_t shift(vec2d(0, 0));
          bool    ok = tiles_already_warped_ ?

                                             refine_one_point_fft(*null_log(),
                                                               shift,

                                                               // fixed:
                                                               neighbor_tile,
                                                               neighbor_mask,

                                                               // moving:
                                                               tile,
                                                               mask,

                                                               center,
                                                               minimum_overlap_,

                                                               img[0].GetPointer(),
                                                               msk[0].GetPointer(),

                                                               img[1].GetPointer(),
                                                               msk[1].GetPointer())
                                             :

                                             refine_one_point_fft(*null_log(),
                                                               shift,

                                                               // fixed:
                                                               neighbor_tile,
                                                               neighbor_mask,

                                                               // moving:
                                                               tile,
                                                               mask,

                                                               neighbor_xform,
                                                               transform,

                                                               center,
                                                               minimum_overlap_,

                                                               sz,
                                                               sp,

                                                               img_large.GetPointer(),
                                                               msk_large.GetPointer(),

                                                               img[0].GetPointer(),
                                                               msk[0].GetPointer(),

                                                               img[1].GetPointer(),
                                                               msk[1].GetPointer());
          if (!ok)
          {
            continue;
          }

          dx->SetPixel(index, shift[0]);
          dy->SetPixel(index, shift[1]);
          db->SetPixel(index, 1);
        }
      }
    }
  }

  the_log_t &                                            log_;
  unsigned int                                           thread_offset_;
  unsigned int                                           thread_stride_;
  std::vector<itk::GridTransform::Pointer> &             transform_;
  const std::vector<img_ptr_t> &                         warped_tile_;
  const std::vector<msk_ptr_t> &                         warped_mask_;
  const std::vector<the_dynamic_array_t<unsigned int>> & neighbors_;
  const bool                                             tiles_already_warped_;
  const unsigned int                                     neighborhood_;
  const double                                           minimum_overlap_; // neighborhood overlap
  const bool                                             keep_first_tile_fixed_;
  std::vector<intermediate_result_t> &                   results_;
};


//----------------------------------------------------------------
// update_tile_mesh_t
//
class update_tile_mesh_t : public the_transaction_t
{
public:
  update_tile_mesh_t(the_log_t &                                log,
                     unsigned int                               tile_index,
                     bool                                       keep_first_tile_fixed,
                     unsigned int                               median_filter_radius,
                     std::vector<itk::GridTransform::Pointer> & transform,
                     std::vector<intermediate_result_t> &       results)
    : log_(log)
    , tile_index_(tile_index)
    , keep_first_tile_fixed_(keep_first_tile_fixed)
    , median_filter_radius_(median_filter_radius)
    , transform_(transform)
    , results_(results)
  {}

  // virtual:
  void
  execute(the_thread_interface_t * thread)
  {
    WRAP(the_terminator_t terminator("update_tile_mesh_t::execute"));

    log_ << tile_index_ << " mesh update" << endl;

    // shortcuts:
    itk::GridTransform *   transform = transform_[tile_index_];
    the_grid_transform_t & gt = transform->transform_;
    const unsigned int     mesh_size = gt.grid_.mesh_.size();

    // shortcuts to intermediate mesh refinement results for this tile:
    intermediate_result_t & results = results_[tile_index_];
    const unsigned int      num_neighbors = results.dx_.size();

    std::vector<vec2d_t> shift(mesh_size, vec2d(0, 0));
    std::vector<double>  mass(mesh_size, 0);

    for (unsigned int neighbor = 0; neighbor < num_neighbors; neighbor++)
    {
      WRAP(terminator.terminate_on_request());

      // shortcuts:
      image_t::Pointer & dx = results.dx_[neighbor];
      image_t::Pointer & dy = results.dy_[neighbor];
      image_t::Pointer & db = results.db_[neighbor];

      std::vector<vec2d_t> neighbor_pull(mesh_size, vec2d(0, 0));
      regularize_displacements(neighbor_pull, mass, dx, dy, db, median_filter_radius_);

      // blend the displacement vectors:
      for (unsigned int i = 0; i < mesh_size; i++)
      {
        shift[i] += neighbor_pull[i];
      }
    }

    // FIXME: if (num_neighbors > 1)
    if (!keep_first_tile_fixed_)
    {
      // normalize:
      for (unsigned int i = 0; i < mesh_size; i++)
      {
        double scale = 1 / (1 + mass[i]);
        shift[i] *= scale;
      }
    }

    gt.grid_.update(&(shift[0]));
    transform->setup(gt);
  }

  the_log_t &                                log_;
  const unsigned int                         tile_index_;
  const bool                                 keep_first_tile_fixed_;
  const unsigned int                         median_filter_radius_;
  std::vector<itk::GridTransform::Pointer> & transform_;
  std::vector<intermediate_result_t> &       results_;
};


//----------------------------------------------------------------
// refine_mosaic
//
// use FFT to refine the grid transforms directly:
//
template <typename image_t, typename mask_t>
void
refine_mosaic(the_log_t &                                         log,
              std::vector<itk::GridTransform::Pointer> &          transform,
              const std::vector<typename image_t::ConstPointer> & tile,
              const std::vector<typename mask_t::ConstPointer> &  mask,
              const unsigned int &                                neighborhood,
              const bool &                                        prewarp_tiles = true,
              const double &                                      minimum_overlap = 0.25,
              const unsigned int &                                median_radius = 1,
              const unsigned int &                                num_passes = 1,
              const bool &                                        keep_first_tile_fixed = false)
{
  typedef itk::GridTransform itk_grid_transform_t;

  // shortcut:
  unsigned int num_tiles = tile.size();
  if (num_tiles < 2)
    return;

  std::vector<typename image_t::ConstPointer> ConstTile(tile.size());
  ConstTile.assign(tile.begin(), tile.end());

  // image space bounding boxes:
  std::vector<pnt2d_t> image_min;
  std::vector<pnt2d_t> image_max;
  calc_image_bboxes<image_t>(ConstTile, image_min, image_max);

  std::vector<itk_grid_transform_t::ConstPointer> ConstTransform(transform.size());
  ConstTransform.assign(transform.begin(), transform.end());

  // mosaic space bounding boxes:
  std::vector<pnt2d_t> mosaic_min;
  std::vector<pnt2d_t> mosaic_max;
  calc_mosaic_bboxes<pnt2d_t, itk_grid_transform_t>(ConstTransform, image_min, image_max, mosaic_min, mosaic_max, 16);

  unsigned int start = keep_first_tile_fixed ? 1 : 0;

  // find overlapping neighbors for each tile:
  std::vector<the_dynamic_array_t<unsigned int>> neighbors(num_tiles);
  for (unsigned int i = start; i < num_tiles; i++)
  {
    the_aa_bbox_t ibox;
    ibox << p3x1_t(mosaic_min[i][0], mosaic_min[i][1], 0) << p3x1_t(mosaic_max[i][0], mosaic_max[i][1], 0);

    for (unsigned int j = 0; j < num_tiles; j++)
    {
      if (i == j)
        continue;

      the_aa_bbox_t jbox;
      jbox << p3x1_t(mosaic_min[j][0], mosaic_min[j][1], 0) << p3x1_t(mosaic_max[j][0], mosaic_max[j][1], 0);

      if (!ibox.intersects(jbox))
        continue;

      neighbors[i].push_back(j);
    }
  }

  typedef itk::ImageDuplicator<image_t> image_duplicator_t;
  typedef itk::ImageDuplicator<mask_t>  mask_duplicator_t;

  std::vector<typename image_t::Pointer> warped_tile(num_tiles);
  std::vector<typename mask_t::Pointer>  warped_mask(num_tiles);

  if (keep_first_tile_fixed)
  {
    typename image_duplicator_t::Pointer image_duplicator = image_duplicator_t::New();
    typename mask_duplicator_t::Pointer  mask_duplicator = mask_duplicator_t::New();

    image_duplicator->SetInputImage(tile[0]);
    mask_duplicator->SetInputImage(mask[0]);

    image_duplicator->Update();
    mask_duplicator->Update();

    warped_tile[0] = image_duplicator->GetOutput();
    warped_mask[0] = mask_duplicator->GetOutput();
  }

  for (unsigned int pass = 0; pass < num_passes; pass++)
  {
    log << "--------------------------- pass " << pass << " ---------------------------" << endl;

    if (prewarp_tiles)
    {
      // warp the tiles:
      for (unsigned int i = start; i < num_tiles; i++)
      {
        log << setw(4) << i << ". warping image tile" << endl;
        warped_tile[i] = warp<image_t>((typename image_t::ConstPointer)tile[i], transform[i].GetPointer());

        if (mask[i].GetPointer() != nullptr)
        {
          log << "      warping image tile mask" << endl;
          warped_mask[i] = warp<mask_t>((typename mask_t::ConstPointer)mask[i], transform[i].GetPointer());
        }
      }
    }

    std::vector<std::vector<vec2d_t>> shift(num_tiles);
    for (unsigned int i = start; i < num_tiles; i++)
    {
      the_grid_transform_t & gt = transform[i]->transform_;
      unsigned int           mesh_size = gt.grid_.mesh_.size();

      std::vector<std::vector<vec2d_t>> shift_i(neighbors[i].size());
      std::vector<double>               mass(mesh_size, 0);

      for (unsigned int k = 0; k < neighbors[i].size(); k++)
      {
        unsigned int j = neighbors[i][k];
        log << "matching " << i << ":" << j << endl;

        calc_displacements<image_t, mask_t>(*null_log(),
                                            shift_i[k],
                                            mass,

                                            // tiles-already-warped flag:
                                            prewarp_tiles,

                                            // fixed:
                                            warped_tile[j],
                                            warped_mask[j],
                                            transform[j].GetPointer(),

                                            // moving:
                                            warped_tile[i],
                                            warped_mask[i],
                                            transform[i].GetPointer(),

                                            neighborhood,
                                            minimum_overlap,
                                            median_radius);
      }

      // blend the displacement vectors:
      shift[i].assign(mesh_size, vec2d(0, 0));

      for (unsigned int j = 0; j < shift_i.size(); j++)
      {
        for (unsigned int k = 0; k < mesh_size; k++)
        {
          shift[i][k] += shift_i[j][k];
        }
      }

      if (!keep_first_tile_fixed)
      {
        // normalize:
        for (unsigned int k = 0; k < mesh_size; k++)
        {
          double scale = 1 / (1 + mass[k]);
          shift[i][k] *= scale;
        }
      }
    }

    for (unsigned int i = start; i < num_tiles; i++)
    {
      the_grid_transform_t & gt = transform[i]->transform_;
      gt.grid_.update(&(shift[i][0]));
      transform[i]->setup(gt);
    }
  }
}


//----------------------------------------------------------------
// warp_tile_transaction_t
//
template <typename image_t, typename mask_t>
class warp_tile_transaction_t : public the_transaction_t
{
public:
  warp_tile_transaction_t(the_log_t &                              log,
                          unsigned int                             tile_index,
                          itk::GridTransform::Pointer &            transform,
                          const typename image_t::ConstPointer &   tile,
                          const typename mask_t::ConstPointer &    mask,
                          std::vector<typename image_t::Pointer> & warped_tile,
                          std::vector<typename mask_t::Pointer> &  warped_mask)
    : log_(log)
    , tile_index_(tile_index)
    , transform_(transform)
    , tile_(tile)
    , mask_(mask)
    , warped_tile_(warped_tile)
    , warped_mask_(warped_mask)
  {}

  // virtual:
  void
  execute(the_thread_interface_t * thread)
  {
    WRAP(the_terminator_t terminator("warp_tile_transaction_t"));

    log_ << setw(4) << tile_index_ << ". warping image tile" << endl;
    warped_tile_[tile_index_] = warp<image_t>(tile_, transform_.GetPointer());

    if (mask_.GetPointer() != nullptr)
    {
      log_ << setw(4) << tile_index_ << ". warping image tile mask" << endl;
      warped_mask_[tile_index_] = warp<mask_t>(mask_, transform_.GetPointer());
    }
  }

  the_log_t &                              log_;
  unsigned int                             tile_index_;
  itk::GridTransform::Pointer &            transform_;
  const typename image_t::ConstPointer &   tile_;
  const typename mask_t::ConstPointer &    mask_;
  std::vector<typename image_t::Pointer> & warped_tile_;
  std::vector<typename mask_t::Pointer> &  warped_mask_;
};


//----------------------------------------------------------------
// refine_one_tile_t
//
template <typename image_t, typename mask_t>
class refine_one_tile_t : public the_transaction_t
{
public:
  refine_one_tile_t(the_log_t &                                            log,
                    unsigned int                                           tile_index,
                    std::vector<itk::GridTransform::Pointer> &             transform,
                    const std::vector<typename image_t::Pointer> &         warped_tile,
                    const std::vector<typename mask_t::Pointer> &          warped_mask,
                    const std::vector<the_dynamic_array_t<unsigned int>> & neighbors,
                    const bool &                                           tiles_already_warped,
                    const unsigned int &                                   neighborhood_size,
                    const double &                                         minimum_overlap,      // neighbrhood overlap
                    const unsigned int &                                   median_filter_radius, // for outliers
                    const bool &                                           keep_first_tile_fixed,
                    std::vector<std::vector<vec2d_t>> &                    shift)
    :

    log_(log)
    , tile_index_(tile_index)
    , transform_(transform)
    , warped_tile_(warped_tile)
    , warped_mask_(warped_mask)
    , neighbors_(neighbors)
    , tiles_already_warped_(tiles_already_warped)
    , neighborhood_(neighborhood_size)
    , minimum_overlap_(minimum_overlap)
    , median_radius_(median_filter_radius)
    , keep_first_tile_fixed_(keep_first_tile_fixed)
    , shift_(shift)
  {}

  // virtual:
  void
  execute(the_thread_interface_t * thread)
  {
    WRAP(the_terminator_t terminator("refine_one_tile_t"));

    the_grid_transform_t & gt = transform_[tile_index_]->transform_;
    unsigned int           mesh_size = gt.grid_.mesh_.size();
    std::vector<double>    mass(mesh_size, 0);

    unsigned int                      num_neighbors = neighbors_[tile_index_].size();
    std::vector<std::vector<vec2d_t>> shift_i(num_neighbors);

    for (unsigned int k = 0; k < num_neighbors; k++)
    {
      unsigned int j = neighbors_[tile_index_][k];
      log_ << "matching " << tile_index_ << ":" << j << endl;

      calc_displacements<image_t, mask_t>(*null_log(),
                                          shift_i[k],
                                          mass,

                                          // tiles-already-warped flag:
                                          tiles_already_warped_,

                                          // fixed:
                                          warped_tile_[j],
                                          warped_mask_[j],
                                          transform_[j].GetPointer(),

                                          // moving:
                                          warped_tile_[tile_index_],
                                          warped_mask_[tile_index_],
                                          transform_[tile_index_].GetPointer(),

                                          neighborhood_,
                                          minimum_overlap_,
                                          median_radius_);
    }

    // blend the displacement vectors:
    shift_[tile_index_].assign(mesh_size, vec2d(0, 0));

    for (unsigned int j = 0; j < num_neighbors; j++)
    {
      for (unsigned int k = 0; k < mesh_size; k++)
      {
        shift_[tile_index_][k] += shift_i[j][k];
      }
    }

    if (!keep_first_tile_fixed_)
    {
      // normalize:
      for (unsigned int k = 0; k < mesh_size; k++)
      {
        double scale = 1 / (1 + mass[k]);
        shift_[tile_index_][k] *= scale;
      }
    }
  }

  the_log_t &                                            log_;
  unsigned int                                           tile_index_;
  std::vector<itk::GridTransform::Pointer> &             transform_;
  const std::vector<typename image_t::Pointer> &         warped_tile_;
  const std::vector<typename mask_t::Pointer> &          warped_mask_;
  const std::vector<the_dynamic_array_t<unsigned int>> & neighbors_;
  const bool                                             tiles_already_warped_;
  const unsigned int                                     neighborhood_;
  const double                                           minimum_overlap_; // neighbrhood overlap
  const unsigned int                                     median_radius_;   // for outliers
  const bool                                             keep_first_tile_fixed_;
  std::vector<std::vector<vec2d_t>> &                    shift_;
};


//----------------------------------------------------------------
// refine_mosaic_mt
//
template <typename image_t, typename mask_t>
void
refine_mosaic_mt(the_log_t &                                         log, // text output stream
                 std::vector<itk::GridTransform::Pointer> &          transform,
                 const std::vector<typename image_t::ConstPointer> & tile,
                 const std::vector<typename mask_t::ConstPointer> &  mask,
                 const unsigned int &                                neighborhood_size,
                 const bool &                                        prewarp_tiles,   // FIXME: always true?
                 const double &                                      minimum_overlap, // neighbrhood overlap
                 const unsigned int &                                median_radius,   // for outliers
                 const unsigned int &                                num_passes,
                 const bool &                                        keep_first_tile_fixed, // FIXME: stos only?
                 const double &                                      displacement_threshold,
                 unsigned int                                        num_threads) // max concurrent threads
{
  if (num_threads == 1)
  {
    refine_mosaic<image_t, mask_t>(log,
                                   transform,
                                   tile,
                                   mask,
                                   neighborhood_size,
                                   prewarp_tiles,
                                   minimum_overlap,
                                   median_radius,
                                   num_passes,
                                   keep_first_tile_fixed);
    return;
  }

  typedef itk::GridTransform itk_grid_transform_t;

  // shortcut:
  unsigned int num_tiles = tile.size();
  if (num_tiles < 2)
    return;

  log << "num tiles:  " << num_tiles;
  // image space bounding boxes:
  std::vector<typename image_t::ConstPointer> ConstTile;
  std::vector<pnt2d_t>                        image_min;
  std::vector<pnt2d_t>                        image_max;

  ConstTile.reserve(tile.size());
  ConstTile.assign(tile.begin(), tile.end());

  std::vector<typename itk_grid_transform_t::ConstPointer> ConstTransforms;
  ConstTransforms.reserve(transform.size());

  ConstTransforms.assign(transform.begin(), transform.end());

  calc_image_bboxes<image_t>(ConstTile, image_min, image_max);

  // mosaic space bounding boxes:
  std::vector<pnt2d_t> mosaic_min;
  std::vector<pnt2d_t> mosaic_max;
  calc_mosaic_bboxes<pnt2d_t, itk_grid_transform_t>(ConstTransforms, image_min, image_max, mosaic_min, mosaic_max, 16);

  // Relative to the image size.
  // double threshold =
  //  displacement_threshold * std::abs(mosaic_min[0][0] - mosaic_max[0][0]);

  // Relative to a single pixel.
  double threshold = displacement_threshold;

  unsigned int start = keep_first_tile_fixed ? 1 : 0;

  // find overlapping neighbors for each tile:
  std::vector<the_dynamic_array_t<unsigned int>> neighbors(num_tiles);
  for (unsigned int i = start; i < num_tiles; i++)
  {
    the_aa_bbox_t ibox;
    ibox << p3x1_t(mosaic_min[i][0], mosaic_min[i][1], 0) << p3x1_t(mosaic_max[i][0], mosaic_max[i][1], 0);

    for (unsigned int j = 0; j < num_tiles; j++)
    {
      if (i == j)
        continue;

      the_aa_bbox_t jbox;
      jbox << p3x1_t(mosaic_min[j][0], mosaic_min[j][1], 0) << p3x1_t(mosaic_max[j][0], mosaic_max[j][1], 0);

      if (!ibox.intersects(jbox))
        continue;

      neighbors[i].push_back(j);
    }
  }

  std::vector<typename image_t::Pointer> warped_tile(num_tiles);
  std::vector<typename mask_t::Pointer>  warped_mask(num_tiles);

  double last_average = std::numeric_limits<double>::max();

  typedef itk::ImageDuplicator<image_t> image_duplicator_t;
  typedef itk::ImageDuplicator<mask_t>  mask_duplicator_t;

  // Initialize "warped" tiles.
  // If warping is actually requested do it on each pass
  for (int i = 0; i < num_tiles; i++)
  {
    typename image_duplicator_t::Pointer image_duplicator = image_duplicator_t::New();
    typename mask_duplicator_t::Pointer  mask_duplicator = mask_duplicator_t::New();

    image_duplicator->SetInputImage(tile[i]);
    image_duplicator->Update();
    warped_tile[i] = image_duplicator->GetOutput();

    if (mask[i].GetPointer() != nullptr)
    {
      mask_duplicator->SetInputImage(mask[i]);
      mask_duplicator->Update();
      warped_mask[i] = mask_duplicator->GetOutput();
    }
  }

  the_thread_pool_t thread_pool(num_threads);
  thread_pool.set_idle_sleep_duration(50); // 50 usec

  for (unsigned int pass = 0; pass < num_passes; pass++)
  {
    double major_percent = 0.15 + 0.8 * ((double)pass / (double)num_passes);
    double next_major = 0.15 + 0.8 * ((double)(pass + 1) / (double)num_passes);
    set_major_progress(major_percent);

    log << "--------------------------- pass " << pass << " ---------------------------" << endl;

    if (prewarp_tiles)
    {
      // warp the tiles -- split into a set of transactions and executed:

      std::list<the_transaction_t *> schedule;
      for (unsigned int i = start; i < num_tiles; i++)
      {
        warp_tile_transaction_t<image_t, mask_t> * t = new warp_tile_transaction_t<image_t, mask_t>(
          log, i, transform[i], tile[i], mask[i], warped_tile, warped_mask);
        schedule.push_back(t);
      }

      thread_pool.push_back(schedule);
      thread_pool.pre_distribute_work();
      suspend_itk_multithreading_t suspend_itk_mt;
      thread_pool.start();
      thread_pool.wait();
    }

    set_minor_progress(0.2, next_major);

#if 1
    // calculating displacements:
    std::vector<std::vector<vec2d_t>> shift(num_tiles);

    // this is the original coarse-scale parallelization:
    unsigned int num_tiles_distributed = (num_tiles - start) - (num_tiles - start) % num_threads;

    unsigned int num_tiles_remaining = (num_tiles - start) - num_tiles_distributed;

    std::list<the_transaction_t *> schedule;
    for (unsigned int i = 0; i < num_tiles_distributed; i++)
    {
      unsigned int index = start + i;

      typedef refine_one_tile_t<image_t, mask_t> tile_transaction_t;

      tile_transaction_t * t = new tile_transaction_t(log,
                                                      index,
                                                      transform,
                                                      warped_tile,
                                                      warped_mask,
                                                      neighbors,
                                                      prewarp_tiles,
                                                      neighborhood_size,
                                                      minimum_overlap,
                                                      median_radius,
                                                      keep_first_tile_fixed,
                                                      shift);
      schedule.push_back(t);
    }

    thread_pool.push_back(schedule);
    // thread_pool.pre_distribute_work();
    suspend_itk_multithreading_t suspend_itk_mt;
    thread_pool.start();
    thread_pool.wait();

    set_minor_progress(0.9, next_major);

    // this is the broken fine-scale parallelization:
    for (unsigned int i = 0; i < num_tiles_remaining; i++)
    {
      unsigned int index = start + num_tiles_distributed + i;

      the_grid_transform_t & gt = transform[index]->transform_;
      unsigned int           mesh_size = gt.grid_.mesh_.size();

      std::vector<std::vector<vec2d_t>> shift_i(neighbors[index].size());
      std::vector<double>               mass(mesh_size, 0);

      for (unsigned int k = 0; k < neighbors[index].size(); k++)
      {
        unsigned int j = neighbors[index][k];
        log << "matching " << index << ":" << j << endl;

        calc_displacements<image_t, mask_t>(*null_log(),
                                            shift_i[k],
                                            mass,

                                            // tiles-already-warped flag:
                                            prewarp_tiles,

                                            // fixed:
                                            (typename image_t::ConstPointer)warped_tile[j],
                                            (typename mask_t::ConstPointer)warped_mask[j],
                                            (base_transform_t::ConstPointer)transform[j].GetPointer(),

                                            // moving:
                                            (typename image_t::ConstPointer)warped_tile[index],
                                            (typename mask_t::ConstPointer)warped_mask[index],
                                            (itk::GridTransform::ConstPointer)transform[index].GetPointer(),

                                            neighborhood_size,
                                            minimum_overlap,
                                            median_radius);
      }

      // blend the displacement vectors:
      shift[index].assign(mesh_size, vec2d(0, 0));

      for (unsigned int j = 0; j < shift_i.size(); j++)
      {
        for (unsigned int k = 0; k < mesh_size; k++)
        {
          shift[index][k] += shift_i[j][k];
        }
      }

      if (!keep_first_tile_fixed)
      {
        // normalize:
        for (unsigned int k = 0; k < mesh_size; k++)
        {
          double scale = 1 / (1 + mass[k]);
          shift[index][k] *= scale;
        }
      }
    }

    // once all displacement calculations are
    // finished the transform grids can be updated and we
    // can move on to the next pass:

    for (unsigned int i = start; i < num_tiles; i++)
    {
      the_grid_transform_t & gt = transform[i]->transform_;
      gt.grid_.update(&(shift[i][0]));
      transform[i]->setup(gt);
    }

#else
    // this is the "improved" fine scale parallelization:

    // setup intermediate mesh refinement result structures:
    std::vector<intermediate_result_t> results(num_tiles);
    for (unsigned int i = start; i < num_tiles; i++)
    {
      const unsigned int           num_neighbors = neighbors[i].size();
      const the_grid_transform_t & gt = transform[i]->transform_;
      const unsigned int           mesh_rows = gt.rows_ + 1;
      const unsigned int           mesh_cols = gt.cols_ + 1;

      intermediate_result_t & result = results[i];
      result = intermediate_result_t(num_neighbors, mesh_rows, mesh_cols);
    }

    std::list<the_transaction_t *> schedule;
    for (unsigned int i = 0; i < num_threads; i++)
    {
      typedef calc_intermediate_results_t<typename TImage::Pointer, typename TMask::Pointer> tile_transaction_t;

      tile_transaction_t * t = new tile_transaction_t(log,
                                                      i,
                                                      num_threads,
                                                      transform,
                                                      warped_tile,
                                                      warped_mask,
                                                      neighbors,
                                                      prewarp_tiles,
                                                      neighborhood_size,
                                                      minimum_overlap,
                                                      keep_first_tile_fixed,
                                                      results);
      schedule.push_back(t);
    }

    thread_pool.push_back(schedule);
    thread_pool.pre_distribute_work();
    suspend_itk_multithreading_t suspend_itk_mt;
    thread_pool.start();
    thread_pool.wait();

    for (unsigned int i = start; i < num_tiles; i++)
    {
      update_tile_mesh_t * t = new update_tile_mesh_t(log, i, keep_first_tile_fixed, median_radius, transform, results);
      schedule.push_back(t);
    }

    thread_pool.push_back(schedule);
    thread_pool.pre_distribute_work();
    thread_pool.start();
    thread_pool.wait();
#endif
    double worst = 0.0, avg = 0.0, count = 0.0;
    for (int i = 0; i < (int)shift.size(); i++)
    {
      for (int k = 0; k < (int)shift[i].size(); k++)
      {
        if (std::abs(shift[i][k][0]) > worst)
          worst = std::abs(shift[i][k][0]);
        if (std::abs(shift[i][k][1]) > worst)
          worst = std::abs(shift[i][k][1]);

        avg += std::abs(shift[i][k][0]) + std::abs(shift[i][k][1]);
        count += 2;
      }
    }
    avg /= count;
    cout << pass << "  Average Displacement: " << avg << "   Max Displacement: " << worst << endl;

    // If there's an exact cutoff...
    if (count > 0)
    {
      if (avg <= threshold)
        break;
      else if (avg >= last_average)
        break;
      last_average = avg;
    }
  }
}


#endif // MOSAIC_REFINEMENT_COMMON_HXX_
