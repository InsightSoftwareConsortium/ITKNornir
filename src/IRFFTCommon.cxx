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


// File         : fft_common.cxx
// Author       : Pavel A. Koshevoy
// Created      : 2005/11/10 14:05
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : Helper functions for image alignment (registration)
//                using phase correlation to find the translation vector.

// local includes:
#include "IRFFTCommon.h"
#include "itkIRUtils.h"

// system includes:
#include <functional>

//----------------------------------------------------------------
// DEBUG_COUNTER1
//
#ifdef DEBUG_PDF
unsigned int DEBUG_COUNTER1 = 0;
unsigned int DEBUG_COUNTER2 = 0;
#endif


//----------------------------------------------------------------
// find_maxima_cm
//
// The percentage below refers to the number of pixels that fall
// below the maxima. Thus, the number of pixels above the maxima
// is 1 - percentage. This way it is possible to specify a
// thresholding value without knowing anything about the image.
//
// The given image is thresholded, the resulting clusters/blobs
// are identified/classified, the center of mass of each cluster
// is treated as a maxima in the image.
//
unsigned int
find_maxima_cm(std::list<local_max_t> &     max_list,
               const itk_image_t::Pointer & image,
               const double                 percentage,
               const the_text_t &           prefix,
               const the_text_t &           suffix)
{
  typedef itk::ImageRegionConstIterator<itk_image_t>          iter_t;
  typedef itk::ImageRegionConstIteratorWithIndex<itk_image_t> itex_t;
  typedef itk_image_t::IndexType                              index_t;
  typedef std::list<index_t>                                  cluster_t;

  // local copy of the image that will be destroyed in the process:
  itk_image_t::Pointer peaks = cast<itk_image_t, itk_image_t>(image);

  // save<itk_image_t>(peaks, "Peaks.png");

  // FIXME:
#ifdef DEBUG_PDF
  if (prefix.size() != 0)
  {
    save<native_image_t>(cast<itk_image_t, native_image_t>(remap_min_max<itk_image_t>(peaks, 0.0, 255.0)),
                         prefix + the_text_t("PDF") + suffix);
  }
#endif

  // first find minimax/maxima of the image:
  double v_min = std::numeric_limits<double>::max();
  double v_max = -v_min;

  iter_t iter(peaks, peaks->GetLargestPossibleRegion());
  for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
  {
    const double v = iter.Get();
    v_min = std::min(v_min, v);
    v_max = std::max(v_max, v);
  }

  // calculate the min/max range:
  const double v_rng = v_max - v_min;

  // NaN is the only number which is not equal to itself
  if (v_rng == 0.0 || v_rng != v_rng || v_rng == std::numeric_limits<double>::infinity())
  {
    // there are no peaks in this image:
    return 0;
  }

  // build a histogram:
  const unsigned int bins = 4096;
  unsigned int       pdf[bins] = { 0 };

  for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
  {
    const double       v = iter.Get();
    const unsigned int bin = (unsigned int)(double((v - v_min) / v_rng) * double(bins - 1));
    pdf[bin]++;
  }

  // build the cumulative histogram:
  unsigned int cdf[bins];
  cdf[0] = pdf[0];
  for (unsigned int i = 1; i < bins; i++)
  {
    cdf[i] = cdf[i - 1] + pdf[i];
  }

  // shortcuts:
  itk_image_t::SizeType size = peaks->GetLargestPossibleRegion().GetSize();
  const unsigned int &  w = size[0];
  const unsigned int &  h = size[1];
  const double          wh = double(w * h);

  // find the CDF bin that contains a given percentage of the total image:
  double clip_min = 0.0;
  for (unsigned int i = 1; i < bins; i++)
  {
    clip_min = v_min + (double(i) / double(bins - 1)) * v_rng;
    if (double(cdf[i]) >= percentage * wh)
      break;
  }

  // threshold the peaks:
  double background = clip_min - v_rng * 1e-3;
  peaks = threshold<itk_image_t>(peaks, clip_min, v_max, background, v_max);
  peaks = remap_min_max<itk_image_t>(peaks, 0.0, 1.0);
  background = 0.0;

  // FIXME:
#ifdef DEBUG_CLUSTERS
  if (prefix.size() != 0)
  {
    save<native_image_t>(cast<itk_image_t, native_image_t>(remap_min_max<itk_image_t>(peaks, 0.0, 255.0)),
                         prefix + the_text_t("clusters") + suffix);
  }
#endif

  // classify the clusters:
  static const int stencil[][2] = { // 4 connected:
                                    { 0, -1 },
                                    { -1, 0 },
                                    { 0, 1 },
                                    { 1, 0 },

                                    // 8 connected:
                                    { -1, -1 },
                                    { 1, 1 },
                                    { -1, 1 },
                                    { 1, -1 }
  };

  the_dynamic_array_t<cluster_t>      clusters;
  the_dynamic_array_t<cluster_bbox_t> bboxes;
  std::vector<unsigned int>           cluster_map(w * h);
  cluster_map.assign(w * h, ~0);

  itex_t itex(peaks, peaks->GetLargestPossibleRegion());
  for (itex.GoToBegin(); !itex.IsAtEnd(); ++itex)
  {
    const double v = itex.Get();

    // skip over the background:
    if (v <= background)
      continue;

    const index_t      index = itex.GetIndex();
    const unsigned int x = index[0];
    const unsigned int y = index[1];

    // iterate over the neighborhood, collect the blob ids
    // of the neighbors:
    // James A: Changed to vector instead of a list... much faster
    std::vector<unsigned int> neighbors;
    neighbors.reserve(8);

    for (unsigned int k = 0; k < 8; k++)
    {
      int u = x + stencil[k][0];
      int v = y + stencil[k][1];
      if ((unsigned int)(u) >= w || (unsigned int)(v) >= h)
        continue;

      unsigned int cluster_id = cluster_map[u * h + v];
      if (cluster_id != (unsigned int)(~0))
      {
        push_back_unique(neighbors, cluster_id);
      }
    }

    if (neighbors.empty())
    {
      // make a new cluster:
      clusters.append(cluster_t());
      bboxes.append(cluster_bbox_t());

      unsigned int id = clusters.end_index(true);
      cluster_map[x * h + y] = id;
      clusters[id].push_back(index);
      bboxes[id].update(x, y);
    }
    else
    {
      // add this pixel to the cluster:
      unsigned int id = *(neighbors.begin());
      cluster_map[x * h + y] = id;
      clusters[id].push_back(index);
      bboxes[id].update(x, y);

      if (neighbors.size() > 1)
      {
        // merge the clusters into one (the first one):
        // James A: Changed to vector instead of a list... much faster
        std::vector<unsigned int>::iterator bi = ++(neighbors.begin());
        for (; bi != neighbors.end(); ++bi)
        {
          unsigned int old_id = *bi;
          bboxes[old_id].reset();

          while (!clusters[old_id].empty())
          {
            index_t ij = remove_head(clusters[old_id]);

            cluster_map[ij[0] * h + ij[1]] = id;
            clusters[id].push_back(ij);
            bboxes[id].update(ij[0], ij[1]);
          }
        }
      }
    }
  }

  // merge the clusters that are broken up across the periodic boundary:
  for (unsigned int i = 0; i < clusters.size(); i++)
  {
    cluster_t & cluster = clusters[i];
    if (cluster.empty())
      continue;

    for (std::list<index_t>::iterator j = cluster.begin(); j != cluster.end(); ++j)
    {
      const index_t index = *j;
      unsigned int  x = (index[0] + w) % w;
      unsigned int  y = (index[1] + h) % h;

      for (unsigned int k = 0; k < 8; k++)
      {
        int dx = stencil[k][0];
        int dy = stencil[k][1];

        // adjust for periodicity:
        int u = (x + dx + w) % w;
        int v = (y + dy + h) % h;

        unsigned int cluster_id = cluster_map[u * h + v];
        if (cluster_id == i || cluster_id == (unsigned int)(~0))
          continue;

        // figure out which boundaries this cluster was broken accross:
        cluster_bbox_t & ba = bboxes[i];
        cluster_bbox_t & bb = bboxes[cluster_id];

        bool merge_x = ((bb.max_[0] - ba.min_[0] > int(w / 2)) || (ba.max_[0] - bb.min_[0] > int(w / 2)));

        bool merge_y = ((bb.max_[1] - ba.min_[1] > int(h / 2)) || (ba.max_[1] - bb.min_[1] > int(h / 2)));

        int shift_x = (!merge_x) ? 0 : (ba.min_[0] <= 0) ? -w : w;
        int shift_y = (!merge_y) ? 0 : (ba.min_[1] <= 0) ? -h : h;

        // merge the clusters into one (the first one):
        cluster_t & neighbor = clusters[cluster_id];

        bb.reset();
        while (!neighbor.empty())
        {
          index_t ij = remove_head(neighbor);
          cluster_map[ij[0] * h + ij[1]] = i;

          ij[0] += shift_x;
          ij[1] += shift_y;
          clusters[i].push_back(ij);
          ba.update(ij[0], ij[1]);
        }
      }
    }
  }

  // FIXME:
#ifdef DEBUG_MARKERS
  itk_image_t::Pointer markers = make_image<itk_image_t>(size, background);
#endif

  // calculate the center of mass for each cluster:
  unsigned int num_peaks = 0;
  for (unsigned int i = 0; i < clusters.size(); i++)
  {
    const cluster_t & cluster = clusters[i];
    if (cluster.empty())
      continue;

    double mx = 0.0;
    double my = 0.0;
    double mt = 0.0;

    for (std::list<index_t>::const_iterator j = cluster.begin(); j != cluster.end(); ++j)
    {
      index_t ij = *j;
      double  x = double(ij[0]);
      double  y = double(ij[1]);

      // adjust index for periodicity:
      if (x < 0)
        ij[0] += w;
      if (x >= w)
        ij[0] -= w;
      if (y < 0)
        ij[1] += h;
      if (y >= h)
        ij[1] -= h;

      double m = peaks->GetPixel(ij);
      mx += m * x;
      my += m * y;
      mt += m;
    }

    double cm_x = mx / mt;
    double cm_y = my / mt;
    double m = mt / double(cluster.size());

    // FIXME:
#ifdef DEBUG_MARKERS
    mark<itk_image_t>(markers, pnt2d(cm_x, cm_y), m, 2, '+');
#endif

    max_list.push_back(local_max_t(m, cm_x, cm_y, cluster.size()));
    num_peaks++;
  }

  // FIXME:
#ifdef DEBUG_MARKERS
  save<native_image_t>(cast<itk_image_t, native_image_t>(remap_min_max<itk_image_t>(markers, 0.0, 255.0)),
                       prefix + the_text_t("markings") + suffix,
                       false);
#endif

  // sort the max points so that the best candidate is first:
  max_list.sort(std::greater<local_max_t>());

  return num_peaks;
}


//----------------------------------------------------------------
// threshold_maxima
//
// Discard maxima whose mass is below a given threshold ratio
// of the total mass of all maxima:
//
void
threshold_maxima(std::list<local_max_t> & max_list, const double threshold)
{
  double total_mass = 0.0;
  for (std::list<local_max_t>::iterator i = max_list.begin(); i != max_list.end(); ++i)
  {
    double mass = double((*i).area_) * ((*i).value_);
    total_mass += mass;
  }

  std::list<local_max_t> new_list;
  double                 threshold_mass = threshold * total_mass;
  for (std::list<local_max_t>::iterator i = max_list.begin(); i != max_list.end(); ++i)
  {
    double mass = double((*i).area_) * (*i).value_;
    if (mass < threshold_mass)
      continue;

    new_list.push_back(*i);
  }

  max_list = new_list;
}


//----------------------------------------------------------------
// reject_negligible_maxima
//
// Discard maxima that are worse than the best maxima by a factor
// greater than the given threshold ratio:
//
unsigned int
reject_negligible_maxima(std::list<local_max_t> & max_list, const double threshold)
{
  double best_mass = 0.0;
  for (std::list<local_max_t>::iterator i = max_list.begin(); i != max_list.end(); ++i)
  {
    double mass = (*i).value_;
    best_mass = std::max(best_mass, mass);
  }

  unsigned int           new_size = 0;
  std::list<local_max_t> new_list;
  while (!max_list.empty())
  {
    local_max_t lm = remove_head(max_list);
    double      mass = lm.value_;
    if (best_mass / mass > threshold)
      continue;

    new_list.push_back(lm);
    new_size++;
  }

  max_list.splice(max_list.end(), new_list);
  return new_size;
}


//----------------------------------------------------------------
// reject_negligible_overlap
//
void
reject_negligible_overlap(std::list<overlap_t> & ol, const double threshold)
{
  double best_overlap = 0.0;
  for (std::list<overlap_t>::iterator i = ol.begin(); i != ol.end(); ++i)
  {
    double overlap = (*i).overlap_;
    best_overlap = std::max(best_overlap, overlap);
  }

  std::list<overlap_t> new_list;
  for (std::list<overlap_t>::iterator i = ol.begin(); i != ol.end(); ++i)
  {
    double overlap = (*i).overlap_;
    if (overlap == 0.0)
      continue;
    if (best_overlap / overlap > threshold)
      continue;

    new_list.push_back(*i);
  }

  ol = new_list;
}


//----------------------------------------------------------------
// find_correlation
//
template <>
unsigned int
find_correlation(std::list<local_max_t> & max_list,
                 const itk_image_t *      fi,
                 const itk_image_t *      mi,

                 // low pass filter parameters
                 // (resampled data requires less smoothing):
                 double       lp_filter_r,
                 double       lp_filter_s,
                 const double overlap_min,
                 const double overlap_max)
{
  itk_image_t::SizeType max_sz = calc_padding<itk_image_t>(fi, mi);

  typedef itk_image_t::SizeType   sz_t;
  typedef itk_image_t::RegionType rn_t;

  rn_t fi_region = fi->GetLargestPossibleRegion();
  sz_t fi_size = fi_region.GetSize();

  rn_t mi_region = mi->GetLargestPossibleRegion();
  sz_t mi_size = mi_region.GetSize();

  itk_image_t::ConstPointer z0 =
    (fi_size[0] == max_sz[0] && fi_size[1] == max_sz[1]) ? fi : pad<itk_image_t>(fi, max_sz);
  itk_image_t::ConstPointer z1 =
    (mi_size[0] == max_sz[0] && mi_size[1] == max_sz[1]) ? mi : pad<itk_image_t>(mi, max_sz);

  fft_data_t f0;
  fft(z0, f0);
  f0.apply_lp_filter(lp_filter_r, lp_filter_s);

  fft_data_t f1;
  fft(z1, f1);
  f1.apply_lp_filter(lp_filter_r, lp_filter_s);

  const unsigned int & nx = f0.nx();
  const unsigned int & ny = f0.ny();
  fft_data_t           P(nx, ny);

  // Blank out areas outside the min or max overlap
  unsigned int xLeftBorderStart = std::ceil(nx * overlap_min);
  unsigned int xRightBorderStart = nx - xLeftBorderStart;
  unsigned int yLowBorderStart = std::ceil(ny * overlap_min);
  unsigned int yHighBorderStart = ny - yLowBorderStart;


  unsigned int xLeftBorderCutoff = std::ceil(nx * overlap_max);
  unsigned int xRightBorderCutoff = nx - xLeftBorderCutoff;
  unsigned int yLowBorderCutoff = std::ceil(ny * overlap_max);
  unsigned int yHighBorderCutoff = ny - yLowBorderCutoff;

  for (unsigned int x = 0; x < nx; x++)
  {
    for (unsigned int y = 0; y < ny; y++)
    {
      // TODO: Skip pixels we know can't possibly be overlapping
      /*
      if(!(y <= yLowBorderCutoff ||
           y >= yHighBorderCutoff ||
           x <= xLeftBorderCutoff ||
           x >= xRightBorderCutoff))
        {
            P(x, y) = 0;
            continue;
        }
      */
#if 1
      // Girod-Kuo, normalized cross power spectrum,
      // corresponds to phase correlation in spatial domain:
      fft_complex_t p10 = f1(x, y) * std::conj(f0(x, y));
      P(x, y) = _div(p10, _add(std::sqrt(p10 * std::conj(p10)), 1e-8f));
#else
      // cross power spectrum,
      // corresponds to cross correlation in spatial domain:
      P(x, y) = f1(x, y) * std::conj(f0(x, y));
#endif
    }
  }

  // resampled data produces less noisy PDF and requires less smoothing:
  P.apply_lp_filter(lp_filter_r * 0.8, lp_filter_s);

  // calculate the displacement probability density function:
  fft_data_t ifft_P;

#ifndef NDEBUG // get around an annoying compiler warning:
  bool ok =
#endif
    ifft(P, ifft_P);
  assert(ok);

  itk_image_t::Pointer PDF = ifft_P.real();

  itk_image_t::PixelType min;
  itk_image_t::PixelType max;
  image_min_max<itk_image_t>(PDF.GetPointer(), min, max);

  // save<native_image_t>(cast<itk_image_t, native_image_t>(remap_min_max<itk_image_t>(PDF, 0.0, 255.0)),
  // "Prefill.tif");

  // Set areas that cannot be a match to zero
  // TODO: //Need to fill with an ellipse to get mins/maxs correct

  //	if(xLeftBorderCutoff < xRightBorderCutoff &&
  //	   yLowBorderCutoff < yHighBorderCutoff)
  //	{
  //		fill<itk_image_t>(PDF, xLeftBorderCutoff, yLowBorderCutoff, xRightBorderCutoff-xLeftBorderCutoff,
  //yHighBorderCutoff-yLowBorderCutoff, min);
  //	}
  //
  //  #ifdef DEBUG
  //	save<native_image_t>(cast<itk_image_t, native_image_t>(remap_min_max<itk_image_t>(PDF, 0.0, 255.0)),
  //"Prefill.tif");
  //  #endif

  int                    PixelsInOverlapZone = 0;
  itk_image_t::IndexType iPixel;
  for (iPixel[1] = 0; iPixel[1] <= ny / 2; iPixel[1]++)
  {
    vec2d_t pt;
    for (iPixel[0] = 0; iPixel[0] <= nx / 2; iPixel[0]++)
    {
      pt[0] = iPixel[0];
      pt[1] = iPixel[1];

      double overlap = OverlapPercent(fi_size, pt);
      if (overlap >= overlap_min && overlap <= overlap_max)
      {
        PixelsInOverlapZone += 4;
        continue;
      }

      pt[0] = nx - iPixel[0];
      pt[1] = iPixel[1];

      overlap = OverlapPercent(fi_size, pt);
      if (overlap >= overlap_min && overlap <= overlap_max)
      {
        PixelsInOverlapZone += 4;
        continue;
      }

      pt[0] = iPixel[0];
      pt[1] = ny - iPixel[1];

      overlap = OverlapPercent(fi_size, pt);
      if (overlap >= overlap_min && overlap <= overlap_max)
      {
        PixelsInOverlapZone += 4;
        continue;
      }

      pt[0] = nx - iPixel[0];
      pt[1] = ny - iPixel[1];

      overlap = OverlapPercent(fi_size, pt);
      if (overlap >= overlap_min && overlap <= overlap_max)
      {
        PixelsInOverlapZone += 4;
        continue;
      }

      // If we got here the pixel can't be overlapping
      itk_image_t::IndexType iSetPixel = iPixel;
      PDF->SetPixel(iSetPixel, min);

      iSetPixel[0] = (nx - 1) - iPixel[0];
      PDF->SetPixel(iSetPixel, min);

      iSetPixel[0] = iPixel[0];
      iSetPixel[1] = (ny - 1) - iPixel[1];
      PDF->SetPixel(iSetPixel, min);

      iSetPixel[0] = (nx - 1) - iPixel[0];
      iSetPixel[1] = (ny - 1) - iPixel[1];
      PDF->SetPixel(iSetPixel, min);
      /*
      if(overlap < overlap_min)
        PDF->SetPixel(iPixel, min);

      else if(overlap > overlap_max)
        PDF->SetPixel(iPixel, min);
      */
    }
  }

  /*
    itk_image_t::IndexType ix;
    ix[0] = 0;
    ix[1] = 0;
    PDF->SetPixel(ix, min);
    ix[0] = nx-1;
    ix[1] = 0;
    PDF->SetPixel(ix, min);
    ix[0] = 0;
    ix[1] = ny-1;
    PDF->SetPixel(ix, min);
    ix[0] = nx-1;
    ix[1] = ny-1;
    PDF->SetPixel(ix, min);
  */
  //	fill<itk_image_t>(PDF, 0, 0, xLeftBorderStart, ny, min);
  //	ill<itk_image_t>(PDF, xRightBorderStart, 0, nx-xRightBorderStart, ny, min);
  //	fill<itk_image_t>(PDF, xLeftBorderStart, 0, xRightBorderStart-xLeftBorderStart, yLowBorderStart, min);
  //	fill<itk_image_t>(PDF, xLeftBorderStart, yHighBorderStart, xRightBorderStart-xLeftBorderStart,
  //ny-yHighBorderStart, min);

  //  #ifdef DEBUG
  //	save<native_image_t>(cast<itk_image_t, native_image_t>(remap_min_max<itk_image_t>(PDF, 0.0, 255.0)),
  //"Postfill.tif");
  //  #endif

  // look for the maxima in the PDF,
  // TODO: Should we count the dead regions towards our histogram because they'll be examined?
  // double area = double(max_sz[0] * max_sz[1]);
  // if(xLeftBorderCutoff < xRightBorderCutoff &&
  // yLowBorderCutoff < yHighBorderCutoff)
  // area -= ((xRightBorderCutoff-xLeftBorderCutoff) * (yHighBorderCutoff-yLowBorderCutoff));
  double area = PixelsInOverlapZone;

  // a minimum of 5 pixels and a maximum of 64 pixels may be attributed
  // to local maxima in the image:
  // double fraction = std::min(64.0 / area, std::max(5.0 / area, 1e-3));
  double fraction = std::min(64.0 / area, std::max(5.0 / area, 1e-2));

  // the entire image should never be treated as a maxima cluster:
  assert(fraction < 1.0);

  // find the maxima clusters:
  return find_maxima_cm(max_list, PDF, 1.0 - fraction);
}
