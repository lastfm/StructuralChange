#ifndef _STRUCTURAL_CHANGE_HPP
#define _STRUCTURAL_CHANGE_HPP

/* -*- c-basic-offset: 3 indent-tabs-mode: nil -*-  vi:set ts=8 sts=3 sw=3: */

/*
 
   Structural Change
   Authors: Matthias Mauch, Marcus Holland-Moritz at Last.fm, 2011
   
   This is and implementation of the Structural Change meta feature
   described in the article "Structural Change on Multiple Time Scales 
   as a Correlate of Musical Complexity" by Matthias Mauch and Mark Levy,
   published at the 12th International Conference on Music Information
   Retrieval (ISMIR 2011), available at 
 
      http://ismir2011.ismir.net/papers/PS4-3.pdf
 
   Contact: {matthias,mark}@last.fm
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.  See the file
   COPYING included with this distribution for more information.
  
*/

#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>

namespace fm { namespace last { namespace audio {

/*

   The DivergencePolicy classes below implement the respective functionality
   for several divergence functions that can be used to compare the left
   and right windows during the Structural Change calculation.

   The reason that these are classes outside the main StructuralChange class
   is that the Mahalanobis distance needs additional information (the 
   inverse covariance matrix), and this way the same function 
   StructuralChange::calculate can be used for all cases; you can also
   implement your own divergence (distance) measure without touching this code.

*/

class CorrelationDivergencePolicy // implements a divergence measure from the correlation of two vectors
{ 
public:
   float
   operator() (const std::vector<float>& a, const std::vector<float>& b) const
   {
      size_t num_element = a.size();

      float a_mean = 0;
      float b_mean = 0;
      float a_shifted;
      float b_shifted;
      bool a_valid = false;
      bool b_valid = false;

      for (size_t i_element = 0; i_element < num_element; ++i_element) // get the mean
      {
         a_mean += a[i_element];
         b_mean += b[i_element];

         // elements must not be all the same
         if (i_element > 0) 
         {
            if (a[i_element] != a[i_element - 1]) a_valid = true;
            if (b[i_element] != b[i_element - 1]) b_valid = true;
         }
      }

      if (a_valid && b_valid)
      {
         a_mean = a_mean / num_element;
         b_mean = b_mean / num_element;

         float above = 0; // referring to the definition of r = ... on the Wikipedia page

         float below_a = 0;
         float below_b = 0;

         for (size_t i_element = 0; i_element < num_element; ++i_element)  // rest of the loops for stdev calculation
         {
            a_shifted = a[i_element] - a_mean;
            b_shifted = b[i_element] - b_mean;
            above += a_shifted * b_shifted;
            below_a += a_shifted * a_shifted;
            below_b += b_shifted * b_shifted;
         }

         return  0.5 + ( - 0.5 * above / (std::sqrt(below_a) * std::sqrt(below_b)) ); // the 0.5 + ( - 0.5 * ... bit is to squeeze it the right way into [0,1]
      }
      else
      {
         return 0.0;
      }
   }
};


class JensonShannonDivergencePolicy // implements the Jenson Shannon divergence
{
public:
   float
   operator() (std::vector<float> a, std::vector<float> b) const
   {
      size_t num_element = a.size();

      float a_sum = 0;
      float b_sum = 0;
      bool a_valid = false;
      bool b_valid = false;

      for (size_t i_element = 0; i_element < num_element; ++i_element)  // get the sums
      {
         a_sum += a[i_element];
         b_sum += b[i_element];

         if (a[i_element] > 0.0) a_valid = true;
         if (b[i_element] > 0.0) b_valid = true;

         // have to do checks separately because the condition is: each element non-negative AND at least one positive!
         if (a[i_element] < 0.0)
         {
            std::cerr << "ERROR: numbers have to be greater than 0.\n";
            a_valid = false;
            break;
         }

         if (b[i_element] < 0.0)
         {
            std::cerr << "ERROR: numbers have to be greater than 0.\n";
            b_valid = false;
            break;
         }
      }

      if (a_valid && b_valid)
      {
         std::vector<float> m(num_element); // will be the mean of a and b
         float js_divergence = 0;

         for (size_t i_element = 0; i_element < num_element; ++i_element)
         { 
            a[i_element] /= a_sum;
            b[i_element] /= b_sum;

            m[i_element] = 0.5 * (a[i_element] + b[i_element]);

            if (a[i_element] > 0) js_divergence += a[i_element] * log(a[i_element] / m[i_element]);
            if (b[i_element] > 0) js_divergence += b[i_element] * log(b[i_element] / m[i_element]);
         }

         return js_divergence * 0.5;
      } 
      else 
      {
         return 0.0;
      }
   }
};

class MahalanobisDivergencePolicy // implements the Mahalanobis distance
{
public:
   MahalanobisDivergencePolicy(const std::vector<std::vector<float> >& inv_cov)
      : m_inv_cov(inv_cov)
   {
   }

   float
   operator() (std::vector<float> a, std::vector<float> b) const
   {
      size_t num_element = a.size();

      if (num_element > m_inv_cov.size())
      {
         num_element = m_inv_cov.size();
      }

      float outsum = 0;

      for (size_t i = 0; i < num_element; ++i)
      {
         float currval = 0;

         for (size_t j = 0; j < num_element; ++j) 
         {
             currval += m_inv_cov[i][j] * (a[j]-b[j]);
         }

         outsum += currval * (a[i]-b[i]);
      }

      return sqrt(outsum);
   }

private:
   const std::vector<std::vector<float> > m_inv_cov;
};

class EuclideanDivergencePolicy // Euclidean distance
{
public:
   float
   operator() (std::vector<float> a, std::vector<float> b) const
   {
      size_t num_element = a.size();
      
      float outsum = 0;

      for (size_t i = 0; i < num_element; ++i)
      {
         outsum += (a[i] - b[i]) * (a[i] - b[i]);
      }

      return sqrt(outsum);
   }
};

template <typename>
struct FeaturesAccessPolicy;
 
template <>
struct FeaturesAccessPolicy< std::vector<float> >
{
   typedef std::vector<float> features_t;

   static const std::vector<float>& get(const features_t& v)
   {
      return v;
   }

   static std::vector<float>& get(features_t& v)
   {
      return v;
   }

   template <typename T>
   static void copy_meta(features_t&, const T&)
   {
   }
};

template <typename Dest, typename Source>
struct FeaturesCopyMetaPolicy
{
   static void copy_meta(Dest&, const Source&)
   {
   }
};

class StructuralChange 
{
public:
   StructuralChange(size_t num_dim)
      : m_num_dim(num_dim)
   {
   }
   
   // This method takes the main computational load, it modifies the vector rv via reference such
   // that its elements contain Structural Change values.
   template <class OutputFeatures, class InputFeatures, class DivergencePolicy>
   void
   calculate(std::vector<OutputFeatures>& rv, const std::vector<InputFeatures>& input_vector, DivergencePolicy& divergence) const
   {
      size_t num_frame = input_vector.size();

      if (num_frame == 0)
      {
         rv.clear();
         return;
      }

      size_t num_element = FeaturesAccessPolicy<InputFeatures>::get(input_vector[0]).size();

      std::vector<OutputFeatures> change_matrix(num_frame);
      std::vector<float> cumulative_vector (num_element);
      std::vector<std::vector<float> > cumulative_matrix;

      // Accumulating the values over frames once here saves many, many additions later.
      cumulative_matrix.push_back(std::vector<float>(num_element)); // initial row of zeros.

      for (size_t i_frame = 0; i_frame < num_frame; ++i_frame)
      {
         for (size_t i_element = 0; i_element < num_element; ++i_element)
         {
            cumulative_vector[i_element] += FeaturesAccessPolicy<InputFeatures>::get(input_vector[i_frame])[i_element];
         }
         cumulative_matrix.push_back(cumulative_vector);
         FeaturesAccessPolicy<OutputFeatures>::get(change_matrix[i_frame]).resize(m_num_dim);
         FeaturesCopyMetaPolicy<OutputFeatures, InputFeatures>::copy_meta(change_matrix[i_frame], input_vector[i_frame]);
      }

      // get window boundaries for all timescales and frame
      std::vector<std::vector<std::vector<int> > > window_boundaries;
      init_window_boundaries(window_boundaries, num_frame);

      std::vector<float> mean_L (num_element);
      std::vector<float> mean_R (num_element);

      float mean_div; // this is needed to make Sonic Annotator mean and median calculations meaningful
      int count_valid;

      for (size_t i_timescale = 0; i_timescale < m_num_dim; ++i_timescale)  // loop over all different times scales
      { 
         // LEFT and RIGHT are placeholders in places where no valid structural change an be calculated
         // (if either of the windows are too small); they will be replaced by values that make just that
         // allow Sonic Annotator to take mean/median values over *all* frames.
         const float LEFT  = -1.0f;
         const float RIGHT = -2.0f;

         mean_div = 0;
         count_valid = 0;

         // loop over all frames
         for (size_t i_boundary = 0; i_boundary < window_boundaries[i_timescale].size(); ++i_boundary) 
         { 
            if ( window_boundaries[i_timescale][i_boundary][4] == 0 )                     // normal case
            {
               int current_num_frames = window_boundaries[i_timescale][i_boundary][1] - window_boundaries[i_timescale][i_boundary][0];

               // loop over dimensions of the input vector and get local mean in the windows
               for (size_t i_element = 0; i_element < num_element; ++i_element)                              
               { 
                  mean_L[i_element] = (cumulative_matrix[window_boundaries[i_timescale][i_boundary][1]][i_element] -
                                  cumulative_matrix[window_boundaries[i_timescale][i_boundary][0]][i_element]) / current_num_frames;
                  mean_R[i_element] = (cumulative_matrix[window_boundaries[i_timescale][i_boundary][3]][i_element] -
                                  cumulative_matrix[window_boundaries[i_timescale][i_boundary][2]][i_element]) / current_num_frames;
               }

               FeaturesAccessPolicy<OutputFeatures>::get(change_matrix[i_boundary])[i_timescale] = divergence(mean_L, mean_R);
               mean_div += FeaturesAccessPolicy<OutputFeatures>::get(change_matrix[i_boundary])[i_timescale]; 
               ++count_valid;
            } 
            else if ( window_boundaries[i_timescale][i_boundary][4] == -1 )              // we're to the left of the valid part
            {
               FeaturesAccessPolicy<OutputFeatures>::get(change_matrix[i_boundary])[i_timescale] = LEFT;
            } 
            else if ( window_boundaries[i_timescale][i_boundary][4] == -2 )              // on the right of the valid part
            {
               FeaturesAccessPolicy<OutputFeatures>::get(change_matrix[i_boundary])[i_timescale] = RIGHT;
            } 
            else FeaturesAccessPolicy<OutputFeatures>::get(change_matrix[i_boundary])[i_timescale] = 0.0;
         }

         if (count_valid > 0) mean_div /= count_valid;

         // Set artificial values for invalid regions. This results in correct mean and median summaries in Sonic Annotator.
         for (size_t i_frame = 0; i_frame < num_frame; ++i_frame)
         {
            if (FeaturesAccessPolicy<OutputFeatures>::get(change_matrix[i_frame])[i_timescale] == RIGHT)
               FeaturesAccessPolicy<OutputFeatures>::get(change_matrix[i_frame])[i_timescale] = 3 * mean_div;
            else if (FeaturesAccessPolicy<OutputFeatures>::get(change_matrix[i_frame])[i_timescale] == LEFT)
               FeaturesAccessPolicy<OutputFeatures>::get(change_matrix[i_frame])[i_timescale] = (-1) * mean_div;
         }
      }

      rv.swap(change_matrix);
   }

   // Default method.
   template <class OutputFeatures, class InputFeatures>
   void
   calculate(std::vector<OutputFeatures>& rv, const std::vector<InputFeatures>& input_vector) const
   {
      JensonShannonDivergencePolicy dp;
      calculate(rv, input_vector, dp);
   }

private:
   const size_t m_num_dim;

   void
   init_window_boundaries(std::vector<std::vector<std::vector<int> > >& window_boundaries, size_t num_frame) const
   {
      std::vector<std::vector<std::vector<int> > > tmp_window_boundaries;

      tmp_window_boundaries.resize(m_num_dim);
   
      for (size_t i_timescale = 0; i_timescale < m_num_dim; ++i_timescale)
      {
         size_t window_width = 1U << i_timescale;  // smallest window is always one frame
   
         tmp_window_boundaries[i_timescale].resize(num_frame);
   
         for (size_t i_frame = 0; i_frame < num_frame; ++i_frame)
         {
            std::vector<int>& window_boundary = tmp_window_boundaries[i_timescale][i_frame];
   
            // The following is a bit quirky because it gets the cumulative matrix's indices.
            // (The cumulative_matrix has an additional vector of zeros at position 0.)
   
            // left and right boundaries
            int ll = (i_frame + 1 > window_width) ? i_frame - window_width : 0;
            int rr = (i_frame + window_width < num_frame) ? i_frame + window_width : num_frame;
   
            window_boundary.push_back(ll);
            window_boundary.push_back(i_frame);
            window_boundary.push_back(i_frame); // redundant here, could become meaningful for different window type
            window_boundary.push_back(rr);
   
            if (2*static_cast<int>(window_width) == (rr - ll))      // normal case
            {
               window_boundary.push_back(0);
            }
            else if (window_width == (rr - i_frame))                // left window too short
            {
               window_boundary.push_back(-1);
            }
            else if (window_width == (i_frame - ll))                // left window too short
            {
               window_boundary.push_back(-2);
            }
            else                                                    // both windows too short
            {
               window_boundary.push_back(-3);
            }
         }
      }

      window_boundaries.swap(tmp_window_boundaries);
   }
};

}}}

#endif
