#ifndef _VAMP_FEATURES_ACCESS_HPP
#define _VAMP_FEATURES_ACCESS_HPP

/* -*- c-basic-offset: 3 indent-tabs-mode: nil -*-  vi:set ts=8 sts=3 sw=3: */

/*

   Vamp Features Access
   Author: Matthias Mauch, Marcus Holland-Moritz at Last.fm, 2011

   Some templates to allow using the Structural Change implementation with 
   Vamp feature vectors.

*/

#include <vector>
#include <vamp-sdk/Plugin.h>

namespace fm { namespace last { namespace audio {

template <>
struct FeaturesAccessPolicy< _VampPlugin::Vamp::Plugin::Feature >
{
   typedef _VampPlugin::Vamp::Plugin::Feature features_t;

   static const std::vector<float>& get(const features_t& v)
   {
      return v.values;
   }

   static std::vector<float>& get(features_t& v)
   {
      return v.values;
   }
};

template <typename Source>
struct FeaturesCopyMetaPolicy<_VampPlugin::Vamp::Plugin::Feature, Source>
{
   static void copy_meta(_VampPlugin::Vamp::Plugin::Feature& lhs, const Source&)
   {
      lhs.hasTimestamp = false;
   }
};

template <>
struct FeaturesCopyMetaPolicy<_VampPlugin::Vamp::Plugin::Feature, _VampPlugin::Vamp::Plugin::Feature>
{
   static void copy_meta(_VampPlugin::Vamp::Plugin::Feature& lhs, const _VampPlugin::Vamp::Plugin::Feature& rhs)
   {
      lhs.hasTimestamp = rhs.hasTimestamp;
      if (rhs.hasTimestamp)
         lhs.timestamp = rhs.timestamp;
   }
};

}}}

#endif
