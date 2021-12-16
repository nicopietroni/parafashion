#ifndef IGL_OUTER_HULL_H
#define IGL_OUTER_HULL_H
#include "../igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Compute the "outer hull" of a potentially non-manifold mesh (V,F) whose
  // intersections have been "resolved" (e.g. using `cork` or
  // `igl::selfintersect`). The outer hull is defined to be all facets
  // (regardless of orientation) for which there exists some path from infinity
  // to the face without intersecting any other facets. For solids, this is the
  // surface of the solid. In general this includes any thin "wings" or "flaps".
  // This implementation largely follows Section 3.6 of "Direct repair of
  // self-intersecting meshes" [Attene 2014].
  //
  // Inputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices into V
  //   N  #F by 3 list of per-face normals
  // Outputs:
  //   G  #G by 3 list of output triangle indices into V
  //   J  #G list of indices into F
  //   flip  #F list of whether facet was added to G **and** flipped orientation
  //     (false for faces not added to G)
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedN,
    typename DerivedG,
    typename DerivedJ,
    typename Derivedflip>
  IGL_INLINE void outer_hull(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    const Eigen::PlainObjectBase<DerivedN> & N,
    Eigen::PlainObjectBase<DerivedG> & G,
    Eigen::PlainObjectBase<DerivedJ> & J,
    Eigen::PlainObjectBase<Derivedflip> & flip);
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedG,
    typename DerivedJ,
    typename Derivedflip>
  IGL_INLINE void outer_hull(
    const Eigen::PlainObjectBase<DerivedV> & V,
    const Eigen::PlainObjectBase<DerivedF> & F,
    Eigen::PlainObjectBase<DerivedG> & G,
    Eigen::PlainObjectBase<DerivedJ> & J,
    Eigen::PlainObjectBase<Derivedflip> & flip);
}

#ifndef IGL_STATIC_LIBRARY
#  include "outer_hull.cpp"
#endif
#endif
