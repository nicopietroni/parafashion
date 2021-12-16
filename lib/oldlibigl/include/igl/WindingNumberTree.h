// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_WINDINGNUMBERTREE_H
#define IGL_WINDINGNUMBERTREE_H
#include <list>
#include <map>
#include <Eigen/Dense>
#include "WindingNumberMethod.h"

namespace igl
{
  // This is only need to fill in references, it should never actually be touched
  // and shouldn't cause race conditions. (This is a hack, but I think it's "safe")
  static Eigen::MatrixXd dummyV;
  // Space partitioning tree for computing winding number hierarchically.
  //
  // Templates:
  //   Point  type for points in space, e.g. Eigen::Vector3d
  template <typename Point>
  class WindingNumberTree
  {
    public:
      // Method to use (see enum above)
      //static double min_max_w;
      static std::map< 
        std::pair<const WindingNumberTree*,const WindingNumberTree*>, double>
          cached;
    protected:
      WindingNumberMethod method;
      const WindingNumberTree * parent;
      std::list<WindingNumberTree * > children;
      //// List of boundary edges (recall edges are vertices in 2d)
      //const Eigen::MatrixXi boundary;
      // Base mesh vertices
      Eigen::MatrixXd & V;
      // Base mesh vertices with duplicates removed
      Eigen::MatrixXd SV;
      // Facets in this bounding volume
      Eigen::MatrixXi F;
      // Tesselated boundary curve
      Eigen::MatrixXi cap;
      // Upper Bound on radius of enclosing ball
      double radius;
      // (Approximate) center (of mass)
      Point center;
    public:
      inline WindingNumberTree();
      // For root
      inline WindingNumberTree(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F);
      // For chilluns 
      inline WindingNumberTree(
        const WindingNumberTree<Point> & parent,
        const Eigen::MatrixXi & F);
      inline virtual ~WindingNumberTree();
      inline void delete_children();
      inline virtual void set_mesh(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F);
      // Set method
      inline void set_method( const WindingNumberMethod & m);
    public:
      inline const Eigen::MatrixXd & getV() const;
      inline const Eigen::MatrixXi & getF() const;
      inline const Eigen::MatrixXi & getcap() const;
      // Grow the Tree recursively
      inline virtual void grow();
      // Determine whether a given point is inside the bounding 
      //
      // Inputs:
      //   p  query point 
      // Returns true if the point p is inside this bounding volume
      inline virtual bool inside(const Point & p) const;
      // Compute the (partial) winding number of a given point p
      // According to method
      //  
      // Inputs:
      //   p  query point 
      // Returns winding number 
      inline double winding_number(const Point & p) const;
      // Same as above, but always computes winding number using exact method
      // (sum over every facet)
      inline double winding_number_all(const Point & p) const;
      // Same as above, but always computes using sum over tesslated boundary
      inline double winding_number_boundary(const Point & p) const;
      //// Same as winding_number above, but if max_simple_abs_winding_number is
      //// less than some threshold min_max_w just return 0 (colloquially the "fast
      //// multipole method)
      ////
      ////
      //// Inputs:
      ////   p  query point 
      ////   min_max_w  minimum max simple w to be processed
      //// Returns approximate winding number
      //double winding_number_approx_simple(
      //  const Point & p, 
      //  const double min_max_w);
      // Print contents of Tree
      //
      // Optional input:
      //   tab  tab to show depth
      inline void print(const char * tab="");
      // Determine max absolute winding number
      //
      // Inputs:
      //   p  query point 
      // Returns max winding number of 
      inline virtual double max_abs_winding_number(const Point & p) const; 
      // Same as above, but stronger assumptions on (V,F). Assumes (V,F) is a
      // simple polyhedron
      inline virtual double max_simple_abs_winding_number(const Point & p) const;
      // Compute or read cached winding number for point p with respect to mesh
      // in bounding box, recursing according to approximation criteria
      //
      // Inputs:
      //   p  query point 
      //   that  WindingNumberTree containing mesh w.r.t. which we're computing w.n.
      // Returns cached winding number
      inline virtual double cached_winding_number(const WindingNumberTree & that, const Point & p) const;
  };
}

// Implementation

#include "WindingNumberTree.h"
#include "winding_number.h"
#include "triangle_fan.h"
#include "exterior_edges.h"

#include <igl/PI.h>
#include <igl/remove_duplicate_vertices.h>

#include <iostream>
#include <limits>

//template <typename Point>
//WindingNumberMethod WindingNumberTree<Point>::method = EXACT_WINDING_NUMBER_METHOD;
//template <typename Point>
//double WindingNumberTree<Point>::min_max_w = 0;
template <typename Point>
std::map< std::pair<const igl::WindingNumberTree<Point>*,const igl::WindingNumberTree<Point>*>, double>
  igl::WindingNumberTree<Point>::cached;

template <typename Point>
inline igl::WindingNumberTree<Point>::WindingNumberTree():
  method(EXACT_WINDING_NUMBER_METHOD),
  parent(NULL),
  V(igl::dummyV),
  SV(),
  F(),
  //boundary(igl::boundary_facets<Eigen::MatrixXi,Eigen::MatrixXi>(F))
  cap(),
  radius(std::numeric_limits<double>::infinity()),
  center(0,0,0)
{
}

template <typename Point>
inline igl::WindingNumberTree<Point>::WindingNumberTree(
  const Eigen::MatrixXd & _V,
  const Eigen::MatrixXi & _F):
  method(EXACT_WINDING_NUMBER_METHOD),
  parent(NULL),
  V(igl::dummyV),
  SV(),
  F(),
  //boundary(igl::boundary_facets<Eigen::MatrixXi,Eigen::MatrixXi>(F))
  cap(),
  radius(std::numeric_limits<double>::infinity()),
  center(0,0,0)
{
  set_mesh(_V,_F);
}

template <typename Point>
inline void igl::WindingNumberTree<Point>::set_mesh(
    const Eigen::MatrixXd & _V,
    const Eigen::MatrixXi & _F)
{
  using namespace std;
  // Remove any exactly duplicate vertices
  // Q: Can this ever increase the complexity of the boundary?
  // Q: Would we gain even more by remove almost exactly duplicate vertices?
  Eigen::MatrixXi SF,SVI,SVJ;
  igl::remove_duplicate_vertices(_V,_F,0.0,SV,SVI,SVJ,F);
  triangle_fan(igl::exterior_edges(F),cap);
  V = SV;
}

template <typename Point>
inline igl::WindingNumberTree<Point>::WindingNumberTree(
  const igl::WindingNumberTree<Point> & parent,
  const Eigen::MatrixXi & _F):
  method(parent.method),
  parent(&parent),
  V(parent.V),
  SV(),
  F(_F),
  cap(triangle_fan(igl::exterior_edges(_F)))
{
}

template <typename Point>
inline igl::WindingNumberTree<Point>::~WindingNumberTree()
{
  delete_children();
}

template <typename Point>
inline void igl::WindingNumberTree<Point>::delete_children()
{
  using namespace std;
  // Delete children
  typename list<WindingNumberTree<Point>* >::iterator cit = children.begin();
  while(cit != children.end())
  {
    // clear the memory of this item
    delete (* cit);
    // erase from list, returns next element in iterator
    cit = children.erase(cit);
  }
}
      
template <typename Point>
inline void igl::WindingNumberTree<Point>::set_method(const WindingNumberMethod & m)
{
  this->method = m;
  for(auto child : children)
  {
    child->set_method(m);
  }
}

template <typename Point>
inline const Eigen::MatrixXd & igl::WindingNumberTree<Point>::getV() const
{
  return V;
}

template <typename Point>
inline const Eigen::MatrixXi & igl::WindingNumberTree<Point>::getF() const
{
  return F;
}

template <typename Point>
inline const Eigen::MatrixXi & igl::WindingNumberTree<Point>::getcap() const
{
  return cap;
}

template <typename Point>
inline void igl::WindingNumberTree<Point>::grow()
{
  // Don't grow
  return;
}

template <typename Point>
inline bool igl::WindingNumberTree<Point>::inside(const Point & /*p*/) const
{
  return true;
}

template <typename Point>
inline double igl::WindingNumberTree<Point>::winding_number(const Point & p) const
{
  using namespace std;
  //cout<<"+"<<boundary.rows();
  // If inside then we need to be careful
  if(inside(p))
  {
    // If not a leaf then recurse
    if(children.size()>0)
    {
      // Recurse on each child and accumulate
      double sum = 0;
      for(
        typename list<WindingNumberTree<Point>* >::const_iterator cit = children.begin();
        cit != children.end();
        cit++)
      {
        switch(method)
        {
          case EXACT_WINDING_NUMBER_METHOD:
            sum += (*cit)->winding_number(p);
            break;
          case APPROX_SIMPLE_WINDING_NUMBER_METHOD:
          case APPROX_CACHE_WINDING_NUMBER_METHOD:
            //if((*cit)->max_simple_abs_winding_number(p) > min_max_w)
            //{
              sum += (*cit)->winding_number(p);
            //}
            break;
          default:
            assert(false);
            break;
        }
      }
      return sum;
    }else
    {
      return winding_number_all(p);
    }
  }else{
    // Otherwise we can just consider boundary
    // Q: If we using the "multipole" method should we also subdivide the
    // boundary case?
    if((cap.rows() - 2) < F.rows())
    {
      switch(method)
      {
        case EXACT_WINDING_NUMBER_METHOD:
          return winding_number_boundary(p);
        case APPROX_SIMPLE_WINDING_NUMBER_METHOD:
        {
          double dist = (p-center).norm();
          // Radius is already an overestimate of inside
          if(dist>1.0*radius)
          {
            return 0;
          }else
          {
            return winding_number_boundary(p);
          }
        }
        case APPROX_CACHE_WINDING_NUMBER_METHOD:
        {
          return parent->cached_winding_number(*this,p);
        }
        default: assert(false);break;
      }
    }else
    {
      // doesn't pay off to use boundary
      return winding_number_all(p);
    }
  }
  return 0;
}

template <typename Point>
inline double igl::WindingNumberTree<Point>::winding_number_all(const Point & p) const
{
  double w = 0;
  igl::winding_number_3(
    V.data(),
    V.rows(),
    F.data(),
    F.rows(),
    p.data(),
    1,
    &w);
  return w;
}

template <typename Point>
inline double igl::WindingNumberTree<Point>::winding_number_boundary(const Point & p) const
{
  using namespace Eigen;
  using namespace std;

  double w = 0;
  // `cap` is already flipped inside out, so we don't need to flip sign of w
  igl::winding_number_3(
    V.data(),
    V.rows(),
    cap.data(),
    cap.rows(),
    &p[0],
    1,
    &w);
  return w;
}

//template <typename Point>
//inline double igl::WindingNumberTree<Point>::winding_number_approx_simple(
//  const Point & p, 
//  const double min_max_w)
//{
//  using namespace std;
//  if(max_simple_abs_winding_number(p) > min_max_w)
//  {
//    return winding_number(p);
//  }else
//  {
//    cout<<"Skipped! "<<max_simple_abs_winding_number(p)<<"<"<<min_max_w<<endl;
//    return 0;
//  }
//}

template <typename Point>
inline void igl::WindingNumberTree<Point>::print(const char * tab)
{
  using namespace std;
  // Print all facets
  cout<<tab<<"["<<endl<<F<<endl<<"]";
  // Print children
  for(
      typename list<WindingNumberTree<Point>* >::iterator cit = children.begin();
      cit != children.end();
      cit++)
  {
    cout<<","<<endl;
    (*cit)->print((string(tab)+"").c_str());
  }
}

template <typename Point>
inline double 
igl::WindingNumberTree<Point>::max_abs_winding_number(const Point & /*p*/) const
{
  return std::numeric_limits<double>::infinity();
}

template <typename Point>
inline double 
igl::WindingNumberTree<Point>::max_simple_abs_winding_number(
  const Point & /*p*/) const
{
  using namespace std;
  return numeric_limits<double>::infinity();
}

template <typename Point>
inline double igl::WindingNumberTree<Point>::cached_winding_number(
  const igl::WindingNumberTree<Point> & that,
  const Point & p) const
{
  using namespace std;
  // Simple metric for `is_far`
  //
  //   this             that
  //                   --------
  //   -----          /   |    \ .
  //  /  r  \        /    R     \ .
  // | p !   |      |     !      |
  //  \_____/        \          /
  //                  \________/
  //
  // 
  // a = angle formed by trapazoid formed by raising sides with lengths r and R
  // at respective centers.
  //
  // a = atan2(R-r,d), where d is the distance between centers

  // That should be bigger (what about parent? what about sister?)
  bool is_far = this->radius<that.radius;
  if(is_far)
  {
    double a = atan2(
      that.radius - this->radius,
      (that.center - this->center).norm());
    assert(a>0);
    is_far = (a<PI/8.0);
  }

  if(is_far)
  {
    // Not implemented yet
    pair<const WindingNumberTree*,const WindingNumberTree*> this_that(this,&that);
    // Need to compute it for first time?
    if(cached.count(this_that)==0)
    {
      cached[this_that] = 
        that.winding_number_boundary(this->center);
    }
    return cached[this_that];
  }else if(children.size() == 0)
  {
    // not far and hierarchy ended too soon: can't use cache
    return that.winding_number_boundary(p);
  }else
  {
    for(
      typename list<WindingNumberTree<Point>* >::const_iterator cit = children.begin();
      cit != children.end();
      cit++)
    {
      if((*cit)->inside(p))
      {
        return (*cit)->cached_winding_number(that,p);
      }
    }
    // Not inside any children? This can totally happen because bounding boxes
    // are set to bound contained facets. So sibilings may overlap and their
    // union may not contain their parent (though, their union is certainly a
    // subset of their parent).
    assert(false);
  }
  return 0;
}

// Explicit instanciation
//template class igl::WindingNumberTree<Eigen::Vector3d >;

#endif
