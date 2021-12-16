#include "peel_outer_hull_layers.h"
#include "../per_face_normals.h"
#include "outer_hull.h"
#include <vector>
#include <iostream>
//#define IGL_PEEL_OUTER_HULL_LAYERS_DEBUG
#ifdef IGL_PEEL_OUTER_HULL_LAYERS_DEBUG
#include "../writePLY.h"
#include "../writeDMAT.h"
#include "../STR.h"
#endif

using namespace std;
template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedN,
  typename Derivedodd,
  typename Derivedflip>
IGL_INLINE size_t igl::peel_outer_hull_layers(
  const Eigen::PlainObjectBase<DerivedV > & V,
  const Eigen::PlainObjectBase<DerivedF > & F,
  const Eigen::PlainObjectBase<DerivedN > & N,
  Eigen::PlainObjectBase<Derivedodd > & odd,
  Eigen::PlainObjectBase<Derivedflip > & flip)
{
  using namespace Eigen;
  using namespace std;
  typedef typename DerivedF::Index Index;
  typedef Matrix<typename DerivedF::Scalar,Dynamic,DerivedF::ColsAtCompileTime> MatrixXF;
  typedef Matrix<typename DerivedN::Scalar,Dynamic,DerivedN::ColsAtCompileTime> MatrixXN;
  typedef Matrix<Index,Dynamic,1> MatrixXI;
  typedef Matrix<typename Derivedflip::Scalar,Dynamic,Derivedflip::ColsAtCompileTime> MatrixXflip;
  const Index m = F.rows();
#ifdef IGL_PEEL_OUTER_HULL_LAYERS_DEBUG
  cout<<"peel outer hull layers..."<<endl;
#endif
#ifdef IGL_PEEL_OUTER_HULL_LAYERS_DEBUG
  cout<<"calling outer hull..."<<endl;
  writePLY(STR("peel-outer-hull-input.ply"),V,F);
#endif

#ifdef IGL_PEEL_OUTER_HULL_LAYERS_DEBUG
  cout<<"resize output ..."<<endl;
#endif
  // keep track of iteration parity and whether flipped in hull
  MatrixXF Fr = F;
  MatrixXN Nr = N;
  odd.resize(m,1);
  flip.resize(m,1);
  // Keep track of index map
  MatrixXI IM = MatrixXI::LinSpaced(m,0,m-1);
  // This is O(n * layers)
  bool odd_iter = true;
  MatrixXI P(m,1);
  Index iter = 0;
  while(Fr.size() > 0)
  {
    assert(Fr.rows() == IM.rows());
    // Compute outer hull of current Fr
    MatrixXF Fo;
    MatrixXI Jo;
    MatrixXflip flipr;
#ifdef IGL_PEEL_OUTER_HULL_LAYERS_DEBUG
  cout<<"calling outer hull..."<<endl;
  writePLY(STR("outer-hull-input-"<<iter<<".ply"),V,Fr);
  writeDMAT(STR("outer-hull-input-"<<iter<<".dmat"),Nr);
#endif
    outer_hull(V,Fr,Nr,Fo,Jo,flipr);
#ifdef IGL_PEEL_OUTER_HULL_LAYERS_DEBUG
  writePLY(STR("outer-hull-output-"<<iter<<".ply"),V,Fo);
  cout<<"reindex, flip..."<<endl;
#endif
    assert(Fo.rows() == Jo.rows());
    // all faces in Fo of Fr
    vector<bool> in_outer(Fr.rows(),false);
    for(Index g = 0;g<Jo.rows();g++)
    {
      odd(IM(Jo(g))) = odd_iter;
      P(IM(Jo(g))) = iter;
      in_outer[Jo(g)] = true;
      flip(IM(Jo(g))) = flipr(Jo(g));
    }
    // Fr = Fr - Fo
    // update IM
    MatrixXF prev_Fr = Fr;
    MatrixXN prev_Nr = Nr;
    MatrixXI prev_IM = IM;
    Fr.resize(prev_Fr.rows() - Fo.rows(),F.cols());
    Nr.resize(Fr.rows(),3);
    IM.resize(Fr.rows());
    {
      Index g = 0;
      for(Index f = 0;f<prev_Fr.rows();f++)
      {
        if(!in_outer[f])
        {
          Fr.row(g) = prev_Fr.row(f);
          Nr.row(g) = prev_Nr.row(f);
          IM(g) = prev_IM(f);
          g++;
        }
      }
    }
    odd_iter = !odd_iter;
    iter++;
  }
  return iter;
}

template <
  typename DerivedV,
  typename DerivedF,
  typename Derivedodd,
  typename Derivedflip>
IGL_INLINE size_t igl::peel_outer_hull_layers(
  const Eigen::PlainObjectBase<DerivedV > & V,
  const Eigen::PlainObjectBase<DerivedF > & F,
  Eigen::PlainObjectBase<Derivedodd > & odd,
  Eigen::PlainObjectBase<Derivedflip > & flip)
{
  Eigen::Matrix<typename DerivedV::Scalar,DerivedF::RowsAtCompileTime,3> N;
  per_face_normals(V,F,N);
  return peel_outer_hull_layers(V,F,N,odd,flip);
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template size_t igl::peel_outer_hull_layers<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<bool, -1, 1, 0, -1, 1>, Eigen::Matrix<bool, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&);
template size_t igl::peel_outer_hull_layers<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<bool, -1, 1, 0, -1, 1>, Eigen::Matrix<bool, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&);
#endif
