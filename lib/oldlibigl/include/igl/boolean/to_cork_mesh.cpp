#ifndef IGL_NO_CORK
#include "to_cork_mesh.h"
template <
  typename DerivedV,
  typename DerivedF>
IGL_INLINE void igl::to_cork_mesh(
  const Eigen::PlainObjectBase<DerivedV > & V,
  const Eigen::PlainObjectBase<DerivedF > & F,
  CorkTriMesh & mesh)
{
  using namespace std;
  assert((F.cols() == 0 || F.cols() == 3) && "Facets should be triangles.");
  assert((V.cols() == 0 || V.cols() == 3) && "Vertices should be in 3D.");
  mesh.n_triangles = F.rows();
  mesh.n_vertices = V.rows();
  mesh.vertices = new double[mesh.n_vertices*3];
  mesh.triangles = new uint[mesh.n_triangles*3];
  for(size_t v = 0;v<mesh.n_vertices;v++)
  {
    for(size_t c = 0;c<3;c++)
    {
      mesh.vertices[v*3+c] = V(v,c);
    }
  }
  for(size_t f = 0;f<mesh.n_triangles;f++)
  {
    for(size_t c = 0;c<3;c++)
    {
      mesh.triangles[f*3+c] = F(f,c);
    }
  }
}
#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::to_cork_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, CorkTriMesh&);
template void igl::to_cork_mesh<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, CorkTriMesh&);
#endif
#endif
