#ifndef IGL_NO_CORK
#include "mesh_boolean_cork.h"
#include "to_cork_mesh.h"
#include "from_cork_mesh.h"

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC>
IGL_INLINE void igl::mesh_boolean_cork(
  const Eigen::PlainObjectBase<DerivedVA > & VA,
  const Eigen::PlainObjectBase<DerivedFA > & FA,
  const Eigen::PlainObjectBase<DerivedVB > & VB,
  const Eigen::PlainObjectBase<DerivedFB > & FB,
  const MeshBooleanType & type,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC)
{
  CorkTriMesh A,B,C;
  // pointer to output so it's easy to redirect on degenerate cases
  CorkTriMesh *ret = &C;
  to_cork_mesh(VA,FA,A);
  to_cork_mesh(VB,FB,B);
  switch(type)
  {
    case MESH_BOOLEAN_TYPE_UNION:
      if(A.n_triangles == 0)
      {
        ret = &B;
      }else if(B.n_triangles == 0)
      {
        ret = &A;
      }else
      {
        computeUnion(A,B,ret);
      }
      break;
    case MESH_BOOLEAN_TYPE_INTERSECT:
      if(A.n_triangles == 0 || B.n_triangles == 0)
      {
        ret->n_triangles = 0;
        ret->n_vertices = 0;
      }else
      {
        computeIntersection(A,B,ret);
      }
      break;
    case MESH_BOOLEAN_TYPE_MINUS:
      if(A.n_triangles == 0)
      {
        ret->n_triangles = 0;
        ret->n_vertices = 0;
      }else if(B.n_triangles == 0)
      {
        ret = &A;
      }else
      {
        computeDifference(A,B,ret);
      }
      break;
    case MESH_BOOLEAN_TYPE_XOR:
      if(A.n_triangles == 0)
      {
        ret = &B;
      }else if(B.n_triangles == 0)
      {
        ret = &A;
      }else
      {
        computeSymmetricDifference(A,B,&C);
      }
      break;
    case MESH_BOOLEAN_TYPE_RESOLVE:
      resolveIntersections(A,B,&C);
      break;
    default:
      assert(false && "Unknown type");
      return;
  }
  from_cork_mesh(*ret,VC,FC);
  freeCorkTriMesh(&A);
  freeCorkTriMesh(&B);
  freeCorkTriMesh(&C);
}
#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::mesh_boolean_cork<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::mesh_boolean_cork<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, igl::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&);
#endif

#endif
