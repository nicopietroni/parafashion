#include "field_computation.h"

#include <wrap/io_trimesh/export_field.h>
#include <wrap/igl/smooth_field.h>
#include <mesh_type.h>

template <class TriMeshType>
void FieldComputation<TriMeshType>::ComputeField(TriMeshType &mesh,
                             FieldMode FMode,
                             ScalarType CurvatureFidelity,
                             bool DebugMsg)
{
    typedef FieldSmoother<TriMeshType> FieldSmootherType;
    typename FieldSmootherType::SmoothParam Param;
    Param.align_borders=true;
    if (FMode==FMBoundary)
    {
        //Param.alpha_curv=CurvatureFidelity;
        Param.SmoothM=SMNPoly;
    }
    else
    {
        Param.alpha_curv=CurvatureFidelity;
        Param.SmoothM=SMMiq;
    }

    //add features
    size_t NumFeatures=0;
    for (size_t i=0;i<mesh.face.size();i++)
    {
        for (size_t j=0;j<3;j++)
        {
            if (!mesh.face[i].IsFaceEdgeS(j))continue;
            CoordType Dir=mesh.face[i].P0(j)-mesh.face[i].P1(j);
            Dir.Normalize();
            Param.AddConstr.push_back(std::pair<int,CoordType>(i,Dir));
            NumFeatures++;
            break;//one constraint per face
        }
    }
    if (DebugMsg)
        std::cout<<"Added "<<NumFeatures<<" features"<<std::endl;

    FieldSmootherType::SmoothDirections(mesh,Param);

    vcg::tri::CrossField<TriMeshType>::UpdateSingularByCross(mesh);
    vcg::tri::CrossField<TriMeshType>::SetVertCrossVectorFromFace(mesh);
}

//Manual instantiation:
template class FieldComputation<CMesh>;