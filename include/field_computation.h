#ifndef FIELD_COMPUTATION
#define FIELD_COMPUTATION

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/simplex/face/topology.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
//#include <wrap/gl/trimesh.h>
#include <wrap/io_trimesh/export_field.h>
#include <iostream>
#include <fstream>
#include <vcg/complex/algorithms/attribute_seam.h>
//#include <vcg/complex/algorithms/crease_cut.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <animation_manager.h>

enum FieldMode{FMBoundary,FMCurvature,FMCurvatureFrames,FMCurvatureOnly};

template <class TriMeshType>
class FieldComputation
{
    typedef typename TriMeshType::FaceType FaceType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;

public:

    static void ComputeField(TriMeshType &mesh,
                             AnimationManager<TriMeshType> &AManag,
                             FieldMode &FMode,
                             ScalarType CurvatureFidelity,
                             bool DebugMsg=false);

};

#include "field_computation.h"
#include "smooth_field_directional.h"
#include <wrap/io_trimesh/export_field.h>
#include <tracing/mesh_type.h>

//template <class TriMeshType>
//void FieldComputation<TriMeshType>::ComputeField(TriMeshType &mesh,
//                             FieldMode FMode,
//                             ScalarType CurvatureFidelity,
//                             bool DebugMsg)
//{
//    //typedef vcg::tri::FieldSmoother<TriMeshType> FieldSmootherType;
//    typedef DirectionalFieldSmoother<TriMeshType> FieldSmootherType;
//    typename FieldSmootherType::SmoothParam Param;
//    Param.align_borders=true;
////    if (FMode==FMBoundary)
////    {
////        //Param.alpha_curv=CurvatureFidelity;
////        Param.SmoothM=vcg::tri::SMNPoly;
////    }
////    else
////    {
////        Param.alpha_curv=CurvatureFidelity;
////        Param.SmoothM=vcg::tri::SMMiq;
////    }

//    //add features
//    size_t NumFeatures=0;
//    for (size_t i=0;i<mesh.face.size();i++)
//    {
//        for (size_t j=0;j<3;j++)
//        {
//            if (!mesh.face[i].IsFaceEdgeS(j))continue;
//            CoordType Dir=mesh.face[i].P0(j)-mesh.face[i].P1(j);
//            Dir.Normalize();
//            Param.AddConstr.push_back(std::pair<int,CoordType>(i,Dir));
//            NumFeatures++;
//            break;//one constraint per face
//        }
//    }
//    if (DebugMsg)
//        std::cout<<"Added "<<NumFeatures<<" features"<<std::endl;

//    FieldSmootherType::SmoothDirections(mesh,Param);

//    vcg::tri::CrossField<TriMeshType>::UpdateSingularByCross(mesh);
//    vcg::tri::CrossField<TriMeshType>::SetVertCrossVectorFromFace(mesh);
//}


template <class TriMeshType>
void FieldComputation<TriMeshType>::ComputeField(TriMeshType &mesh,
                                                 AnimationManager<TriMeshType> &AManag,
                                                 FieldMode &FMode,
                                                 ScalarType SmoothFactor,
                                                 bool DebugMsg)
{

    //typedef vcg::tri::FieldSmoother<TriMeshType> FieldSmootherType;
    typedef DirectionalFieldSmoother<TriMeshType> FieldSmootherType;
    typename FieldSmootherType::SmoothParam Param;
    Param.align_borders=true;

    if (FMode==FMCurvatureOnly)
    {
        FieldSmootherType::InitByCurvature(mesh,3);
        vcg::tri::UpdateQuality<TriMeshType>::VertexFromFace(mesh);
        //vcg::tri::UpdateColor<TriMeshType>::PerFaceConstant(mesh,vcg::Color4b::LightGray);
        vcg::tri::UpdateColor<TriMeshType>::PerFaceQualityRamp(mesh);
//        for (size_t i=0;i<mesh.face.size();i++)
//        {
//            mesh.face[i].PD2()*=mesh.face[i].Q();
//            mesh.face[i].PD1()*=mesh.face[i].Q();
//        }
        return;
    }

    if (FMode==FMCurvatureFrames)
    {
        if (AManag.NumFrames()<2)
        {
            std::cout<<"WARNING ONLY ONE FRAME DEFINED"<<std::endl;
            FMode=FMCurvature;
        }
        else
        {

            AManag.UpdateAnimationMesh();
            //AManag.UpdateProjectionBasis();
            AManag.InitTargetDirectionsOnMesh();
            AManag.TransferDirOnMesh(mesh);
            Param.use_predefined_field=true;
            Param.alpha_curv=SmoothFactor;
            Param.curv_thr=0.5;
        }
    }
    if (FMode==FMBoundary)
    {
        Param.alpha_curv=0;
        Param.curv_thr=0;
    }
    if (FMode==FMCurvature)
    {
        Param.alpha_curv=SmoothFactor;
        Param.curv_thr=0.3;
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

    //AManag.UpdateAnimationMesh();
    //AManag.UpdateProjectionBasisIfNeeded();
}

//Manual instantiation:
template class FieldComputation<TraceMesh>;


#endif
