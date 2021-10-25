#ifndef FIELD_COMPUTATION
#define FIELD_COMPUTATION

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/simplex/face/topology.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
//#include <wrap/gl/trimesh.h>
#include <wrap/igl/smooth_field.h>
#include <wrap/io_trimesh/export_field.h>
#include <iostream>
#include <fstream>
#include <vcg/complex/algorithms/attribute_seam.h>
#include <vcg/complex/algorithms/crease_cut.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
enum FieldMode{FMBoundary,FMCurvature};

template <class TriMeshType>
class FieldComputation
{
    typedef typename TriMeshType::FaceType FaceType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;

public:

    static void ComputeField(TriMeshType &mesh,
                             FieldMode FMode,
                             ScalarType CurvatureFidelity,
                             bool DebugMsg=false)
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

};

#endif
