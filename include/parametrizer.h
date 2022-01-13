#ifndef PARAMETRIZER
#define PARAMETRIZER

#include <wrap/igl/arap_parametrization.h>
#include <wrap/igl/lscm_parametrization.h>
#include <vcg/complex/algorithms/parametrization/distortion.h>
#include <vcg/space/outline2_packer.h>
#include <vcg/complex/algorithms/parametrization/uv_utils.h>
#include <tracing/patch_manager.h>
#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <param/cloth_param.h>

//enum ParamType{Arap,LSQMap};

//template <class TriMeshType>
//void ComputeUV(TriMeshType &mesh)
//{
//    if (PType==Arap)
//        vcg::tri::OptimizeUV_ARAP(mesh,100,0,true);
//    else
//        vcg::tri::OptimizeUV_LSCM(mesh,TriMeshType::VertexType::SELECTED);
//}

template <class TriMeshType>
bool ClothParametrize(TriMeshType &mesh,typename TriMeshType::ScalarType error=0.05)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    vcg::tri::MeshToMatrix<TriMeshType>::GetTriMeshData(mesh, F, V);

    ClothParam cloth(V, F, error);
    bool success = cloth.paramAttempt();

    if (!success) {
        cloth.printStretchStats();
    }
    
    Eigen::MatrixXd V_uv;
    V_uv = cloth.getV2d();

    for (int i=0; i<(int)mesh.vert.size(); i++)
    {
        mesh.vert[i].T().P()[0] = V_uv(i,0);
        mesh.vert[i].T().P()[1] = V_uv(i,1);
    }
    return success;
}

template <class TriMeshType>
class Parametrizer
{
    typedef typename TriMeshType::FaceType FaceType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;

public:

    static void Parametrize(TriMeshType &mesh,ScalarType BorderPatch=0)
    {
        std::vector<size_t> StartF;
        for (size_t i=0;i<mesh.face.size();i++)
            StartF.push_back(i);

        std::vector<std::vector<size_t> > Partitions;
        RetrievePatchesFromSelEdges(mesh,StartF,Partitions);

        //save previous quality
        std::vector<ScalarType> OldQ;
        for (size_t i=0;i<mesh.face.size();i++)
            OldQ.push_back(mesh.face[i].Q());

        //set quality on original face
        for (size_t i=0;i<mesh.face.size();i++)
            mesh.face[i].Q()=i;

        std::vector<TriMeshType*> PatchMeshes;
        for (size_t i=0;i<Partitions.size();i++)
        {
            PatchMeshes.push_back(new TriMeshType);
            PatchManager<TriMeshType>::GetMeshFromPatch(mesh,i,Partitions,(*PatchMeshes.back()),true);
            (*PatchMeshes.back()).UpdateAttributes();
            //vcg::tri::InitializeArapWithLSCM((*PatchMeshes.back()),0);
            //vcg::tri::OptimizeUV_ARAP((*PatchMeshes.back()),100,0,true);
            ClothParametrize<TriMeshType>((*PatchMeshes.back()),0.0);

//            vcg::Box2<ScalarType> uv_box=vcg::tri::UV_Utils<TriMeshType>::PerVertUVBox(*PatchMeshes.back());
//            std::cout<<"UV Box Dim Y:"<<uv_box.DimY()<<std::endl;
//            std::cout<<"UV Box Dim X:"<<uv_box.DimX()<<std::endl;

        }
        PatchManager<TriMeshType>::ArrangeUVPatches(PatchMeshes,BorderPatch,false);

        for (size_t i=0;i<PatchMeshes.size();i++)
            for (size_t j=0;j<PatchMeshes[i]->face.size();j++)
            {
                size_t IndexFOrig=PatchMeshes[i]->face[j].Q();
                assert(IndexFOrig<mesh.face.size());
                for (size_t k=0;k<3;k++)
                {
                    vcg::Point2<ScalarType> UV=PatchMeshes[i]->face[j].V(k)->T().P();
                    mesh.face[IndexFOrig].WT(k).P()=UV;
                }
            }

        for (size_t i=0;i<PatchMeshes.size();i++)
            delete(PatchMeshes[i]);

        //restore quality
        for (size_t i=0;i<mesh.face.size();i++)
            mesh.face[i].Q()=OldQ[i];
    }

};

#endif
