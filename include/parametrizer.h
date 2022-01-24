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
#include <param/multi_patch_param.h>
#include <param/metrics.h>

enum ParamMode{PMConformal,PMArap,PMCloth};

//enum ParamType{Arap,LSQMap};

//template <class TriMeshType>
//void ComputeUV(TriMeshType &mesh)
//{
//    if (PType==Arap)
//        vcg::tri::OptimizeUV_ARAP(mesh,100,0,true);
//    else
//        vcg::tri::OptimizeUV_LSCM(mesh,TriMeshType::VertexType::SELECTED);
//}
//template <class TriMeshType>
//void Decimate()
//{

//}

template <class TriMeshType>
bool ClothParametrize(TriMeshType &mesh,
                      std::vector<typename TriMeshType::ScalarType> &StretchU,
                      std::vector<typename TriMeshType::ScalarType> &StretchV,
                      typename TriMeshType::ScalarType error,
                      bool SelfIntCheck,
                      bool &DoIntersectUV)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    vcg::tri::MeshToMatrix<TriMeshType>::GetTriMeshData(mesh, F, V);
    
    ClothParam cloth(V, F, error);

    DoIntersectUV=false;

    if (SelfIntCheck)
        cloth.enableIntersectionCheck();
    else
        cloth.disableIntersectionCheck();

    bool success = cloth.paramAttempt();
    
    Eigen::MatrixXd V_uv;
    V_uv = cloth.getV2d();
    for (int i=0; i<(int)mesh.vert.size(); i++)
    {
        mesh.vert[i].T().P()[0] = V_uv(i,0);
        mesh.vert[i].T().P()[1] = V_uv(i,1);
    }

    Eigen::VectorXd stretch_u_vec,stretch_v_vec;
    measureStretchScore(V_uv,V, F,stretch_u_vec,stretch_v_vec);

    StretchU.clear();
    StretchV.clear();

    for (size_t i=0;i<mesh.face.size();i++)
    {
        StretchU.push_back(stretch_u_vec(i));
        StretchV.push_back(stretch_v_vec(i));
    }

    if (SelfIntCheck)
        DoIntersectUV=cloth.checkSelfIntersect();

    

    //    Eigen::VectorXd stretch_u;
    //    Eigen::VectorXd stretch_v;
    //    ClothParam.getStretchStats(stretch_u,stretch_v);
    //    StretchU.clear();StretchV.clear();
    //    for (size_t i=0;i<mesh.face.size();i++)
    //    {
    //        StretchU.push_back(stretch_u(i));
    //        StretchV.push_back(stretch_v(i));
    //    }
    return success;
}

template <class TriMeshType>
class Parametrizer
{
    typedef typename TriMeshType::FaceType FaceType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;
    typedef vcg::face::Pos<FaceType> PosType;
    
    //std::vector<typename TriMeshType::ScalarType> StretchU,StretchV;

    static void SetQasDistorsion()
    {
        //        assert(StretchU.size()==mesh.face.size());
        //        assert(StretchV.size()==mesh.face.size());
        //        for (size_t i=0;i<mesh.face.size();i++)
        //        {
        //            ScalarType diffU=fabs(StretchU[i]-1);
        //            ScalarType diffV=fabs(StretchV[i]-1);
        //            mesh.face[i].Q()=(diffU > diffV) ? StretchU[i] : StretchV[i];
        //        }
    }
    
    //    static void FindDartTipFromPos(TriMeshType &mesh,
    //                                   std::vector<std::vector<PosType> > &PosSeq,
    //                                   std::vector<bool> &IsDartTipVert)
    //    {
    //        //then count the number of times each sequence hits a vertex
    //        std::vector<size_t> EndPathNum(mesh.vert.size(),0);
    //        for (size_t i=0;i<PosSeq.size();i++)
    //            for (size_t j=0;j<PosSeq[i].size();j++)
    //            {
    //                PosType CurrP=PosSeq[i][j];
    //                if (CurrP.IsBorder())continue;

    //                VertexType *v0=CurrP.V();
    //                VertexType *v1=CurrP.VFlip();

    //                size_t IndexV0=vcg::tri::Index(mesh,v0);
    //                size_t IndexV1=vcg::tri::Index(mesh,v1);

    //                //count only once
    //                if (IndexV0>IndexV1)continue;

    //                EndPathNum[IndexV0]++;
    //                EndPathNum[IndexV1]++;
    //            }

    //        IsDartTipVert=std::vector<bool>(mesh.vert.size(),false);

    //        //then find the darts endpoint
    //        for (size_t i=0;i<EndPathNum.size();i++)
    //        {
    //            if (mesh.vert[i].IsB())continue;
    //            if (EndPathNum[i]!=1)continue;
    //            IsDartTipVert[i]=true;
    //        }
    //    }
    
    static void SplitIntoSubMeshes(TriMeshType &mesh,
                                   std::vector<TriMeshType*> &SubMeshes,
                                   std::vector<std::pair<int,int> >  &OriginalToSub)
    {
        //deallocate
        for (size_t i=0;i<SubMeshes.size();i++)
            delete(SubMeshes[i]);
        SubMeshes.clear();
        OriginalToSub.clear();
        
        //retrieve the partitions
        std::vector<std::vector<size_t> > Partitions;
        std::vector<size_t> StartF;
        for (size_t i=0;i<mesh.face.size();i++)
            StartF.push_back(i);

        RetrievePatchesFromSelEdges(mesh,StartF,Partitions);

        //save previous quality
        std::vector<ScalarType> OldQ;
        for (size_t i=0;i<mesh.face.size();i++)
            OldQ.push_back(mesh.face[i].Q());
        
        //set quality on original face
        for (size_t i=0;i<mesh.face.size();i++)
            mesh.face[i].Q()=i;
        
        OriginalToSub=std::vector<std::pair<int,int> >(mesh.face.size(),std::pair<int,int>(-1,-1));
        for (size_t i=0;i<Partitions.size();i++)
        {
            SubMeshes.push_back(new TriMeshType);
            PatchManager<TriMeshType>::GetMeshFromPatch(mesh,i,Partitions,(*SubMeshes.back()),true);
            (*SubMeshes.back()).UpdateAttributes();
            
            //            //then set the indexes
            //            OriginalFIdx.resize(OriginalFIdx.size()+1);
            for (size_t j=0;j<(*SubMeshes.back()).face.size();j++)
            {
                int IndexF=(*SubMeshes.back()).face[j].Q();
                assert(IndexF>=0);
                assert(IndexF<mesh.face.size());
                assert(OriginalToSub[IndexF]==(std::pair<int,int>(-1,-1)));
                OriginalToSub[IndexF]=std::pair<int,int>(i,j);
            }
        }
        
        //restore quality
        for (size_t i=0;i<mesh.face.size();i++)
            mesh.face[i].Q()=OldQ[i];
    }
    
    //    struct SeamData
    //    {
    //        int IndexSeam;

    //    };
    
    //    static void SplitForGlobalParam(TriMeshType &mesh,
    //                                    std::vector<TriMeshType*> &SubMeshes,
    //    {

    //    }

    static void MergeAcrossBoundarySeams(TriMeshType &mesh)
    {
        std::vector<std::vector<bool> > PerFaceEdgeB;
        PerFaceEdgeB.resize(mesh.face.size(),std::vector<bool>(3,false));
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
                PerFaceEdgeB[i][j]=vcg::face::IsBorder(mesh.face[i],j);

        vcg::tri::Clean<TriMeshType>::RemoveDuplicateVertex(mesh);
        vcg::tri::Clean<TriMeshType>::RemoveUnreferencedVertex(mesh);
        vcg::tri::Allocator<TriMeshType>::CompactEveryVector(mesh);
        mesh.UpdateAttributes();

        //if no border anymore means it is a seam
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                bool CurrIsB=vcg::face::IsBorder(mesh.face[i],j);
                if (PerFaceEdgeB[i][j]!=CurrIsB)
                    mesh.face[i].SetFaceEdgeS(j);
            }
    }

    static void GetSeamDataFromPos(TriMeshType &mesh,
                                   PosType &CurrPos,
                                   std::vector<TriMeshType*> &SubMeshes,
                                   const std::vector<std::pair<int,int> > &OriginalToSub,
                                   int &IndexSubM0,int &IndexSubM1,
                                   size_t &IndexV0,size_t &IndexV1,
                                   bool OtherV)
    {
        //get the pos and its opposite
        PosType OppPos=CurrPos;
        OppPos.FlipF();

        int IndexF0=vcg::tri::Index(mesh,CurrPos.F());
        int IndexF1=vcg::tri::Index(mesh,OppPos.F());

        //get the corresponding submesh
        IndexSubM0=OriginalToSub[IndexF0].first;
        IndexSubM1=OriginalToSub[IndexF1].first;

        TriMeshType *SubM0=SubMeshes[IndexSubM0];
        TriMeshType *SubM1=SubMeshes[IndexSubM1];

        //get the index of subfaces
        int IndexSubF0=OriginalToSub[IndexF0].second;
        int IndexSubF1=OriginalToSub[IndexF1].second;

        //check
        assert(IndexSubF0<(int)SubM0->face.size());
        assert(IndexSubF1<(int)SubM1->face.size());

        //get the position on the original mesh
        CoordType PosOrig=CurrPos.V()->P();
        if (OtherV)
            PosOrig=CurrPos.VFlip()->P();

        //get the vertex indexes on the sub and their position
        FaceType *SubF0=&SubM0->face[IndexSubF0];
        int IndexE0=CurrPos.E();

        //second position on the edge
        if (!OtherV)
            IndexV0=vcg::tri::Index(*SubM0,SubF0->V1(IndexE0));
        else
            IndexV0=vcg::tri::Index(*SubM0,SubF0->V0(IndexE0));

        assert(IndexV0<SubM0->vert.size());
        CoordType PosSub0=SubM0->vert[IndexV0].P();

        FaceType *SubF1=&SubM1->face[IndexSubF1];
        int IndexE1=OppPos.E();

        if (!OtherV)
            IndexV1=vcg::tri::Index(*SubM1,SubF1->V0(IndexE1));
        else
            IndexV1=vcg::tri::Index(*SubM1,SubF1->V1(IndexE1));

        assert(IndexV1<SubM1->vert.size());
        CoordType PosSub1=SubM1->vert[IndexV1].P();

        //check
        assert(PosOrig==PosSub0);
        assert(PosOrig==PosSub1);
    }

public:

    static void SplitForGlobalParam(TriMeshType &mesh,
                                    std::vector<TriMeshType*> &SubMeshes,
                                    std::vector<std::pair<int,int> > &OriginalToSub,
                                    std::vector<std::pair<int,int> > &MeshToMesh,
                                    std::vector<std::vector<std::pair<int,int> > > &VertToVert,
                                    std::vector<int > &DartTipVert)
    {
        OriginalToSub.clear();
        MeshToMesh.clear();
        VertToVert.clear();
        DartTipVert.clear();

        //get Pos Seq
        std::vector<std::vector<PosType> > PosSeq;
        RetrievePosSeqFromSelEdges(mesh,PosSeq);
        std::cout<<"There are "<<PosSeq.size()<<" seams"<<std::endl;

        //        //find the darts
        //        std::vector<bool> VertIsDartTip;
        //        FindDartTipFromPos(mesh,PosSeq,VertIsDartTip);
        //        assert(VertIsDartTip.size()==mesh.vert.size());
        //        std::set<CoordType> DartSet;
        //        for (size_t i=0;i<mesh.vert.size();i++)
        //        {
        //            if (!VertIsDartTip[i])continue;
        //            DartSet.insert(mesh.vert[i].P());
        //        }

        //split into submeshes
        SplitIntoSubMeshes(mesh,SubMeshes,OriginalToSub);

        //allocate
        MeshToMesh.resize(PosSeq.size(),std::pair<int,int>(-1,-1));
        VertToVert.resize(PosSeq.size());
        DartTipVert.resize(PosSeq.size(),-1);

        //then go over all the pos

        for (size_t i=0;i<PosSeq.size();i++)
        {
            assert(PosSeq[i].size()>0);

            PosType Pos0=PosSeq[i][0];

            int IndexSubM0,IndexSubM1;
            size_t IndexV0,IndexV1;
            GetSeamDataFromPos(mesh,Pos0,SubMeshes,
                               OriginalToSub,
                               IndexSubM0,IndexSubM1,
                               IndexV0,IndexV1,true);

            //            TriMeshType *SubM0=SubMeshes[IndexSubM0];
            //            TriMeshType *SubM1=SubMeshes[IndexSubM1];

            MeshToMesh[i]=std::pair<int,int>(IndexSubM0,IndexSubM1);

            //add the firs pair

            //check if is a Dart Tip
            if ((IndexSubM0==IndexSubM1)&&(IndexV0==IndexV1))
            {
                //CoordType Pos=SubM0->vert[IndexV0].P();
                //assert(DartSet.count(Pos)==1);
                DartTipVert[i]=IndexV0;
            }
            else
            {
                //add to the sequence at that point
                VertToVert[i].push_back(std::pair<int,int>(IndexV0,IndexV1));
            }

            for (size_t j=0;j<PosSeq[i].size();j++)
            {
                //get the pos and its opposite
                PosType CurrPos=PosSeq[i][j];
                assert(!CurrPos.IsBorder());
                int IndexSubMt0,IndexSubMt1;
                size_t IndexV0,IndexV1;
                GetSeamDataFromPos(mesh,CurrPos,SubMeshes,
                                   OriginalToSub,
                                   IndexSubMt0,IndexSubMt1,
                                   IndexV0,IndexV1,false);

                assert(IndexSubMt0==IndexSubM0);
                assert(IndexSubMt1==IndexSubM1);


                //check if is a Dart Tip
                if ((IndexSubM0==IndexSubM1)&&(IndexV0==IndexV1))
                {
                    //                    CoordType Pos=SubM0->vert[IndexV0].P();
                    //                    assert(DartSet.count(Pos)==1);
                    DartTipVert[i]=IndexV0;
                }
                else
                {
                    //add to the sequence at that point
                    VertToVert[i].push_back(std::pair<int,int>(IndexV0,IndexV1));
                }
            }
        }

        //final data check
        assert(VertToVert.size()==PosSeq.size());
        assert(VertToVert.size()==MeshToMesh.size());
        assert(DartTipVert.size()==MeshToMesh.size());

        for (size_t i=0;i<VertToVert.size();i++)
        {
            int IndexM0=MeshToMesh[i].first;
            int IndexM1=MeshToMesh[i].second;

            assert(IndexM0>=0);
            assert(IndexM0<(int)SubMeshes.size());

            assert(IndexM1>=0);
            assert(IndexM1<(int)SubMeshes.size());
            for (size_t j=0;j<VertToVert[i].size();j++)
            {
                int IndexV0=VertToVert[i][j].first;
                int IndexV1=VertToVert[i][j].second;

                assert(IndexV0>=0);
                assert(IndexV0<(int)SubMeshes[IndexM0]->vert.size());

                assert(IndexV1>=0);
                assert(IndexV1<(int)SubMeshes[IndexM1]->vert.size());

                CoordType Pos0=SubMeshes[IndexM0]->vert[IndexV0].P();
                CoordType Pos1=SubMeshes[IndexM1]->vert[IndexV1].P();
                assert(Pos0==Pos1);
            }
        }
    }

    
    //    static void Parametrize(TriMeshType &mesh,
    //                            ParamMode UVMode,
    //                            ScalarType BorderPatch=0)
    //    {
    //        std::vector<size_t> StartF;
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            StartF.push_back(i);

    //        std::vector<std::vector<size_t> > Partitions;
    //        RetrievePatchesFromSelEdges(mesh,StartF,Partitions);

    //        //save previous quality
    //        std::vector<ScalarType> OldQ;
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            OldQ.push_back(mesh.face[i].Q());

    //        //set quality on original face
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            mesh.face[i].Q()=i;

    //        std::vector<TriMeshType*> PatchMeshes;
    //        ScalarType A=0;
    //        for (size_t i=0;i<Partitions.size();i++)
    //        {
    //            PatchMeshes.push_back(new TriMeshType);
    //            PatchManager<TriMeshType>::GetMeshFromPatch(mesh,i,Partitions,(*PatchMeshes.back()),true);
    //            (*PatchMeshes.back()).UpdateAttributes();

    //            if (UVMode==PMConformal)
    //                vcg::tri::InitializeArapWithLSCM((*PatchMeshes.back()),0);
    //            if (UVMode==PMArap)
    //                vcg::tri::OptimizeUV_ARAP((*PatchMeshes.back()),100,0,true);
    //            if (UVMode==PMCloth)
    //                ClothParametrize<TriMeshType>((*PatchMeshes.back()),0.0);

    //            A=vcg::tri::UV_Utils<TriMeshType>::PerVertUVArea(*PatchMeshes.back());
    //        }
    //        ScalarType AbsDelta=math::Sqrt(A)*BorderPatch;
    //        PatchManager<TriMeshType>::ArrangeUVPatches(PatchMeshes,AbsDelta,false);

    //        for (size_t i=0;i<PatchMeshes.size();i++)
    //            for (size_t j=0;j<PatchMeshes[i]->face.size();j++)
    //            {
    //                size_t IndexFOrig=PatchMeshes[i]->face[j].Q();
    //                assert(IndexFOrig<mesh.face.size());
    //                for (size_t k=0;k<3;k++)
    //                {
    //                    vcg::Point2<ScalarType> UV=PatchMeshes[i]->face[j].V(k)->T().P();
    //                    mesh.face[IndexFOrig].WT(k).P()=UV;
    //                }
    //            }

    //        for (size_t i=0;i<PatchMeshes.size();i++)
    //            delete(PatchMeshes[i]);

    //        //restore quality
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            mesh.face[i].Q()=OldQ[i];
    //    }
    
    static void Parametrize(TriMeshType &mesh,
                            ParamMode UVMode,
                            std::vector<TriMeshType*> &SubMeshes,
                            std::vector<std::pair<int,int> > &MeshToMesh,
                            std::vector<std::vector<std::pair<int,int> > > &VertToVert,
                            std::vector<int> &DartTipVert,
                            bool continuity_seams,
                            bool continuity_darts,
                            ScalarType BorderPatch=0)
    {
        MergeAcrossBoundarySeams(mesh);
        //        std::vector<TriMeshType*> SubMeshes;
        std::vector<std::pair<int,int> > OriginalToSub;
        //        std::vector<std::pair<int,int> > MeshToMesh;
        //        std::vector<std::vector<std::pair<int,int> > > VertToVert;
        //        std::vector<int> DartTipVert;
        SplitForGlobalParam(mesh,SubMeshes,OriginalToSub,MeshToMesh,VertToVert,DartTipVert);

        ScalarType A=0;

        if (UVMode==PMConformal)
        {
            for (size_t i=0;i<SubMeshes.size();i++)
            {
                (*SubMeshes[i]).UpdateAttributes();
                vcg::tri::InitializeArapWithLSCM((*SubMeshes[i]),0);
                A+=vcg::tri::UV_Utils<TriMeshType>::PerVertUVArea(*SubMeshes[i]);
            }
        }

        if (UVMode==PMArap)
        {
            for (size_t i=0;i<SubMeshes.size();i++)
            {
                (*SubMeshes[i]).UpdateAttributes();
                vcg::tri::OptimizeUV_ARAP((*SubMeshes[i]),100,0,true);
                A+=vcg::tri::UV_Utils<TriMeshType>::PerVertUVArea(*SubMeshes[i]);
            }
        }

        if (UVMode==PMCloth)
        {

            //            for (size_t i=0;i<SubMeshes.size();i++)
            //            {
            //                (*SubMeshes[i]).UpdateAttributes();
            //                std::vector<typename TriMeshType::ScalarType> StretchU,StretchV;
            //                bool DoIntersectUV;
            //                ClothParametrize<TriMeshType>((*SubMeshes[i]),StretchU,StretchV,0.0,false,DoIntersectUV);
            //                A+=vcg::tri::UV_Utils<TriMeshType>::PerVertUVArea(*SubMeshes[i]);
            //            }

            //ge the seams
            //std::cout<<"0"<<std::endl;

            std::vector<Seam> seams;

            std::vector<std::vector<int> > vec_dart_tips;
            vec_dart_tips.resize(SubMeshes.size());
            std::vector<std::vector<std::vector<std::pair<int, int>>>> vec_dart_duplicates;
            vec_dart_duplicates.resize(SubMeshes.size());

            //std::cout<<"1"<<std::endl;

            for (size_t i=0;i<MeshToMesh.size();i++)
            {
                if (DartTipVert[i]!=-1)//in this case is a dart
                {
                    int IndexM0=MeshToMesh[i].first;
                    int IndexM1=MeshToMesh[i].second;
                    assert(IndexM0==IndexM1);
                    size_t IndexTip=DartTipVert[i];
                    if (continuity_darts)
                    {
                        vec_dart_tips[IndexM0].push_back(IndexTip);
                        vec_dart_duplicates[IndexM0].push_back(VertToVert[i]);
                    }
                }else//in this case is a dart
                {
                    Seam s;
                    s.patch1_id=MeshToMesh[i].first;
                    s.patch2_id=MeshToMesh[i].second;
                    s.corres=VertToVert[i];
                    if (continuity_seams)
                        seams.push_back(s);
                }

            }
            //std::cout<<"2"<<std::endl;

            //get the meshes
            std::vector<Eigen::MatrixXd> vec_V_3d;
            std::vector<Eigen::MatrixXi> vec_F;
            vec_V_3d.resize(SubMeshes.size());
            vec_F.resize(SubMeshes.size());


            for (size_t i=0;i<SubMeshes.size();i++)
            {
                vcg::tri::MeshToMatrix<TriMeshType>::GetTriMeshData(*SubMeshes[i], vec_F[i], vec_V_3d[i]);
            }


            std::cout<<"Global Param"<<std::endl;
            std::vector<Eigen::MatrixXd> vec_V_2d;
            finalParamMultiPatch(vec_V_3d, vec_F,vec_dart_duplicates,vec_dart_tips,seams,vec_V_2d);

            assert(vec_V_2d.size()==SubMeshes.size());
            for (int i=0; i<(int)SubMeshes.size(); i++)
            {
                for (int j=0; j<(int)SubMeshes[i]->vert.size(); j++)
                {
                    SubMeshes[i]->vert[j].T().P()[0] = vec_V_2d[i](j,0);
                    SubMeshes[i]->vert[j].T().P()[1] = vec_V_2d[i](j,1);
                }
                A+=vcg::tri::UV_Utils<TriMeshType>::PerVertUVArea(*SubMeshes[i]);
            }
        }

        ScalarType AbsDelta=math::Sqrt(A)*BorderPatch;
        PatchManager<TriMeshType>::ArrangeUVPatches(SubMeshes,AbsDelta,false);

        for (size_t i=0;i<mesh.face.size();i++)
        {
            size_t IndexSubM=OriginalToSub[i].first;
            size_t IndexFSub=OriginalToSub[i].second;
            assert(IndexSubM<SubMeshes.size());
            assert(IndexFSub<SubMeshes[IndexSubM]->face.size());
            for (size_t k=0;k<3;k++)
            {
                vcg::Point2<ScalarType> UV=SubMeshes[IndexSubM]->face[IndexFSub].V(k)->T().P();
                mesh.face[i].WT(k).P()=UV;
            }
        }

        //        for (size_t i=0;i<SubMeshes.size();i++)
        //            delete(SubMeshes[i]);
    }

    //    static void Parametrize(TriMeshType &mesh,
    //                            ParamMode UVMode,
    //                            ScalarType BorderPatch=0)
    //    {

    //        std::vector<size_t> StartF;
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            StartF.push_back(i);

    //        std::vector<std::vector<size_t> > Partitions;
    //        RetrievePatchesFromSelEdges(mesh,StartF,Partitions);

    //        //save previous quality
    //        std::vector<ScalarType> OldQ;
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            OldQ.push_back(mesh.face[i].Q());

    //        //set quality on original face
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            mesh.face[i].Q()=i;

    //        std::vector<TriMeshType*> PatchMeshes;
    //        ScalarType A=0;
    //        for (size_t i=0;i<Partitions.size();i++)
    //        {
    //            PatchMeshes.push_back(new TriMeshType);
    //            PatchManager<TriMeshType>::GetMeshFromPatch(mesh,i,Partitions,(*PatchMeshes.back()),true);
    //            (*PatchMeshes.back()).UpdateAttributes();

    //            if (UVMode==PMConformal)
    //                vcg::tri::InitializeArapWithLSCM((*PatchMeshes.back()),0);
    //            if (UVMode==PMArap)
    //                vcg::tri::OptimizeUV_ARAP((*PatchMeshes.back()),100,0,true);
    //            if (UVMode==PMCloth)
    //                ClothParametrize<TriMeshType>((*PatchMeshes.back()),0.0);

    //            A=vcg::tri::UV_Utils<TriMeshType>::PerVertUVArea(*PatchMeshes.back());
    //        }
    //        ScalarType AbsDelta=math::Sqrt(A)*BorderPatch;
    //        PatchManager<TriMeshType>::ArrangeUVPatches(PatchMeshes,AbsDelta,false);

    //        for (size_t i=0;i<PatchMeshes.size();i++)
    //            for (size_t j=0;j<PatchMeshes[i]->face.size();j++)
    //            {
    //                size_t IndexFOrig=PatchMeshes[i]->face[j].Q();
    //                assert(IndexFOrig<mesh.face.size());
    //                for (size_t k=0;k<3;k++)
    //                {
    //                    vcg::Point2<ScalarType> UV=PatchMeshes[i]->face[j].V(k)->T().P();
    //                    mesh.face[IndexFOrig].WT(k).P()=UV;
    //                }
    //            }

    //        for (size_t i=0;i<PatchMeshes.size();i++)
    //            delete(PatchMeshes[i]);

    //        //restore quality
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            mesh.face[i].Q()=OldQ[i];
    //    }
};

#endif
