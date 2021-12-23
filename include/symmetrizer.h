#ifndef SYMMETRIZER
#define SYMMETRIZER

#include <vcg/space/plane3.h>
#include <vcg/space/fitting3.h>
#include <vcg/space/point_matching.h>
#include <vcg/math/matrix33.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/clean.h>
#include <tracing/patch_tracer.h>
#include <vcg/complex/algorithms/clean.h>

template <class TriMeshType>
void RefineMesh(TriMeshType &mesh,const vcg::Plane3<typename TriMeshType::ScalarType> &plane)
{
    //typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;
    typedef typename vcg::tri::QualityMidPointFunctor<TriMeshType> QualityMidP;
    typedef typename vcg::tri::QualityEdgePredicate<TriMeshType> QualityEdgeP;

    ScalarType Thr=0;//mesh.bbox.Diag()*0.001;

    vcg::tri::UpdateQuality<TriMeshType>::VertexFromPlane(mesh, plane);
    QualityMidP SlicingFunc((ScalarType)Thr);
    QualityEdgeP SlicingPred((ScalarType)Thr,0.0001);

    vcg::tri::RefineE<TriMeshType, QualityMidP, QualityEdgeP > (mesh, SlicingFunc, SlicingPred, false);

    mesh.UpdateAttributes();
}

template <class TriMeshType>
void MirrorMesh(TriMeshType &mesh,const vcg::Plane3<typename TriMeshType::ScalarType> &plane)
{
//    typedef typename TriMeshType::CoordType CoordType;
//    typedef typename TriMeshType::ScalarType ScalarType;
    typename TriMeshType::ScalarType alpha=0.00001;
    TriMeshType OtherSide;
    vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(OtherSide,mesh);
    for (size_t i=0;i<OtherSide.vert.size();i++)
    {
        //check if on the plane let it be
        typename TriMeshType::CoordType proj=plane.Projection(OtherSide.vert[i].P());
        if ((proj-OtherSide.vert[i].P()).Norm()<alpha)
            continue;
        OtherSide.vert[i].P()=plane.Mirror(OtherSide.vert[i].P());
    }

    for (size_t i=0;i<OtherSide.face.size();i++)
    {
        std::swap(OtherSide.face[i].V(0),OtherSide.face[i].V(1));
        bool IsSelE1=OtherSide.face[i].IsFaceEdgeS(1);
        bool IsSelE2=OtherSide.face[i].IsFaceEdgeS(2);

        if (IsSelE1)
            OtherSide.face[i].SetFaceEdgeS(2);
        else
            OtherSide.face[i].ClearFaceEdgeS(2);

        if (IsSelE2)
            OtherSide.face[i].SetFaceEdgeS(1);
        else
            OtherSide.face[i].ClearFaceEdgeS(1);
    }

    //invert face orientation on the other side
    vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(mesh,OtherSide);
    mesh.UpdateAttributes();

//    vcg::tri::Clean<TriMeshType>::RemoveDuplicateVertex(mesh);
//    mesh.UpdateAttributes();
//    vcg::tri::io::ExporterPLY<TriMeshType>::Save(mesh,"test.ply");
//    exit(0);
}

template <class TriMeshType>
void DeleteHalfMesh(TriMeshType &mesh,const vcg::Plane3<typename TriMeshType::ScalarType> &plane)
{
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;
    std::vector<vcg::Triangle3<ScalarType> > HalfFaces;
    for (size_t i=0;i<mesh.face.size();i++)
    {
        CoordType P0=mesh.face[i].P(0);
        CoordType P1=mesh.face[i].P(1);
        CoordType P2=mesh.face[i].P(2);
        CoordType AvgPos=(P0+P1+P2)/3;
        CoordType ProjPos=plane.Projection(AvgPos);
        CoordType Dir=AvgPos-ProjPos;
        Dir.Normalize();
        if ((plane.Direction()*Dir)<0)continue;
        HalfFaces.push_back(vcg::Triangle3<ScalarType>(P0,P1,P2));
    }
    mesh.Clear();
    for (size_t i=0;i<HalfFaces.size();i++)
        vcg::tri::Allocator<TriMeshType>::AddFace(mesh,HalfFaces[i].P(0),HalfFaces[i].P(1),HalfFaces[i].P(2));

    vcg::tri::Clean<TriMeshType>::RemoveDuplicateVertex(mesh);
    vcg::tri::Clean<TriMeshType>::RemoveUnreferencedVertex(mesh);
    vcg::tri::Allocator<TriMeshType>::CompactEveryVector(mesh);
}

template <class TriMeshType>
void SelectPatchBorderEdges(TriMeshType &mesh)
{
    vcg::tri::UpdateFlags<TriMeshType>::FaceClearFaceEdgeS(mesh);
    for (size_t i=0;i<mesh.face.size();i++)
    {
        size_t Q0=mesh.face[i].Q();
        for (size_t j=0;j<3;j++)
        {
            if (vcg::face::IsBorder(mesh.face[i],j))
            {
                mesh.face[i].SetFaceEdgeS(j);
                continue;
            }
            size_t Q1=mesh.face[i].FFp(j)->Q();
            if (Q0!=Q1)
                mesh.face[i].SetFaceEdgeS(j);
        }
    }
}

template <class TriMeshType>
void MakePartitionOnQConsistent(TriMeshType &mesh)
{
    std::vector<size_t> StartF;
    for (size_t i=0;i<mesh.face.size();i++)
        StartF.push_back(i);

    std::vector<std::vector<size_t> > Partitions;
    RetrievePatchesFromSelEdges(mesh,StartF,Partitions);

    for (size_t i=0;i<Partitions.size();i++)
        for (size_t j=0;j<Partitions[i].size();j++)
            mesh.face[Partitions[i][j]].Q()=i;

}

template <class TriMeshType>
void SymmetrizeMesh(TriMeshType &mesh,const vcg::Plane3<typename TriMeshType::ScalarType> &plane)
{
    RefineMesh(mesh,plane);
    DeleteHalfMesh(mesh,plane);
    MirrorMesh(mesh,plane);
}

template <class MeshType>
class SymmetrizeDeformation
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;

public:

    static vcg::Plane3<ScalarType> SymmetryPlane()
    {
        vcg::Plane3<ScalarType> Pl;
        Pl.Init(CoordType(0,0,0),CoordType(1,0,0));
        return Pl;
    }

    static CoordType SymmetricVect(const CoordType &vect)
    {
        vcg::Plane3<ScalarType> Pl=SymmetryPlane();
        Pl.SetOffset(0);
        CoordType mirroredV=Pl.Mirror(vect);
        return mirroredV;
    }

    static CoordType SymmetricPoint(const CoordType &pos)
    {
        vcg::Plane3<ScalarType> Pl=SymmetryPlane();
        CoordType mirrored=Pl.Mirror(pos);
        return mirrored;
    }

    static void RigidAlign(MeshType &Template,
                           MeshType &Deformed)
    {
        std::vector<CoordType> FixPoints,MovPoints;
        assert(Template.vert.size()==Deformed.vert.size());

        for (size_t i=0;i<Template.vert.size();i++)
        {
            FixPoints.push_back(Template.vert[i].P());
            MovPoints.push_back(Deformed.vert[i].P());
        }

        //then find the alignment
        vcg::Matrix44<ScalarType> Rigid;

        //compute rigid match
        vcg::ComputeRigidMatchMatrix<ScalarType>(FixPoints,MovPoints,Rigid);

        for (size_t i=0;i<Deformed.vert.size();i++)
            Deformed.vert[i].P()=Rigid*Deformed.vert[i].P();
    }

    static void FindSymmetricMap(MeshType &Template,
                                 std::vector<size_t> &OppFace,
                                 std::vector<CoordType> &OppBary)
    {
        vcg::GridStaticPtr<FaceType,ScalarType> MeshGrid;
        MeshGrid.Set(Template.face.begin(),Template.face.end());

        ScalarType MaxD=Template.bbox.Diag();
        ScalarType MinD;
        for (size_t i=0;i<Template.vert.size();i++)
        {
            CoordType Pos=Template.vert[i].P();
            CoordType MirrPos=SymmetricPoint(Pos);
            CoordType closestPt,normI,Bary;
            FaceType *f=vcg::tri::GetClosestFaceBase(Template,MeshGrid,MirrPos,
                                                     MaxD,MinD,closestPt,normI,Bary);
            assert(f!=NULL);
            size_t IndexF=vcg::tri::Index(Template,f);
            OppFace.push_back(IndexF);
            OppBary.push_back(Bary);
        }
    }

    static void SymmetrizationStep(MeshType &Template,
                                   MeshType &Deformed,
                                   const std::vector<size_t> &OppFace,
                                   const std::vector<CoordType> &OppBary,
                                   ScalarType Damp=0.5)
    {
        assert(Template.vert.size()==Deformed.vert.size());
        std::vector<CoordType> Displ0,Displ1;
        for (size_t i=0;i<Template.vert.size();i++)
            Displ0.push_back(Deformed.vert[i].P()-Template.vert[i].P());

        Displ1=Displ0;

        for (size_t i=0;i<Template.vert.size();i++)
        {
            size_t OppF=OppFace[i];
            assert(OppF<Template.face.size());
            CoordType OppB=OppBary[i];
            size_t IndexMirrV0=vcg::tri::Index(Template,Template.face[OppF].V(0));
            size_t IndexMirrV1=vcg::tri::Index(Template,Template.face[OppF].V(1));
            size_t IndexMirrV2=vcg::tri::Index(Template,Template.face[OppF].V(2));
            CoordType MirrorD= Displ0[IndexMirrV0]*OppB.X()+
                               Displ0[IndexMirrV1]*OppB.Y()+
                               Displ0[IndexMirrV2]*OppB.Z();
            MirrorD=SymmetricVect(MirrorD);
            Displ1[i]=Displ0[i]*Damp+MirrorD*(1-Damp);
        }

        for (size_t i=0;i<Deformed.vert.size();i++)
            Deformed.vert[i].P()=Template.vert[i].P()+Displ1[i];
    }

    static void SymmetrizeDistortion(MeshType &Template,
                                     MeshType &Deformed,
                                     size_t step_num=10,
                                     ScalarType Damp=0.5)
    {
        std::vector<size_t> OppFace;
        std::vector<CoordType> OppBary;
        FindSymmetricMap(Template,OppFace,OppBary);
        for (size_t i=0;i<step_num;i++)
            SymmetrizationStep(Template,Deformed,OppFace,OppBary,Damp);
    }

};

template <class TriMeshType>
class Symmetrizer
{
    typedef typename TriMeshType::FaceType FaceType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;

    TriMeshType &deformed_mesh;
    TriMeshType &reference_mesh;

//    std::vector<size_t> MapFace;
//    std::vector<CoordType> MapBary;

    //    CoordType SymmetricPos(const CoordType &Pos)
    //    {return(SymmPlane.Mirror(Pos));}


//    void ComputeMapPos()
//    {
//        MapBary.clear();
//        MapFace.clear();
//        vcg::GridStaticPtr<FaceType,ScalarType> grid;
//        grid.Set(reference_mesh.face.begin(),reference_mesh.face.end());
//        for (size_t i=0;i<reference_mesh.vert.size();i++)
//        {
//            //get the position
//            CoordType Pos=reference_mesh.vert[i].P();

//            //go on the other side
//            Pos=SymmPlane().Mirror(Pos);

//            //then find the map
//            ScalarType MaxD=reference_mesh.bbox.Diag();
//            ScalarType MinD;
//            CoordType closest,norm,bary;
//            FaceType *f=vcg::tri::GetClosestFaceBase(reference_mesh,grid,Pos,MaxD,MinD,closest,norm,bary);
//            size_t IndexF=vcg::tri::Index(reference_mesh,f);
//            for (size_t k=0;k<3;k++)
//                bary.V(k)=std::max((ScalarType)0,bary.V(k));
//            MapBary.push_back(bary);
//            MapFace.push_back(IndexF);
//        }
//    }

//    void SymmetrizeDeformationStep(std::vector<CoordType> &DefField)
//    {
//        assert(DefField.size()==reference_mesh.vert.size());

//        std::vector<CoordType> SumDispl(reference_mesh.vert.size(),CoordType(0,0,0));
//        std::vector<ScalarType> SumWeight(reference_mesh.vert.size(),0);

//        assert(MapBary.size()==reference_mesh.vert.size());
//        assert(MapFace.size()==reference_mesh.vert.size());
//        for (size_t i=0;i<reference_mesh.vert.size();i++)
//        {
//            CoordType VertDispl=DefField[i];
//            size_t IndexF=MapFace[i];
//            CoordType Bary=MapBary[i];
//            assert(IndexF<reference_mesh.face.size());
//            FaceType *f=&reference_mesh.face[IndexF];
//            for (size_t j=0;j<3;j++)
//            {
//                size_t IndexV=vcg::tri::Index(reference_mesh,f->V(j));
//                ScalarType W=Bary.V(j);
//                SumDispl[IndexV]+=(VertDispl*W);
//                SumWeight[IndexV]+=W;
//            }
//        }

//        for (size_t i=0;i<SumDispl.size();i++)
//        {
//            if (SumWeight[i]>0)
//                SumDispl[i]/=SumWeight[i];
//            else
//                SumDispl[i]=DefField[i];
//        }

//        for (size_t i=0;i<DefField.size();i++)
//            DefField[i]=DefField[i]*0.5+SumDispl[i]*0.5;
//    }

//    void ApplyDeformationField()
//    {
//        for (size_t i=0;i<deformed_mesh.vert.size();i++)
//            deformed_mesh.vert[i].P()=reference_mesh.vert[i].P()+DefField[i];
//    }

public:

    //std::vector<CoordType> DefField;
    //vcg::Plane3<ScalarType> SymmPlane;

    static vcg::Plane3<ScalarType> SymmPlane()
    {return SymmetrizeDeformation<TriMeshType>::SymmetryPlane();}

    //    void InitDeformedFromHalfMesh()
    //    {
    //        MirrorMesh(reference_mesh,SymmPlane);
    //    }

//    void MirrorReferenceMesh()
//    {
//        //RefineMesh(reference_mesh,SymmPlane);
//        //DeleteHalfMesh(reference_mesh,SymmPlane);
//        //MirrorMesh(reference_mesh,SymmPlane);
//        SymmetrizeMesh(reference_mesh,SymmPlane);
//    }

//    void ComputeSymmetricDeformationField(size_t num_steps=10)
//    {
//        ComputeMapPos();
//        DefField.clear();
//        DefField.resize(reference_mesh.vert.size(),CoordType(0,0,0));
//        for (size_t i=0;i<reference_mesh.vert.size();i++)
//            DefField[i]=deformed_mesh.vert[i].P()-reference_mesh.vert[i].P();

//        for (size_t s=0;s<num_steps;s++)
//            SymmetrizeDeformationStep(DefField);
//    }

    void SymmetrizeDeformedMesh()
    {
        SymmetrizeDeformation<TriMeshType>::RigidAlign(reference_mesh,deformed_mesh);
        //RigidAlign();

//        ComputeSymmetricDeformationField(0);

//        ApplyDeformationField();
//        vcg::tri::io::ExporterPLY<TriMeshType>::Save(deformed_mesh,"test0.ply");
        SymmetrizeDeformation<TriMeshType>::SymmetrizeDistortion(reference_mesh,deformed_mesh,50);
//        vcg::tri::io::ExporterPLY<TriMeshType>::Save(deformed_mesh,"test1.ply");

        //then apply the deformation
        SymmetrizeMesh(deformed_mesh,SymmPlane());

        vcg::tri::Clean<TriMeshType>::RemoveUnreferencedVertex(deformed_mesh);
        vcg::tri::Allocator<TriMeshType>::CompactEveryVector(deformed_mesh);

        deformed_mesh.UpdateAttributes();

        for (size_t i=0;i<deformed_mesh.face.size()/2;i++)
            deformed_mesh.face[i].C()=vcg::Color4b(72,209,204,255);

        for (size_t i=deformed_mesh.face.size()/2;i<deformed_mesh.face.size();i++)
            deformed_mesh.face[i].C()=vcg::Color4b(220,220,220,255);

    }

//    //void Symmetrize
//    void RigidAlign()
//    {
//        assert(deformed_mesh.vert.size()==reference_mesh.vert.size());

//        //then find the alignment
//        vcg::Matrix44<ScalarType> Rigid;

//        //compute rigid match
//        std::vector<CoordType> FixPoints,MovPoints;
//        for (size_t i=0;i<deformed_mesh.vert.size();i++)
//        {
//            FixPoints.push_back(reference_mesh.vert[i].P());
//            MovPoints.push_back(deformed_mesh.vert[i].P());
//        }
//        vcg::ComputeRigidMatchMatrix<ScalarType>(FixPoints,MovPoints,Rigid);

//        //then apply
//        for (size_t i=0;i<deformed_mesh.vert.size();i++)
//            deformed_mesh.vert[i].P()=Rigid*deformed_mesh.vert[i].P();

//        deformed_mesh.UpdateAttributes();
//    }

    //    void ComputeSymmetricPlane()
    //    {

    //    }

//    void Init()
//    {
//        SymmPlane.Init(deformed_mesh.bbox.Center(),CoordType(1,0,0));
//        //        ExtrinsicPlaneSymmetry<TriMeshType> ExtrSymm(reference_mesh);
//        //        ExtrSymm.Init();
//        //        std::vector<vcg::Plane3<ScalarType> > Planes;
//        //        ExtrSymm.GetPlanes(Planes,1);
//        //        SymmPlane=Planes[0];
//    }

    //    void SortDefFacesFacesBySymmetry()
    //    {
    //       DeleteHalfMesh(deformed_mesh,SymmPlane);
    //       MirrorMesh(deformed_mesh,SymmPlane);
    //    }

    void GetHalfDefMesh(TriMeshType &half_def_mesh)
    {
        vcg::tri::UpdateFlags<TriMeshType>::FaceClearS(deformed_mesh);
        for (size_t i=0;i<deformed_mesh.face.size()/2;i++)
            deformed_mesh.face[i].SetS();

        half_def_mesh.Clear();
        vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(half_def_mesh,deformed_mesh,true);
        half_def_mesh.UpdateAttributes();
    }

    void CopyFromHalfDefMesh(TriMeshType &half_def_mesh)
    {
        deformed_mesh.Clear();
        vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(deformed_mesh,half_def_mesh);
        MirrorMesh(deformed_mesh,SymmPlane());
    }

    void CopyFieldFromHalfDefMesh(TriMeshType &half_def_mesh)
    {
        assert(half_def_mesh.face.size()==deformed_mesh.face.size()/2);
        vcg::tri::UpdateFlags<TriMeshType>::FaceClearS(deformed_mesh);
        for (size_t i=0;i<deformed_mesh.face.size()/2;i++)
        {
            size_t OffsetF=deformed_mesh.face.size()/2;
            CoordType PD1=half_def_mesh.face[i].PD1();
            CoordType PD2=half_def_mesh.face[i].PD2();
            deformed_mesh.face[i].PD1()=PD1;
            deformed_mesh.face[i].PD2()=PD2;
            deformed_mesh.face[OffsetF+i].PD1()=SymmPlane().Mirror(PD1);
            deformed_mesh.face[OffsetF+i].PD2()=SymmPlane().Mirror(PD2);
        }
        vcg::tri::CrossField<TriMeshType>::UpdateSingularByCross(deformed_mesh);
        vcg::tri::CrossField<TriMeshType>::SetVertCrossVectorFromFace(deformed_mesh);

    }

    void CopyPropertiesFromHalfDefMesh(TriMeshType &half_def_mesh)
    {
        typename TriMeshType::ScalarType alpha=0.00001;

        assert(half_def_mesh.face.size()==deformed_mesh.face.size()/2);
        assert(half_def_mesh.vert.size()==deformed_mesh.vert.size()/2);
        vcg::tri::UpdateFlags<TriMeshType>::FaceClearS(deformed_mesh);
        for (size_t i=0;i<deformed_mesh.face.size()/2;i++)
        {
            size_t OffsetF=deformed_mesh.face.size()/2;

            deformed_mesh.face[i].Q()=half_def_mesh.face[i].Q();
            deformed_mesh.face[OffsetF+i].Q()=half_def_mesh.face[i].Q();

            deformed_mesh.face[i].C()=half_def_mesh.face[i].C();
            deformed_mesh.face[OffsetF+i].C()=half_def_mesh.face[i].C();
        }

        for (size_t i=0;i<deformed_mesh.vert.size()/2;i++)
        {
            size_t OffsetF=deformed_mesh.vert.size()/2;
            CoordType Pos=half_def_mesh.vert[i].P();
            deformed_mesh.vert[i].P()=Pos;

            typename TriMeshType::CoordType proj=SymmPlane().Projection(Pos);
            if ((proj-Pos).Norm()<alpha)
                deformed_mesh.vert[OffsetF+i].P()=Pos;
            else
                deformed_mesh.vert[OffsetF+i].P()=SymmPlane().Mirror(Pos);
        }

        //then copy the selected edges
        vcg::tri::UpdateFlags<TriMeshType>::FaceClearFaceEdgeS(deformed_mesh);
        for (size_t i=0;i<deformed_mesh.face.size()/2;i++)
        {
            for (size_t j=0;j<3;j++)
            {
                bool IsSelE=half_def_mesh.face[i].IsFaceEdgeS(j);
                if (!IsSelE)continue;

                deformed_mesh.face[i].SetFaceEdgeS(j);

                size_t OtherE=0;
                if (j==1)OtherE=2;
                if (j==2)OtherE=1;

                size_t OffsetF=deformed_mesh.face.size()/2;
                deformed_mesh.face[OffsetF+i].SetFaceEdgeS(OtherE);
            }
        }

    }

    Symmetrizer(TriMeshType &_deformed_mesh,
                TriMeshType &_reference_mesh):deformed_mesh(_deformed_mesh),reference_mesh(_reference_mesh)
    {}

};

#endif
