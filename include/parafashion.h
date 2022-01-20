#pragma once

#include "symmetrizer.h"
#include "parametrizer.h"
#include "field_computation.h"
#include "trace_path_ui.h"
#include <tracing/tracer_interface.h>
#include <tracing/mesh_type.h>
#include "animation_manager.h"
#include "vcg/complex/algorithms/parametrization/uv_utils.h"

#include <vcg/complex/algorithms/isotropic_remeshing.h>
#include <vcg/complex/algorithms/crease_cut.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>

#include <igl/principal_curvature.h>
#include <vcg/complex/algorithms/update/color.h>

#define PRINT_PARAFASHION_TIMING 
#ifdef PRINT_PARAFASHION_TIMING
#include <chrono>
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;
#endif

enum PatchMode{PMMinTJuncions,PMAvgTJuncions,PMAllTJuncions};

template <class TriMeshType>
bool HasDart(TriMeshType &m)
{
    TriMeshType m1;
    vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(m1,m);
    vcg::tri::Clean<TriMeshType>::RemoveDuplicateVertex(m1);
    vcg::tri::Allocator<TriMeshType>::CompactEveryVector(m1);
    return(m1.vert.size()!=m.vert.size());
}

template <class TriMeshType>
void PreProcessMeshForRemeshing(TriMeshType &m,
                                bool merge,
                                typename TriMeshType::ScalarType damp=0.5,
                                size_t smooth_step=2)
{
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;

    //first select all borders
    typename TriMeshType::ScalarType AvgL=0;
    size_t Num=0;
    for (size_t i=0;i<m.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            AvgL+=(m.face[i].P0(j)-m.face[i].P1(j)).Norm();
            Num++;
            if (!vcg::face::IsBorder(m.face[i],j))continue;
            m.face[i].SetFaceEdgeS(j);
        }

    //merge vertices
    if (merge)
        vcg::tri::Clean<TriMeshType>::RemoveDuplicateVertex(m);

    //then smooth only selected vertices
    for (size_t s=0;s<smooth_step;s++)
    {
        std::vector<CoordType> PrevPos(m.vert.size(),CoordType(0,0,0));
        for (size_t i=0;i<m.vert.size();i++)
            PrevPos[i]=m.vert[i].P();

        //smooth the border first
        std::vector<CoordType> AvgPos(m.vert.size(),CoordType(0,0,0));
        std::vector<size_t> Num(m.vert.size(),0);
        vcg::tri::UpdateSelection<TriMeshType>::VertexAll(m);
        for (size_t i=0;i<m.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!m.face[i].IsFaceEdgeS(j))continue;
                size_t IndexV0=vcg::tri::Index(m,m.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(m,m.face[i].V1(j));
                m.vert[IndexV0].ClearS();
                m.vert[IndexV1].ClearS();
                CoordType P0=m.face[i].P0(j);
                CoordType P1=m.face[i].P1(j);
                AvgPos[IndexV0]+=P1;
                AvgPos[IndexV1]+=P0;
                Num[IndexV0]++;
                Num[IndexV1]++;
            }
        //average borders
        for (size_t i=0;i<m.vert.size();i++)
        {
            if (Num[i]==0)continue;

            assert(!m.vert[i].IsS());

            CoordType TargetP=AvgPos[i]/Num[i];
            m.vert[i].P()=TargetP;//m.vert[i].P()*0.5+TargetP*0.5;
        }
        //smooth the rest
        vcg::tri::Smooth<TriMeshType>::VertexCoordLaplacian(m,1,true);

        //apply damp
        for (size_t i=0;i<m.vert.size();i++)
            m.vert[i].P()=PrevPos[i]*damp+m.vert[i].P()*(1-damp);
    }
}

template <class TriMeshType>
void Remesh(TriMeshType &m,typename TriMeshType::ScalarType reduction=2)
{
    typedef typename TriMeshType::ScalarType ScalarType;

    PreProcessMeshForRemeshing(m,true);

    //    bool hasdart=HasDart(m);
    //    if (hasdart)
    //        vcg::tri::io::ExporterPLY<TriMeshType>::Save(m,"before_rem.ply");

    //get average edge
    typename TriMeshType::ScalarType AvgL=0;
    size_t Num=0;
    for (size_t i=0;i<m.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            AvgL+=(m.face[i].P0(j)-m.face[i].P1(j)).Norm();
            Num++;
        }
    AvgL/=Num;

    //then remesh
    typename vcg::tri::IsotropicRemeshing<TriMeshType>::Params para;
    para.iter = 5;
    para.splitFlag    = false;
    para.swapFlag     = true;
    para.collapseFlag = true;
    para.smoothFlag   = true;
    para.projectFlag  = false;
    para.selectedOnly = false;
    para.adapt=false;
    //para.aspectRatioThr = 0.3;


    //para.maxSurfDist = m.bbox.Diag() / 2500.;
    para.surfDistCheck = false;
    para.userSelectedCreases = true;


    ScalarType edgeL = AvgL*reduction;
    para.SetTargetLen(edgeL);

    //std::cout << "Before Remeshing - faces: " << m.FN() << " quality: " <<  computeAR(m) << std::endl;
    vcg::tri::IsotropicRemeshing<TriMeshType>::Do(m, para);

    //then re-cut
    //vcg::tri::CutMeshAlongSelectedFaceEdges(m);
    VertSplitter<TriMeshType>::SplitAlongEdgeSel(m);

    vcg::tri::Clean<TriMeshType>::RemoveUnreferencedVertex(m);
    vcg::tri::Allocator<TriMeshType>::CompactEveryVector(m);

    m.UpdateAttributes();

    //    //save test if has dart
    //    if (hasdart)
    //    {
    //        vcg::tri::io::ExporterPLY<TriMeshType>::Save(m,"after_rem.ply");
    //        exit(0);
    //    }
}

template <class TriMeshType>
void RemeshByDeci(TriMeshType &m,typename TriMeshType::ScalarType reduction=2)
{
    typedef typename TriMeshType::ScalarType ScalarType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename vcg::tri::BasicVertexPair<VertexType> VertPair;

    class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapse< TriMeshType, VertPair, MyTriEdgeCollapse>{
    public:
        inline MyTriEdgeCollapse(const VertPair &p, int mark, BaseParameterClass *pp)
        {
            this->localMark = mark;
            this->pos=p;
            this->_priority = this->ComputePriority(pp);
        }
    };

    PreProcessMeshForRemeshing(m,false);

    vcg::BaseParameterClass BClass;
    vcg::LocalOptimization<TriMeshType> DeciSession(m,&BClass);
    DeciSession.template Init<MyTriEdgeCollapse >();

    //DeciSession.SetTargetSimplices(target_faces);
    size_t targetV=m.vert.size()/reduction;
    std::cout<<"Target:"<<targetV<<std::endl;
    //DeciSession.SetTargetVertices(targetV);
    m.UpdateAttributes();
    while(DeciSession.DoOptimization() && m.vn>targetV ){};

    //    //then re-cut
    //    vcg::tri::CutMeshAlongSelectedFaceEdges(m);

    vcg::tri::Clean<TriMeshType>::RemoveUnreferencedVertex(m);
    vcg::tri::Allocator<TriMeshType>::CompactEveryVector(m);

    m.UpdateAttributes();
}

template <class TriMeshType>
class Parafashion
{
    typedef typename TriMeshType::FaceType FaceType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;
    typedef typename vcg::Point2<ScalarType> UVCoordType;

    TriMeshType &deformed_mesh;
    TriMeshType &reference_mesh;
    AnimationManager<TriMeshType> &AManag;

    TriMeshType half_def_mesh;

    TriMeshType deformed_mesh_step0;
    TriMeshType reference_mesh_step0;

    TriMeshType deformed_mesh_step1;
    TriMeshType half_def_mesh_step1;

    TriMeshType deformed_mesh_step2;
    TriMeshType half_def_mesh_step2;

    std::vector<TriMeshType*> SubMeshes;
    std::vector<std::pair<int,int> > MeshToMesh;
    std::vector<std::vector<std::pair<int,int> > > VertToVert;
    std::vector<int> DartTipVert;

public:

    PatchMode PMode;
    FieldMode FMode;
    ParamMode UVMode;

    bool match_valence;
    bool check_stress;
    bool use_darts;
    bool allow_self_glue;
    bool remove_along_symmetry;
    ScalarType param_boundary;
    size_t max_corners;
    ScalarType max_compression;
    ScalarType max_tension;
    bool remesh_on_test;
    bool CheckUVIntersection;
    bool SmoothBeforeRemove;
    std::vector<typename TraceMesh::CoordType> PatchCornerPos;
    PriorityMode PrioMode;
    bool continuity_seams;
    bool continuity_darts;

    void CleanMeshAttributes()
    {
        vcg::tri::Allocator<TriMeshType>::DeletePerVertexAttribute(half_def_mesh,std::string("Singular"));
        vcg::tri::Allocator<TriMeshType>::DeletePerVertexAttribute(half_def_mesh,std::string("SingularIndex"));

        vcg::tri::Allocator<TriMeshType>::DeletePerVertexAttribute(deformed_mesh,std::string("Singular"));
        vcg::tri::Allocator<TriMeshType>::DeletePerVertexAttribute(deformed_mesh,std::string("SingularIndex"));

        vcg::tri::Allocator<TriMeshType>::DeletePerVertexAttribute(reference_mesh,std::string("Singular"));
        vcg::tri::Allocator<TriMeshType>::DeletePerVertexAttribute(reference_mesh,std::string("SingularIndex"));

    }

    void RestoreInitMesh()
    {
        CleanMeshAttributes();

        deformed_mesh.Clear();
        reference_mesh.Clear();
        vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(deformed_mesh,deformed_mesh_step0);
        vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(reference_mesh,reference_mesh_step0);
        deformed_mesh.InitRPos();
        reference_mesh.InitRPos();
        deformed_mesh.UpdateAttributes();
        reference_mesh.UpdateAttributes();
    }

    template <class MeshType>
    class MeshArapQuality
    {
        typedef typename MeshType::VertexType VertexType;
        typedef typename MeshType::ScalarType  ScalarType;
        typedef typename MeshType::FacePointer FacePointer;


    public:

        static ParamMode &UVMode()
        {
            static ParamMode CurrUVMode=PMConformal;
            return CurrUVMode;
        }

        //{PMConformal,PMArap,PMCloth};

        static ScalarType & MinQ()
        {
            static ScalarType MinV=-0.05;
            return MinV;
        }

        static ScalarType & MaxQ()
        {
            static ScalarType MaxV=0.05;
            return MaxV;
        }

        static bool & RemeshOnTest()
        {
            static bool rem=false;
            return rem;
        }

        static bool &ContinuousCheckSelfInt()
        {
            static bool check=false;
            return check;
        }

        MeshArapQuality(){ContinuousCheckSelfInt()=true;}

        ScalarType operator()(MeshType &m) const
        {

#ifdef PRINT_PARAFASHION_TIMING
            steady_clock::time_point pre_param = steady_clock::now();
#endif
            if (RemeshOnTest())
                Remesh(m);
            //RemeshByDeci(m);

            //            size_t numH=vcg::tri::Clean<TriMeshType>::CountHoles(m);
            //            if (numH!=1)
            //            {
            //                vcg::tri::io::ExporterPLY<TriMeshType>::Save(m,"test_holes.ply");
            //                assert(0);
            //            }

            //assert(numH==1);

            //,PMArap,
            if (UVMode()==PMCloth)
            {
                bool success = ClothParametrize<TriMeshType>(m, MaxQ(),ContinuousCheckSelfInt()); // quality-check param, NOT the final one you see on screen
#ifdef PRINT_PARAFASHION_TIMING
                steady_clock::time_point post_param = steady_clock::now();
                int param_time = duration_cast<microseconds>(post_param - pre_param).count();
                std::cout << "Param time : " << param_time << " [µs]" << std::endl;
#endif
                if (success)
                    return 0;
                else
                    return 1000000.0;
            }

            if (UVMode()==PMConformal)
            {
                vcg::tri::InitializeArapWithLSCM(m,0);
#ifdef PRINT_PARAFASHION_TIMING
                steady_clock::time_point post_param = steady_clock::now();
                int param_time = duration_cast<microseconds>(post_param - pre_param).count();
                std::cout << "Param time : " << param_time << " [µs]" << std::endl;
#endif
                vcg::tri::Distortion<TraceMesh,false>::SetQasDistorsion(m,vcg::tri::Distortion<TraceMesh,false>::EdgeComprStretch);
                ScalarType A=0;
                for (size_t i=0;i<m.face.size();i++)
                {
                    if (m.face[i].Q()<(MinQ()))A+=vcg::DoubleArea(m.face[i]);
                    if (m.face[i].Q()>MaxQ())A+=vcg::DoubleArea(m.face[i]);
                }
                return A;
            }

            if (UVMode()==PMArap)
            {
                vcg::tri::InitializeArapWithLSCM(m,0);
                //vcg::tri::OptimizeUV_ARAP(m,5,0,true);
#ifdef PRINT_PARAFASHION_TIMING
                steady_clock::time_point post_param = steady_clock::now();
                int param_time = duration_cast<microseconds>(post_param - pre_param).count();
                std::cout << "Param time : " << param_time << " [µs]" << std::endl;
#endif
                vcg::tri::Distortion<TraceMesh,false>::SetQasDistorsion(m,vcg::tri::Distortion<TraceMesh,false>::EdgeComprStretch);
                ScalarType A=0;
                for (size_t i=0;i<m.face.size();i++)
                {
                    if (m.face[i].Q()<(MinQ()))A+=vcg::DoubleArea(m.face[i]);
                    if (m.face[i].Q()>MaxQ())A+=vcg::DoubleArea(m.face[i]);
                }
                return A;
            }
            //vcg::tri::OptimizeUV_ARAP(m,5,0,true);
            //vcg::tri::InitializeArapWithLSCM(m,0);
            
            //evaluate the distortion
            //            vcg::tri::Distortion<TraceMesh,false>::SetQasDistorsion(m,vcg::tri::Distortion<TraceMesh,false>::EdgeComprStretch);
            //            ScalarType A=0;
            //            for (size_t i=0;i<m.face.size();i++)
            //            {
            //                if (m.face[i].Q()<(MinQ()))A+=vcg::DoubleArea(m.face[i]);
            //                if (m.face[i].Q()>MaxQ())A+=vcg::DoubleArea(m.face[i]);
            //            }
            //return (A/SumA);

            //return A;
        }
    };

    //    void TracePatch(bool SaveStep=true,
    //                    bool DebugMSG=false)
    //    {
    //        if (SaveStep)
    //        {
    //            deformed_mesh.Clear();
    //            half_def_mesh.Clear();
    //            vcg::tri::Append<TraceMesh,TraceMesh>::Mesh(deformed_mesh,deformed_mesh_step2);
    //            vcg::tri::Append<TraceMesh,TraceMesh>::Mesh(half_def_mesh,half_def_mesh_step2);
    //            deformed_mesh.UpdateAttributes();
    //            half_def_mesh.UpdateAttributes();
    //            half_def_mesh.InitRPos();
    //        }


    //        typedef PatchTracer<TriMeshType,MeshArapQuality<TriMeshType> > PTracerType;

    //        half_def_mesh.UpdateAttributes();

    //        VertexFieldGraph<TriMeshType> VGraph(half_def_mesh);
    //        VGraph.InitGraph(DebugMSG);



    //        PTracerType PTr(VGraph);
    //        PTr.InitTracer(100,DebugMSG);

    //        PTr.MaxVal=max_corners;
    //        PatchGeneralParameters::MaxAdmittable()=max_corners;

    //        PTr.split_on_removal=true;
    //        PTr.away_from_singular=true;
    //        PTr.match_valence=match_valence;
    //        PTr.CClarkability=-1;
    //        half_def_mesh.UpdateAttributes();
    //        bool PreRemoveStep=true;

    //        if (PMode==PMMinTJuncions)
    //        {
    //            PreRemoveStep=true;
    //            PTr.split_on_removal=false;
    //        }
    //        if (PMode==PMAvgTJuncions)
    //        {
    //            PreRemoveStep=true;
    //            PTr.split_on_removal=true;
    //        }
    //        if (PMode==PMAllTJuncions)
    //        {
    //            PreRemoveStep=false;
    //            PTr.split_on_removal=true;
    //        }

    //        if (max_compression<max_tension)
    //        {
    //            MeshArapQuality<TriMeshType>::MaxQ()=max_tension;
    //            MeshArapQuality<TriMeshType>::MinQ()=max_compression;
    //            PTr.check_quality_functor=check_stress;
    //        }

    //        //RecursiveProcess<PTracerType>(PTr,100,true,true,PreRemoveStep,false,false,false,DebugMSG);

    //        //RecursiveProcessForTexturingWithDarts<PTracerType>(PTr,100,true,true,PreRemoveStep,false,false,false,DebugMSG);
    //        RecursiveProcessForTexturingWithDarts<PTracerType>(PTr,100,true,true,PreRemoveStep,false,false,false,DebugMSG);


    //        std::cout<<"Remaining "<<PTr.ChoosenPaths.size()<<" paths"<<std::endl;
    //        //copy partition on quality
    //        for (size_t i=0;i<PTr.FacePartitions.size();i++)
    //            half_def_mesh.face[i].Q()=PTr.FacePartitions[i];

    //        //PTr.ColorByPartitions();
    //        Symmetrizer<TriMeshType> Symm(deformed_mesh,reference_mesh);
    //        Symm.CopyPropertiesFromHalfDefMesh(half_def_mesh);

    //        deformed_mesh.UpdateAttributes();
    //        SelectPatchBorderEdges(deformed_mesh);
    //        MakePartitionOnQConsistent(deformed_mesh);
    //        deformed_mesh.ScatterColorByQualityFace();

    //        PTr.GetVisualCornersPos(PatchCornerPos);

    //    }

    void TracePatch(bool SaveStep=true,
                    bool DebugMSG=false)
    {
        if (SaveStep)
        {
            deformed_mesh.Clear();
            half_def_mesh.Clear();
            vcg::tri::Append<TraceMesh,TraceMesh>::Mesh(deformed_mesh,deformed_mesh_step2);
            vcg::tri::Append<TraceMesh,TraceMesh>::Mesh(half_def_mesh,half_def_mesh_step2);
            deformed_mesh.UpdateAttributes();
            half_def_mesh.UpdateAttributes();
            half_def_mesh.InitRPos();
        }

        typedef PatchTracer<TriMeshType,MeshArapQuality<TriMeshType> > PTracerType;

        half_def_mesh.UpdateAttributes();

        VertexFieldGraph<TriMeshType> VGraph(half_def_mesh);
        VGraph.InitGraph(DebugMSG);


        PTracerType PTr(VGraph);

        PTr.MaxVal=max_corners;
        PatchGeneralParameters::MaxAdmittable()=max_corners;

        PTr.split_on_removal=true;
        PTr.away_from_singular=true;
        PTr.match_valence=match_valence;
        PTr.CClarkability=-1;
        //PTr.FirstBorder=true;
        half_def_mesh.UpdateAttributes();
        bool PreRemoveStep=true;
        PTr.PrioMode=PrioMode;
        PTr.CheckUVIntersection=CheckUVIntersection;

        if (max_compression<max_tension)
        {
            MeshArapQuality<TriMeshType>::UVMode()=UVMode;
            MeshArapQuality<TriMeshType>::MaxQ()=max_tension;
            MeshArapQuality<TriMeshType>::MinQ()=max_compression;
            MeshArapQuality<TriMeshType>::RemeshOnTest()=remesh_on_test;
            PTr.check_quality_functor=check_stress;
        }

        if (PMode==PMMinTJuncions)
        {
            PreRemoveStep=true;
            PTr.split_on_removal=false;
        }
        if (PMode==PMAvgTJuncions)
        {
            PreRemoveStep=true;
            PTr.split_on_removal=true;
        }
        if (PMode==PMAllTJuncions)
        {
            PreRemoveStep=false;
            PTr.split_on_removal=true;
        }


        PTr.InitTracer(100,DebugMSG);
        PTr.AllowRemoveConcave=true;
        //PTr.away_from_singular=false;
        //PTr.AllowRemoveConcave=false;

        std::vector<ScalarType> DartPriority;
        if (use_darts)
            GetVertPriorityByConvexity(DartPriority);

        bool only_needed=true;
        bool check_smooth_folds=false;

        if ((!use_darts)&&(!allow_self_glue))
            RecursiveProcess<PTracerType>(PTr,100,only_needed,true,PreRemoveStep,false,false,check_smooth_folds,SmoothBeforeRemove,DebugMSG);

        if ((use_darts)&&(!allow_self_glue))
            RecursiveProcessWithDarts<PTracerType>(PTr,100,only_needed,true,PreRemoveStep,false,false,check_smooth_folds,DartPriority,SmoothBeforeRemove,DebugMSG);

        if ((!use_darts)&&(allow_self_glue))
            RecursiveProcessForTexturing<PTracerType>(PTr,100,only_needed,true,PreRemoveStep,false,false,check_smooth_folds,SmoothBeforeRemove,DebugMSG);

        if ((use_darts)&&(allow_self_glue))
            RecursiveProcessForTexturingWithDarts<PTracerType>(PTr,100,only_needed,true,PreRemoveStep,false,false,check_smooth_folds,DartPriority,SmoothBeforeRemove,DebugMSG);
        //RecursiveProcessForTexturingWithDarts<PTracerType>(PTr,100,true,false,false,false,false,false,DebugMSG);


        //be sure to have selected all the paths
        std::vector<std::vector<vcg::face::Pos<FaceType> > > PathPos;
        vcg::tri::UpdateFlags<TriMeshType>::FaceClearFaceEdgeS(half_def_mesh);
        GetPathPos(VGraph,PTr.ChoosenPaths,PathPos);
        half_def_mesh.SelectPos(PathPos,true);

        //        std::cout<<"Remaining "<<PTr.ChoosenPaths.size()<<" paths"<<std::endl;
        //        //copy partition on quality
        //        for (size_t i=0;i<PTr.FacePartitions.size();i++)
        //            half_def_mesh.face[i].Q()=PTr.FacePartitions[i];

        //PTr.ColorByPartitions();

        //then copy everything
        Symmetrizer<TriMeshType> Symm(deformed_mesh,reference_mesh);
        Symm.CopyPropertiesFromHalfDefMesh(half_def_mesh);


        deformed_mesh.UpdateAttributes();
        //SelectPatchBorderEdges(deformed_mesh);
        MakePartitionOnQConsistent(deformed_mesh);
        deformed_mesh.ScatterColorByQualityFace();

        PTr.GetVisualCornersPos(PatchCornerPos);

    }

    void  ComputeField(bool SaveStep=true)
    {

        if (SaveStep)
        {
            deformed_mesh.Clear();
            half_def_mesh.Clear();
            vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(deformed_mesh,deformed_mesh_step1);
            vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(half_def_mesh,half_def_mesh_step1);
            deformed_mesh.UpdateAttributes();
            half_def_mesh.UpdateAttributes();
        }

        //compute field on half mesh
        half_def_mesh.UpdateSharpFeaturesFromSelection();

        FieldComputation<TriMeshType>::ComputeField(half_def_mesh,AManag,FMode,1000);//,align_border);
        vcg::tri::CrossField<TriMeshType>::UpdateSingularByCross(half_def_mesh);
        vcg::tri::CrossField<TriMeshType>::SetVertCrossVectorFromFace(half_def_mesh);
        half_def_mesh.InitSingVert();
        half_def_mesh.InitRPos();


        Symmetrizer<TriMeshType> Symm(deformed_mesh,reference_mesh);

        //        //THEN REMOVE IF THEY ARE SOFT
        //        vcg::tri::UpdateFlags<TriMeshType>::FaceClearFaceEdgeS(half_def_mesh);
        //        half_def_mesh.UpdateSharpFeaturesFromSelection();

        PreProcessMesh(half_def_mesh,false);
        Symm.CopyFromHalfDefMesh(half_def_mesh);

        //CHECK
        deformed_mesh.UpdateSharpFeaturesFromSelection();

        Symm.CopyFieldFromHalfDefMesh(half_def_mesh);


        if (SaveStep)
        {
            half_def_mesh_step2.Clear();
            deformed_mesh_step2.Clear();
            vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(half_def_mesh_step2,half_def_mesh);
            vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(deformed_mesh_step2,deformed_mesh);
        }
    }

    void AddSharpConstraints(const std::vector<std::vector<CoordType> > &PickedPoints)
    {
        //get the half mesh
        //std::cout<<"1"<<std::endl;
        Symmetrizer<TriMeshType> Symm(deformed_mesh,reference_mesh);
        Symm.GetHalfDefMesh(half_def_mesh);
        //std::cout<<"2"<<std::endl;
        PathUI<TriMeshType> TPath(half_def_mesh);
        //std::cout<<"Size V:"<<PickedPoints.size()<<std::endl;
        TPath.AddSharpConstraints(PickedPoints);
        //        //std::cout<<"3"<<std::endl;
        //        PreProcessMesh(half_def_mesh,true);
        //std::cout<<"4"<<std::endl;
        //        Symm.CopyFromHalfDefMesh(half_def_mesh);

        //        //CHECK
        //        //half_def_mesh.InitRPos();
        //        //std::cout<<"5"<<std::endl;
        //        deformed_mesh.UpdateSharpFeaturesFromSelection();
    }

    void MakeMeshSymmetric(const std::vector<std::vector<CoordType> > &PickedPoints,
                           bool SaveStep=true)
    {

        //make the mesh symmetric
        Symmetrizer<TriMeshType> Symm(deformed_mesh,reference_mesh);
        Symm.SymmetrizeDeformedMesh();

        AddSharpConstraints(PickedPoints);

        if (SaveStep)
        {
            half_def_mesh_step1.Clear();
            deformed_mesh_step1.Clear();
            vcg::tri::Append<TraceMesh,TraceMesh>::Mesh(half_def_mesh_step1,half_def_mesh);
            vcg::tri::Append<TraceMesh,TraceMesh>::Mesh(deformed_mesh_step1,deformed_mesh);
        }
    }

    void GetUVSeamsPolylines( std::vector<std::vector<UVCoordType> > &UVPolyL,
                              std::vector<vcg::Color4b> &Color,
                              std::vector<UVCoordType> &Dots,
                              vcg::Color4b &ColorDots)
    {
        UVPolyL.clear();
        Color.clear();
        Dots.clear();
        ColorDots=vcg::Color4b::Red;

        for (size_t i=0;i<VertToVert.size();i++)
        {
            TriMeshType *M0=SubMeshes[MeshToMesh[i].first];
            TriMeshType *M1=SubMeshes[MeshToMesh[i].second];

            vcg::Color4b Col=vcg::Color4b::Scatter(VertToVert.size(),i,0.5f,0.95f);
            Color.push_back(Col);
            UVPolyL.resize(UVPolyL.size()+1);
            for (size_t j=0;j<VertToVert[i].size();j++)
            {
                size_t IndexV0=VertToVert[i][j].first;
                UVCoordType UV0=M0->vert[IndexV0].T().P();
                UVPolyL.back().push_back(UV0);
            }
            UVPolyL.resize(UVPolyL.size()+1);
            Color.push_back(Col);
            for (size_t j=0;j<VertToVert[i].size();j++)
            {
                size_t IndexV1=VertToVert[i][j].second;
                UVCoordType UV1=M1->vert[IndexV1].T().P();
                UVPolyL.back().push_back(UV1);
            }

            if (DartTipVert[i]!=-1)
            {
                size_t IndexDart=DartTipVert[i];
                UVCoordType UV0=M0->vert[IndexDart].T().P();
                Dots.push_back(UV0);
            }
        }
    }

    void DoParametrize()
    {
        Parametrizer<TriMeshType>::Parametrize(deformed_mesh,UVMode,
                                               SubMeshes, MeshToMesh,
                                               VertToVert,DartTipVert,
                                               continuity_seams,
                                               continuity_darts,
                                               param_boundary);
        //parametrized=true;
    }

    void Init()
    {
        deformed_mesh_step0.Clear();
        reference_mesh_step0.Clear();
        vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(deformed_mesh_step0,deformed_mesh);
        vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(reference_mesh_step0,reference_mesh);
    }

    void SaveDebugPatches(const std::string &ProjPath)
    {
        std::vector<size_t> StartF;
        for (size_t i=0;i<deformed_mesh.face.size();i++)
            StartF.push_back(i);

        std::vector<std::vector<size_t> > Partitions;
        RetrievePatchesFromSelEdges(deformed_mesh,StartF,Partitions);

        for (size_t i=0;i<Partitions.size();i++)
        {
            std::string numPatch=std::to_string(i);
            std::string Name_3D_mesh=ProjPath+"patch_3D_"+numPatch+".obj";
            std::string Name_UV_mesh=ProjPath+"patch_UV_"+numPatch+".obj";
            std::cout<<"Saving:"<<Name_3D_mesh.c_str()<<std::endl;
            std::cout<<"Saving:"<<Name_UV_mesh.c_str()<<std::endl;

            TriMeshType curr_patch_3D,curr_patch_UV;
            PatchManager<TriMeshType>::GetMeshFromPatch(deformed_mesh,i,Partitions,curr_patch_3D,true);
            curr_patch_3D.UpdateAttributes();


            //copy UV per Wedge into per vert
            vcg::tri::UV_Utils<TriMeshType>::CopyWedgeVertUV(curr_patch_3D);


            //then do the 2D version
            vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(curr_patch_UV,curr_patch_3D);
            //copy UV into coords
            for (size_t j=0;j<curr_patch_UV.vert.size();j++)
            {
                vcg::Point2<ScalarType> UV=curr_patch_UV.vert[j].T().P();
                //std::cout<<"UV:"<<UV.X()<<","<<UV.Y()<<std::endl;
                curr_patch_UV.vert[j].P()=CoordType(UV.X(),UV.Y(),0);
            }
            curr_patch_UV.UpdateAttributes();

            vcg::tri::io::ExporterOBJ<TriMeshType>::Save(curr_patch_3D,Name_3D_mesh.c_str(),vcg::tri::io::Mask::IOM_ALL);
            vcg::tri::io::ExporterOBJ<TriMeshType>::Save(curr_patch_UV,Name_UV_mesh.c_str(),vcg::tri::io::Mask::IOM_ALL);
        }
    }

    void RemoveOnSymmetryPathIfPossible()
    {
        //select along boders, so it is kept as new border when merged
        //std::set<std::pair<CoordType,CoordType> > BorderE;
        for (size_t i=0;i<deformed_mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!vcg::face::IsBorder(deformed_mesh.face[i],j))continue;
                deformed_mesh.face[i].SetFaceEdgeS(j);
            }

        //merge mesh
        vcg::tri::Clean<TraceMesh>::RemoveDuplicateVertex(deformed_mesh);
        vcg::tri::Allocator<TraceMesh>::CompactEveryVector(deformed_mesh);
        deformed_mesh.UpdateAttributes();

        //update the field
        typedef PatchTracer<TriMeshType,MeshArapQuality<TriMeshType> > PTracerType;
        deformed_mesh.UpdateAttributes();
        vcg::tri::CrossField<TriMeshType>::UpdateSingularByCross(deformed_mesh);
        vcg::tri::CrossField<TriMeshType>::SetVertCrossVectorFromFace(deformed_mesh);

        //reinit the graph
        VertexFieldGraph<TriMeshType> VGraph(deformed_mesh);
        VGraph.InitGraph(false);

        //save the selected ones
        std::set<std::pair<CoordType,CoordType> > EdgeS;
        for (size_t i=0;i<deformed_mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!deformed_mesh.face[i].IsFaceEdgeS(j))continue;

                CoordType P0=deformed_mesh.face[i].P0(j);
                CoordType P1=deformed_mesh.face[i].P1(j);

                std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
                EdgeS.insert(Key);
            }

        //        //initialize the tracer
        //        PTracerType PTr(VGraph);
        //        PTr.InitTracer(100,false);

        //initialize the tracer
        PTracerType PTr(VGraph);

        PTr.CClarkability=-1;
        PTr.match_valence=match_valence;

        if (allow_self_glue || use_darts)
        {
            PTr.MinVal=0;
            PTr.CheckQuadrangulationLimits=false;
        }
        PTr.MaxVal=max_corners;
        //PTr.Concave_Need=1;
        PTr.AllowDarts=use_darts;
        PTr.AllowSelfGluedPatch=allow_self_glue;


        if (max_compression<max_tension)
        {
            MeshArapQuality<TriMeshType>::MaxQ()=max_tension;
            MeshArapQuality<TriMeshType>::MinQ()=max_compression;
            MeshArapQuality<TriMeshType>::RemeshOnTest()=remesh_on_test;
            PTr.check_quality_functor=check_stress;
        }

        PTr.InitTracer(100,false);

        //then restore the selected
        for (size_t i=0;i<deformed_mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                CoordType P0=deformed_mesh.face[i].P0(j);
                CoordType P1=deformed_mesh.face[i].P1(j);
                std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
                if (EdgeS.count(Key)==0)continue;
                deformed_mesh.face[i].SetFaceEdgeS(j);
            }

        //reinitialize path from selected
        PTr.ReinitPathFromEdgeSel();

        if (use_darts)
        {
            PTr.SplitIntoSubPaths();
            PTr.SplitIntoIntervals(PTr.ChoosenPaths);
        }

        //then set only removeable the ones that are on symmetry line
        PTr.SetAllUnRemoveable();
        for (size_t i=0;i<PTr.ChoosenPaths.size();i++)
            for (size_t j=0;j<PTr.ChoosenPaths[i].PathNodes.size()-1;j++)
            {
                size_t N0=PTr.ChoosenPaths[i].PathNodes[j];
                size_t N1=PTr.ChoosenPaths[i].PathNodes[j+1];
                CoordType P0=VGraph.NodePos(N0);
                CoordType P1=VGraph.NodePos(N1);
                CoordType AvgP=(P0+P1)/2;
                CoordType Proj=Symmetrizer<TraceMesh>::SymmPlane().Projection(AvgP);
                if ((Proj-AvgP).Norm()<0.00001)
                    PTr.ChoosenPaths[i].Unremovable=false;
            }
        //then do the actual remove
        PTr.RemovePaths();
        //update the selection on the mesh
        SelectMeshPatchBorders(VGraph,PTr.ChoosenPaths);
    }

    void GetVertPriorityByConvexity(std::vector<ScalarType> &Values)
    {
        Values.clear();

        Eigen::MatrixXi F;
        typename vcg::tri::MeshToMatrix<TriMeshType>::MatrixXm Vf;

        Eigen::MatrixXd PD1,PD2,PV1,PV2;
        vcg::tri::MeshToMatrix<TriMeshType>::GetTriMeshData(half_def_mesh,F,Vf);
        Eigen::MatrixXd V = Vf.template cast<double>();

        igl::principal_curvature(V,F,PD1,PD2,PV1,PV2,4,true);

        //then compute convexity
        for (size_t i=0;i<half_def_mesh.vert.size();i++)
        {
            ScalarType V0=PV1(i,0);
            ScalarType V1=PV2(i,0);

            ScalarType Curv=V0*V1;

            Values.push_back(Curv);
        }
    }

    void ColorByConvexity()
    {
        std::vector<ScalarType> Values;
        GetVertPriorityByConvexity(Values);

        //then compute convexity
        std::vector<std::pair<ScalarType,size_t> > val;

        for (size_t i=0;i<half_def_mesh.vert.size();i++)
        {
            val.push_back(std::pair<ScalarType,size_t>(Values[i],i));
        }
        std::sort(val.begin(),val.end());
        //        size_t minI=val.size()*0.2;
        //        size_t maxI=val.size()*0.8;
        //        ScalarType minV=val[minI];
        //        ScalarType maxV=val[maxI];
        //        for (size_t i=0;i<deformed_mesh.vert.size();i++)
        //        {
        //            std::max(deformed_mesh.vert[i].Q(),minV);
        //            std::min(deformed_mesh.vert[i].Q(),maxV);
        //        }
        for (size_t i=0;i<val.size();i++)
        {
            half_def_mesh.vert[val[i].second].Q()=i;
        }
        //vcg::tri::UpdateQuality<TriMeshType>::Ve
        vcg::tri::UpdateColor<TriMeshType>::PerVertexQualityRamp(half_def_mesh);
        vcg::tri::io::ExporterPLY<TriMeshType>::Save(half_def_mesh,"TestCol.ply",vcg::tri::io::Mask::IOM_VERTCOLOR);
    }

    void BatchProcess(const std::vector<std::vector<CoordType> > &PickedPoints,
                      //const std::vector<bool > &Soft,
                      bool writeDebug=false,
                      bool writeTime=true)
    {
        RestoreInitMesh();
        size_t t0=clock();

        MakeMeshSymmetric(PickedPoints,false);

        //        vcg::tri::io::ExporterPLY<TriMeshType>::Save(half_def_mesh,"dede0ply");

        size_t t1=clock();
        ComputeField(false);
        size_t t2=clock();
        //        //TEST, REMOVE CONSTRAINT
        //        vcg::tri::io::ExporterPLY<TriMeshType>::Save(half_def_mesh,"dede1.ply");
        //        //TEST, REMOVE CONSTRAINTS

        TracePatch(false,writeDebug);

        if (remove_along_symmetry)
            RemoveOnSymmetryPathIfPossible();

        size_t t3=clock();
        DoParametrize();
        size_t t4=clock();
        if (writeTime)
        {
            std::cout<<"Time Symmetrize:"<<(t1-t0)/(ScalarType)CLOCKS_PER_SEC<<std::endl;
            std::cout<<"Time Field Computation:"<<(t2-t1)/(ScalarType)CLOCKS_PER_SEC<<std::endl;
            std::cout<<"Time Patch Tracing:"<<(t3-t2)/(ScalarType)CLOCKS_PER_SEC<<std::endl;
            std::cout<<"Time Parametrize:"<<(t4-t3)/(ScalarType)CLOCKS_PER_SEC<<std::endl;
        }
    }

    ~Parafashion()
    {
        for (size_t i=0;i<SubMeshes.size();i++)
            delete(SubMeshes[i]);
    }

    Parafashion(TriMeshType &_deformed_mesh,
                TriMeshType &_reference_mesh,
                AnimationManager<TriMeshType> &_AManag):
        deformed_mesh(_deformed_mesh),
        reference_mesh(_reference_mesh),
        AManag(_AManag)
    {
        PMode=PMAvgTJuncions;
        FMode=FMCurvature;
        PrioMode=PrioModBlend;
        match_valence=false;
        check_stress=true;
        param_boundary=0.03;
        max_corners=6;
        max_compression=-0.1;
        max_tension=0.1;
        continuity_seams=false;
        continuity_darts =false;
        use_darts=true;
        allow_self_glue=true;
        remove_along_symmetry=false;
        remesh_on_test=false;
        UVMode=PMCloth;
        CheckUVIntersection=true;
        SmoothBeforeRemove=true;
    }
};

