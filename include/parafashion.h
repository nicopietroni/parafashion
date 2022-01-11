#pragma once

#include "symmetrizer.h"
#include "parametrizer.h"
#include "field_computation.h"
#include "trace_path_ui.h"
#include <tracing/tracer_interface.h>
#include <tracing/mesh_type.h>
#include "animation_manager.h"

#define PRINT_PARAFASHION_TIMING 
#ifdef PRINT_PARAFASHION_TIMING
#include <chrono>
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;
#endif

enum PatchMode{PMMinTJuncions,PMAvgTJuncions,PMAllTJuncions};

template <class TriMeshType>
class Parafashion
{
    typedef typename TriMeshType::FaceType FaceType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;

    TriMeshType &deformed_mesh;
    TriMeshType &reference_mesh;

    TriMeshType half_def_mesh;

    TriMeshType deformed_mesh_step0;
    TriMeshType reference_mesh_step0;

    TriMeshType deformed_mesh_step1;
    TriMeshType half_def_mesh_step1;

    TriMeshType deformed_mesh_step2;
    TriMeshType half_def_mesh_step2;

public:

    PatchMode PMode;
    bool match_valence;
    bool check_stress;
    bool use_darts;
    bool allow_self_glue;
    bool remove_along_symmetry;

    ScalarType param_boundary;
    size_t max_corners;
    ScalarType max_compression;
    ScalarType max_tension;

    std::vector<typename TraceMesh::CoordType> PatchCornerPos;

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

        MeshArapQuality(){}

        ScalarType operator()(MeshType &m) const
        {
            #ifdef PRINT_PARAFASHION_TIMING
            steady_clock::time_point pre_param = steady_clock::now();
            #endif
            bool success = ClothParametrize<TriMeshType>(m, MaxQ()); // quality-check param, NOT the final one you see on screen


            //vcg::tri::OptimizeUV_ARAP(m,5,0,true);
            //vcg::tri::InitializeArapWithLSCM(m,0);
            
            //evaluate the distortion
            //vcg::tri::Distortion<TraceMesh,false>::SetQasDistorsion(m,vcg::tri::Distortion<TraceMesh,false>::EdgeComprStretch);
            //ScalarType A=0;
            //ScalarType SumA=0;
            /*for (size_t i=0;i<m.face.size();i++)
            {
//                SumA+=vcg::DoubleArea(m.face[i]);
//                A+=m.face[i].Q()*vcg::DoubleArea(m.face[i]);
                if (m.face[i].Q()<(MinQ()))A+=vcg::DoubleArea(m.face[i]);
                if (m.face[i].Q()>MaxQ())A+=vcg::DoubleArea(m.face[i]);
            }*/
            //return (A/SumA);
            #ifdef PRINT_PARAFASHION_TIMING
            steady_clock::time_point post_param = steady_clock::now();
            int param_time = duration_cast<microseconds>(post_param - pre_param).count();
            std::cout << "Param time : " << param_time << " [µs]" << std::endl;
            #endif
            if (success) return 0;
            else return 1000000.0;
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
        half_def_mesh.UpdateAttributes();
        bool PreRemoveStep=true;

        if (max_compression<max_tension)
        {
            MeshArapQuality<TriMeshType>::MaxQ()=max_tension;
            MeshArapQuality<TriMeshType>::MinQ()=max_compression;
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


        if ((!use_darts)&&(!allow_self_glue))
            RecursiveProcess<PTracerType>(PTr,100,true,true,PreRemoveStep,false,false,false,DebugMSG);

//        {
//            std::cout<<"DEDEDE"<<std::endl;
//            RecursiveProcess<PTracerType>(PTr,100,false,false,false,false,false,false,DebugMSG);
//        }

        if ((use_darts)&&(!allow_self_glue))
            RecursiveProcessWithDarts<PTracerType>(PTr,100,true,true,PreRemoveStep,false,false,false,DebugMSG);

        if ((!use_darts)&&(allow_self_glue))
            RecursiveProcessForTexturing<PTracerType>(PTr,100,true,true,PreRemoveStep,false,false,false,DebugMSG);

        if ((use_darts)&&(allow_self_glue))
            RecursiveProcessForTexturingWithDarts<PTracerType>(PTr,100,true,true,PreRemoveStep,false,false,false,DebugMSG);


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

        FieldComputation<TriMeshType>::ComputeField(half_def_mesh,FMBoundary,0.5);//,align_border);
        vcg::tri::CrossField<TriMeshType>::UpdateSingularByCross(half_def_mesh);
        vcg::tri::CrossField<TriMeshType>::SetVertCrossVectorFromFace(half_def_mesh);
        half_def_mesh.InitSingVert();
        half_def_mesh.InitRPos();

        //std::cout<<"3"<<std::endl;
        Symmetrizer<TriMeshType> Symm(deformed_mesh,reference_mesh);

        PreProcessMesh(half_def_mesh,false);
        Symm.CopyFromHalfDefMesh(half_def_mesh);

        //CHECK
        //half_def_mesh.InitRPos();
        //std::cout<<"5"<<std::endl;
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

    void DoParametrize()
    {
        Parametrizer<TriMeshType>::Parametrize(deformed_mesh,param_boundary);
        //parametrized=true;
    }

    void Init()
    {
        deformed_mesh_step0.Clear();
        reference_mesh_step0.Clear();
        vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(deformed_mesh_step0,deformed_mesh);
        vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(reference_mesh_step0,reference_mesh);
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
        PTr.match_valence=false;

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
            PTr.SplitIntoIntervals();

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

    void BatchProcess(const std::vector<std::vector<CoordType> > &PickedPoints,
                      bool writeDebug=false,
                      bool writeTime=true)
    {
        RestoreInitMesh();
        size_t t0=clock();
        MakeMeshSymmetric(PickedPoints,false);
        size_t t1=clock();
        ComputeField(false);
        size_t t2=clock();
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

    Parafashion(TriMeshType &_deformed_mesh,TriMeshType &_reference_mesh):
        deformed_mesh(_deformed_mesh),reference_mesh(_reference_mesh)
    {
        PMode=PMMinTJuncions;
        match_valence=false;
        check_stress=true;
        param_boundary=0.03;
        max_corners=6;
        max_compression=-0.1;
        max_tension=0.1;
        use_darts=true;
        allow_self_glue=true;
        remove_along_symmetry=false;
    }
};

