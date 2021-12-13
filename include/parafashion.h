#ifndef PARAFASHION_H
#define PARAFASHION_H

#include "symmetrizer.h"
#include "parametrizer.h"
#include "field_computation.h"
#include "trace_path_ui.h"
#include <tracer_interface.h>
#include <mesh_type.h>
#include  "animation_manager.h"

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
    ScalarType param_boundary;
    size_t max_corners;
    ScalarType max_compression;
    ScalarType max_tension;

    void CleanMeshAttributes();

    void RestoreInitMesh();

    

    void TracePatch(bool SaveStep=true,
                    bool DebugMSG=false);

    void ComputeField(bool SaveStep=true);

    void AddSharpConstraints(const std::vector<std::vector<CoordType> > &PickedPoints);

    void MakeMeshSymmetric(const std::vector<std::vector<CoordType> > &PickedPoints,
                           bool SaveStep=true);

    void DoParametrize();

    void Init();


    void BatchProcess(const std::vector<std::vector<CoordType> > &PickedPoints,
                      bool writeDebug=false,
                      bool writeTime=true);

    Parafashion(TriMeshType &_deformed_mesh,TriMeshType &_reference_mesh);
};

#endif
