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
                             bool DebugMsg=false);

};

#endif
