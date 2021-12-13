#ifndef FIELD_COMPUTATION
#define FIELD_COMPUTATION

#include <vcg/complex/complex.h>

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
