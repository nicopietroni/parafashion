#ifndef ANIMATION_MANAGER
#define ANIMATION_MANAGER

#include <vector>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/algorithms/closest.h>

#define ANISOTR_PERCENTILE 0.1

template <class FaceType>
class Jacobian{
    typedef typename FaceType::VertexType VertexType;
    typedef typename FaceType::ScalarType ScalarType;
    typedef typename FaceType::CoordType CoordType;
    typedef typename vcg::Point2<ScalarType> Point2X;

    // tangent directions (columns of J)

    static void FromUV( CoordType p0,CoordType p1,CoordType p2,
                        Point2X u0,Point2X u1,Point2X u2,
                        CoordType &u,CoordType &v);


    static void From3DTris(CoordType p0,CoordType p1,CoordType p2,
                           CoordType p0_def,CoordType p1_def,CoordType p2_def,
                           CoordType &u,CoordType &v);

public:

    static void FromFace(const FaceType &F,CoordType &u,CoordType &v);
};


template <class TriMeshType>
class AnimationManager
{
    typedef typename TriMeshType::FaceType FaceType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;

    //per frame per vert position
    std::vector<std::vector<CoordType> > PerFramePos;

    //    //per frame deformation tensor
    //    std::vector<std::vector<CoordType> > PerFrameCurvAnis;

    //the decomposed mesh on which we apply the deformation
    TriMeshType &target_shape;

    //reprojecting onto the copy mesh
    std::vector<size_t> VertFaceIdx;
    std::vector<CoordType> VertFaceBary;

    //reprojecting onto the copy mesh
    std::vector<size_t> FaceFaceIdx;

    //copy of the mesh to reproject on
    TriMeshType animated_template_shape;

    //per frame normal
    std::vector<std::vector<CoordType> > PerFrameNormVect;

    //per frame curvature
    std::vector<std::vector<CoordType> > PerFrameCurvVect;

    //per frame curvature anisotropy
    std::vector<std::vector<ScalarType> > PerFrameCurvAnis;

    //per frame Jacobian first and Second Direction
    std::vector<std::vector<CoordType> > JU,JV;

    //target curvature and weight
    std::vector<CoordType> TargetVect;
    std::vector<ScalarType> TargetAnis;

    void InterpolateFaceField(const size_t &IndexFace,const size_t &IndexFrame,
                              CoordType &InterpCurvDirection,
                              ScalarType &InterpAnisotropyVal);

    void InterpolateFaceStretch(const size_t &IndexFace,const size_t &IndexFrame,
                                CoordType &FaceJU,CoordType &FaceJV);

    void InterpolateFaceNorm(const size_t &IndexFace,const size_t &IndexFrame,
                             CoordType &InterpNorm);

    CoordType InterpolatePos(size_t IndexV,size_t IndexFrame);

    ScalarType percentileAnis;
    ScalarType percentileStretchCompress;
    //ScalarType percentileStretch;

    void UpdateTemplateToFrame(size_t IndexFrame);

    void InitPerFrameCurvature();

    ScalarType getKForStretchCompression(CoordType Vect);

    void InitPerFrameJacobian();

    void ClampAnisotropyByPercentile();

    void ClampStretchByPercentile();

    void UpdateFaceCurvatureField(size_t IndexFrame);

    void UpdateFaceStretchField(size_t IndexFrame);

public:

    ScalarType MaxAnisotropy();

    ScalarType MaxStretchCompress();

    void ColorByAnisotropy();

    void ColorByStretch();

    bool LoadPosFrames(const std::string &path);

    void UpdateProjectionBasis();

    void Init();

    size_t NumFrames()const;

    void UpdateToFrame(size_t IndexFrame,
                       bool UpdateCurvature=false,
                       bool UpdateStretch=false);

    void InitTargetDirectionsAvg();

    void InitTargetDirectionsMax();

    void NormalizeVect(std::vector<ScalarType>  &Values,
                       ScalarType cut_perc=0.1);

    void InitTargetDirectionsOnMesh();

    void UpdateAnimationMesh();

    void TransferDirOnMesh(TriMeshType &target);

    AnimationManager(TriMeshType &_target_shape);
};

#endif
