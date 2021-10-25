#ifndef ANIMATION_MANAGER
#define ANIMATION_MANAGER

#include <vector>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/algorithms/closest.h>
#include <mesh_type.h>
#include <wrap/igl/smooth_field.h>

#define ANISOTR_PERCENTILE 0.1

//struct Jacobian{

//    // tangent directions (columns of J)
//    vec3 u;
//    vec3 v;

//    void setForTriangle( vec3 p0, vec3 p1, vec3 p2,
//                         vec2 u0, vec2 u1, vec2 u2 ){
//        p1-=p0;
//        p2-=p0;
//        u1-=u0;
//        u2-=u0;
//        scalar det = cross(u1,u2);

//        u1/=det;
//        u2/=det;

//        u = p1*gety(u2) - p2*gety(u1);
//        v = p2*getx(u1) - p1*getx(u2) ;
//    }

//    void hammer(scalar degrees){
//        enforceUnitary();
//        boundAngle(degrees);
//    }

//    void print() const{
//        printf("U,V=\n");
//        ::print(u);
//        ::print(v);
//    }

//    void drawGL( vec3 center , scalar size ) const; /* implemented in drawgl.cpp */

//    scalar areaMult() const{
//        return length(cross(u,v));
//    }

//    vec3 n() const{
//        return normalize(cross(u,v));
//    }

//    /* linear operations (for mixing) */
//    void operator += (const Jacobian &other){ u += other.u; v += other.v; }
//    void operator *= (scalar k){ u *= k; v *= k; }
//    Jacobian operator * (scalar k) const { Jacobian res = *this; res*=k; return res; }
//    Jacobian operator + (const Jacobian &other) const { Jacobian res = *this; res += other; return res; }

//private:

//    void enforceUnitary(){
//        u = normalize(u);
//        v = normalize(v);
//    }

//    void boundAngle(scalar degrees){
//        const scalar ANGLE_MIN = 90 - degrees;
//        const scalar ANGLE_MAX = 90 + degrees;

//        const scalar COS_MAX = cos( toRad( ANGLE_MIN/2.0 )  );
//        const scalar SIN_MAX = sin( toRad( ANGLE_MIN/2.0 )  );
//        const scalar COS_MIN = cos( toRad( ANGLE_MAX/2.0 )  );
//        const scalar SIN_MIN = sin( toRad( ANGLE_MAX/2.0 )  );

//        vec3 bisect = normalize(u+v);
//        scalar cos = dot( bisect , v );

//        if ((cos<COS_MIN) || (cos>COS_MAX)){
//            scalar newcos = (cos>COS_MAX) ? COS_MAX : COS_MIN ;
//            scalar newsin = (cos>COS_MAX) ? SIN_MAX : SIN_MIN ;

//            vec3 bisectT = normalize( cross( bisect, cross(u,v) ) );

//            u = bisect*newcos + bisectT*newsin;
//            v = bisect*newcos - bisectT*newsin;
//        }
//    }

// public:

//    void perturb(scalar howmuch){
//        vec3 nn = n();
//        u+= randomMaxUnitVec()*howmuch;
//        v+= randomMaxUnitVec()*howmuch;
//        // put back into plane
//        u = cross( cross(nn,u),nn );
//        v = cross( cross(nn,v),nn );
//    }


//};

template <class FaceType>
class Jacobian{
    typedef typename FaceType::VertexType VertexType;
    typedef typename FaceType::ScalarType ScalarType;
    typedef typename FaceType::CoordType CoordType;
    typedef typename vcg::Point2<ScalarType> Point2X;

    // tangent directions (columns of J)

    static void FromUV( CoordType p0,CoordType p1,CoordType p2,
                        Point2X u0,Point2X u1,Point2X u2,
                        CoordType &u,CoordType &v)
    {
        p1-=p0;
        p2-=p0;
        u1-=u0;
        u2-=u0;

        ScalarType det = u1 ^ u2;

        u1/=det;
        u2/=det;

        u = p1*u2.Y() - p2*u1.Y();
        v = p2*u1.X() - p1*u2.X();
    }


    static void From3DTris(CoordType p0,CoordType p1,CoordType p2,
                           CoordType p0_def,CoordType p1_def,CoordType p2_def,
                           CoordType &u,CoordType &v)
    {
        //retrieve UV Coords
        CoordType bary=(p0_def + p1_def + p2_def)/3;
        p0_def-=bary;
        p1_def-=bary;
        p2_def-=bary;

        CoordType N_source=(p1_def-p0_def)^(p2_def-p0_def);
        CoordType N_dest(0,0,1);
        N_source.Normalize();
        N_dest.Normalize();

        vcg::Matrix33<ScalarType> RotM=vcg::RotationMatrix(N_source,N_dest);
        p0_def=RotM*p0_def;
        p1_def=RotM*p1_def;
        p2_def=RotM*p2_def;

        assert(fabs(p0_def.Z())<0.0001);
        assert(fabs(p1_def.Z())<0.0001);
        assert(fabs(p2_def.Z())<0.0001);

        Point2X UV0(p0_def.X(),p0_def.Y());
        Point2X UV1(p1_def.X(),p1_def.Y());
        Point2X UV2(p2_def.X(),p2_def.Y());

        FromUV(p0,p1,p2,UV0,UV1,UV2,u,v);
    }

public:

    static void FromFace(const FaceType &F,CoordType &u,CoordType &v)
    {
        const VertexType *v0=F.cV(0);
        const VertexType *v1=F.cV(1);
        const VertexType *v2=F.cV(2);

        CoordType P0=v0->P();
        CoordType P1=v1->P();
        CoordType P2=v2->P();

        CoordType RPos0=v0->RPos;
        CoordType RPos1=v1->RPos;
        CoordType RPos2=v2->RPos;

        From3DTris(P0,P1,P2,RPos0,RPos1,RPos2,u,v);
    }
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

//    //reprojecting onto the copy mesh
    std::vector<size_t> VertFaceIdx;
    std::vector<CoordType> VertFaceBary;

    //reprojecting onto the copy mesh
    std::vector<size_t> FaceFaceIdx;

    //copy of the mesh to reproject on
    TriMeshType animated_template_shape;

    //per frame curvature
    std::vector<std::vector<CoordType> > PerFrameCurvVect;

    //per frame curvature anisotropy
    std::vector<std::vector<ScalarType> > PerFrameCurvAnis;

    //per frame Jacobian first and Second Direction
    std::vector<std::vector<CoordType> > JU,JV;

//    void InterpolateVertField(const size_t &IndexV,
//                              const size_t &IndexFrame,
//                              CoordType &CurvDirection,
//                              ScalarType &AnisotropyVal,
//                              CoordType &Norm)
//    {
//        assert(IndexV<target_shape.vert.size());
//        assert(IndexFrame<PerFramePos.size());

//        size_t IndexFace=VertFaceIdx[IndexV];
//        assert(IndexFace<animated_template_shape.face.size());
//        assert(IndexFace<PerFrameCurvVect[IndexFrame].size());
//        Norm=animated_template_shape.face[IndexFace].N();
//        CurvDirection=PerFrameCurvVect[IndexFrame][IndexFace];
//        AnisotropyVal=PerFrameCurvAnis[IndexFrame][IndexFace];
//    }

//    void InterpolateVertJ(const size_t &IndexV,const size_t &IndexFrame,
//                          CoordType &VertJU,CoordType &VertJV,CoordType &Norm)
//    {
//        assert(IndexV<target_shape.vert.size());
//        assert(IndexFrame<PerFramePos.size());

//        size_t IndexFace=VertFaceIdx[IndexV];
//        assert(IndexFace<animated_template_shape.face.size());
//        assert(IndexFace<JU[IndexFrame].size());
//        assert(IndexFace<JV[IndexFrame].size());

//        Norm=animated_template_shape.face[IndexFace].N();
//        VertJU=JU[IndexFrame][IndexFace];
//        VertJV=JV[IndexFrame][IndexFace];
//    }


//    void InterpolateFaceField(const size_t &IndexFace,const size_t &IndexFrame,
//                              CoordType &InterpCurvDirection,
//                              ScalarType &InterpAnisotropyVal)
//    {
//        assert(IndexFace<target_shape.face.size());

//        VertexType *v0=target_shape.face[IndexFace].V(0);
//        VertexType *v1=target_shape.face[IndexFace].V(1);
//        VertexType *v2=target_shape.face[IndexFace].V(2);

//        size_t IndexV0=vcg::tri::Index(target_shape,v0);
//        size_t IndexV1=vcg::tri::Index(target_shape,v1);
//        size_t IndexV2=vcg::tri::Index(target_shape,v2);

//        std::vector<CoordType> CurvDirection;CurvDirection.resize(3);
//        std::vector<ScalarType> AnisotropyVal;AnisotropyVal.resize(3);
//        std::vector<CoordType> Norm;Norm.resize(3);

//        //Normal of each face
//        InterpolateVertField(IndexV0,IndexFrame,CurvDirection[0],AnisotropyVal[0],Norm[0]);
//        InterpolateVertField(IndexV1,IndexFrame,CurvDirection[1],AnisotropyVal[1],Norm[1]);
//        InterpolateVertField(IndexV2,IndexFrame,CurvDirection[2],AnisotropyVal[2],Norm[2]);

//        InterpAnisotropyVal=(AnisotropyVal[0]+AnisotropyVal[1]+AnisotropyVal[2])/3;
//        CoordType TargetN=target_shape.face[IndexFace].N();
//        InterpCurvDirection=vcg::tri::InterpolateNRosy3D(CurvDirection,Norm,AnisotropyVal,4,TargetN);
//    }

    void InterpolateFaceField(const size_t &IndexFace,const size_t &IndexFrame,
                              CoordType &InterpCurvDirection,
                              ScalarType &InterpAnisotropyVal)
    {
        assert(IndexFace<target_shape.face.size());

        assert(IndexFrame<NumFrames());
        assert(IndexFrame<JU.size());
        assert(IndexFrame<JV.size());

        assert(IndexFace<target_shape.face.size());

        assert(IndexFace<FaceFaceIdx.size());
        size_t IndexFAnim=FaceFaceIdx[IndexFace];


        InterpCurvDirection=PerFrameCurvVect[IndexFrame][IndexFAnim];
        InterpAnisotropyVal=PerFrameCurvAnis[IndexFrame][IndexFAnim];
//        FaceJV=JV[IndexFrame][IndexFAnim];

//        VertexType *v0=target_shape.face[IndexFace].V(0);
//        VertexType *v1=target_shape.face[IndexFace].V(1);
//        VertexType *v2=target_shape.face[IndexFace].V(2);

//        size_t IndexV0=vcg::tri::Index(target_shape,v0);
//        size_t IndexV1=vcg::tri::Index(target_shape,v1);
//        size_t IndexV2=vcg::tri::Index(target_shape,v2);

//        std::vector<CoordType> CurvDirection;CurvDirection.resize(3);
//        std::vector<ScalarType> AnisotropyVal;AnisotropyVal.resize(3);
//        std::vector<CoordType> Norm;Norm.resize(3);

//        //Normal of each face
//        InterpolateVertField(IndexV0,IndexFrame,CurvDirection[0],AnisotropyVal[0],Norm[0]);
//        InterpolateVertField(IndexV1,IndexFrame,CurvDirection[1],AnisotropyVal[1],Norm[1]);
//        InterpolateVertField(IndexV2,IndexFrame,CurvDirection[2],AnisotropyVal[2],Norm[2]);

//        InterpAnisotropyVal=(AnisotropyVal[0]+AnisotropyVal[1]+AnisotropyVal[2])/3;
//        CoordType TargetN=target_shape.face[IndexFace].N();
//        InterpCurvDirection=vcg::tri::InterpolateNRosy3D(CurvDirection,Norm,AnisotropyVal,4,TargetN);
    }

    void InterpolateFaceStretch(const size_t &IndexFace,const size_t &IndexFrame,
                                CoordType &FaceJU,CoordType &FaceJV)
    {
        assert(IndexFrame<NumFrames());
        assert(IndexFrame<JU.size());
        assert(IndexFrame<JV.size());

        assert(IndexFace<target_shape.face.size());

        assert(IndexFace<FaceFaceIdx.size());
        size_t IndexFAnim=FaceFaceIdx[IndexFace];

        FaceJU=JU[IndexFrame][IndexFAnim];
        FaceJV=JV[IndexFrame][IndexFAnim];
     }


    CoordType InterpolatePos(size_t IndexV,size_t IndexFrame)
    {
        assert(IndexV<target_shape.vert.size());
        assert(IndexFrame<PerFramePos.size());

        size_t IndexF=VertFaceIdx[IndexV];
        assert(IndexF<animated_template_shape.face.size());

        CoordType Bary=VertFaceBary[IndexV];
        VertexType *V0=animated_template_shape.face[IndexF].V(0);
        VertexType *V1=animated_template_shape.face[IndexF].V(1);
        VertexType *V2=animated_template_shape.face[IndexF].V(2);

        size_t IndexV0=vcg::tri::Index(animated_template_shape,V0);
        size_t IndexV1=vcg::tri::Index(animated_template_shape,V1);
        size_t IndexV2=vcg::tri::Index(animated_template_shape,V2);

        CoordType P0=PerFramePos[IndexFrame][IndexV0];
        CoordType P1=PerFramePos[IndexFrame][IndexV1];
        CoordType P2=PerFramePos[IndexFrame][IndexV2];

        return (P0*Bary.X()+P1*Bary.Y()+P2*Bary.Z());
    }

    ScalarType percentileAnis;
    ScalarType percentileStretchCompress;
    //ScalarType percentileStretch;

    void UpdateTemplateToFrame(size_t IndexFrame)
    {
        assert(IndexFrame<PerFramePos.size());
        assert(PerFramePos[IndexFrame].size()==animated_template_shape.vert.size());

        for (size_t i=0;i<animated_template_shape.vert.size();i++)
            animated_template_shape.vert[i].P()=PerFramePos[IndexFrame][i];

        animated_template_shape.UpdateAttributes();
    }

    void InitPerFrameCurvature()
    {
        if (NumFrames()==0)return;
        PerFrameCurvVect.clear();
        PerFrameCurvAnis.clear();

        PerFrameCurvVect.resize(NumFrames());
        PerFrameCurvAnis.resize(NumFrames());

        std::vector<ScalarType> AnisValue;

        //animated_template_shape.InitRPos();

        for (size_t i=0;i<NumFrames();i++)
        {
            UpdateTemplateToFrame(i);
            vcg::tri::FieldSmoother<TriMeshType>::InitByCurvature(animated_template_shape,4);
            for (size_t j=0;j<animated_template_shape.face.size();j++)
            {
                PerFrameCurvVect[i].push_back(animated_template_shape.face[j].PD1());
                PerFrameCurvAnis[i].push_back(animated_template_shape.face[j].Q());
                AnisValue.push_back(animated_template_shape.face[j].Q());
            }
        }
        animated_template_shape.RestoreRPos();
        animated_template_shape.UpdateAttributes();

        //then get the percentile of anisotropy
        std::sort(AnisValue.begin(),AnisValue.end());
        size_t Index_Percentile=(AnisValue.size()*(1-ANISOTR_PERCENTILE));
//        std::cout<<"Index:"<<Index_Percentile<<std::endl;
//        std::cout<<"Size:"<<AnisValue.size()<<std::endl;

        percentileAnis=AnisValue[Index_Percentile];
 //       std::cout<<"Value:"<<percentileAnis<<std::endl;
    }


    ScalarType getKForStretchCompression(CoordType Vect)
    {
        ScalarType K=Vect.Norm();
        if (K<1){
            K=-1/K;
            K+=1;
            assert(K<=0);
        }
        else
        {
            K-=1;
            assert(K>=0);
        }
        return K;
    }

    void InitPerFrameJacobian()
    {
        if (NumFrames()==0)return;
        JU.clear();
        JV.clear();

        JU.resize(NumFrames());
        JV.resize(NumFrames());

        std::vector<ScalarType> JValue;

        for (size_t i=0;i<NumFrames();i++)
        {
            UpdateTemplateToFrame(i);
            for (size_t j=0;j<animated_template_shape.face.size();j++)
            {
                CoordType U,V;
                Jacobian<FaceType>::FromFace(animated_template_shape.face[j],U,V);
                JU[i].push_back(U);
                JV[i].push_back(V);


                ScalarType KValU=getKForStretchCompression(U);
                ScalarType KValV=getKForStretchCompression(V);
                JValue.push_back(fabs(KValU));
                JValue.push_back(fabs(KValV));
            }
        }
        animated_template_shape.RestoreRPos();
        animated_template_shape.UpdateAttributes();

        std::sort(JValue.begin(),JValue.end());


        size_t Index_Percentile_Stretch_Compress=(JValue.size()*(1-ANISOTR_PERCENTILE));

//        std::cout<<"Index:"<<Index_Percentile<<std::endl;
//        std::cout<<"Size:"<<AnisValue.size()<<std::endl;

        percentileStretchCompress=JValue[Index_Percentile_Stretch_Compress];
    }

    void ClampAnisotropyByPercentile()
    {
        for (size_t i=0;i<target_shape.face.size();i++)
        {
            target_shape.face[i].Q()=std::min(target_shape.face[i].Q(),percentileAnis);

            ScalarType Norm1=target_shape.face[i].PD1().Norm();
            Norm1=std::min(Norm1,percentileAnis);

            ScalarType Norm2=target_shape.face[i].PD2().Norm();
            Norm2=std::min(Norm2,percentileAnis);

            target_shape.face[i].PD1().Normalize();
            target_shape.face[i].PD2().Normalize();
            target_shape.face[i].PD1()*=Norm1;
            target_shape.face[i].PD2()*=Norm2;
        }
    }

    void ClampStretchByPercentile()
    {
        for (size_t i=0;i<target_shape.face.size();i++)
        {
            //then check maximum between compression and tension
            target_shape.face[i].K1()=std::max(target_shape.face[i].K1(),-percentileStretchCompress);
            target_shape.face[i].K1()=std::min(target_shape.face[i].K1(),percentileStretchCompress);
            target_shape.face[i].K2()=std::max(target_shape.face[i].K2(),-percentileStretchCompress);
            target_shape.face[i].K2()=std::min(target_shape.face[i].K2(),percentileStretchCompress);
            target_shape.face[i].Q()=std::max(fabs(target_shape.face[i].K1()),
                                              fabs(target_shape.face[i].K2()));
        }
    }

    void UpdateFaceCurvatureField(size_t IndexFrame)
    {
        assert(IndexFrame<NumFrames());
        for (size_t i=0;i<target_shape.face.size();i++)
        {
            CoordType InterpCurvDirection;
            ScalarType InterpAnisotropyVal;
            InterpolateFaceField(i,IndexFrame,InterpCurvDirection,InterpAnisotropyVal);
            target_shape.face[i].PD1()=InterpCurvDirection;
            target_shape.face[i].PD2()=target_shape.face[i].PD1()^target_shape.face[i].N();
            target_shape.face[i].PD1().Normalize();
            target_shape.face[i].PD2().Normalize();
            target_shape.face[i].PD1()*=InterpAnisotropyVal;
            target_shape.face[i].PD2()*=InterpAnisotropyVal;
            //            target_shape.face[i].K1()=InterpAnisotropyVal;
            //            target_shape.face[i].K2()=InterpAnisotropyVal;
            target_shape.face[i].Q()=InterpAnisotropyVal;
        }
    }


    void UpdateFaceStretchField(size_t IndexFrame)
    {
        assert(IndexFrame<NumFrames());

//        ScalarType minV=std::numeric_limits<ScalarType>::max();
//        ScalarType maxV=std::numeric_limits<ScalarType>::min();

        for (size_t i=0;i<target_shape.face.size();i++)
        {
            CoordType FaceJU,FaceJV;
            InterpolateFaceStretch(i,IndexFrame,FaceJU,FaceJV);
//            std::cout<<"Norm V:"<<FaceJV.Norm()<<std::endl;
//            std::cout<<"Norm U:"<<FaceJU.Norm()<<std::endl;
            target_shape.face[i].PD1()=FaceJU;
            target_shape.face[i].PD2()=FaceJV;

            target_shape.face[i].K1()=getKForStretchCompression(FaceJU);
            target_shape.face[i].K2()=getKForStretchCompression(FaceJV);

//            minV=std::min(minV,target_shape.face[i].K1());
//            minV=std::min(minV,target_shape.face[i].K2());
//            maxV=std::max(maxV,target_shape.face[i].K1());
//            maxV=std::max(maxV,target_shape.face[i].K2());
        }
//        std::cout<<"Max:"<<maxV<<std::endl;
//        std::cout<<"Min:"<<minV<<std::endl;
    }

public:

    ScalarType MaxAnisotropy()
    {
        return percentileAnis;
    }

    ScalarType MaxStretchCompress()
    {
        return percentileStretchCompress;
    }

//    ScalarType MaxAnisotropy()
//    {
//        return percentileAnis;
//    }

//    ScalarType MaxAnisotropy()
//    {
//        return percentileAnis;
//    }

    void ColorByAnisotropy()
    {
        //clamp to the max
        ClampAnisotropyByPercentile();
        //one smooth step
        vcg::tri::UpdateQuality<TriMeshType>::VertexFromFace(target_shape);
        vcg::tri::UpdateQuality<TriMeshType>::FaceFromVertex(target_shape);
        //then update the color
        vcg::tri::UpdateColor<TriMeshType>::PerFaceQualityGray(target_shape,percentileAnis,0);
    }

    void ColorByStretch()
    {
        //clamp to the max
        ClampStretchByPercentile();
        //one smooth step
        vcg::tri::UpdateQuality<TriMeshType>::VertexFromFace(target_shape);
        vcg::tri::UpdateQuality<TriMeshType>::FaceFromVertex(target_shape);
        //then update the color
        vcg::tri::UpdateColor<TriMeshType>::PerFaceQualityGray(target_shape,percentileStretchCompress,0);
    }

    bool LoadPosFrames(const std::string &path)
    {
        FILE *f=fopen(path.c_str(),"rt");
        if (f==NULL)return false;

        int NumFrames;
        fscanf(f,"%d\n",&NumFrames);
        std::cout<<"There are Frames:"<<NumFrames<<std::endl;

        PerFramePos.clear();
        PerFramePos.resize(NumFrames);

        for (size_t i=0;i<PerFramePos.size();i++)
        {
            for (size_t j=0;j<animated_template_shape.vert.size();j++)
            {
                float XPos,YPos,ZPos;
                fscanf(f,"%f,%f,%f\n",&XPos,&YPos,&ZPos);
                PerFramePos[i].push_back(CoordType((ScalarType)XPos,(ScalarType)YPos,(ScalarType)ZPos));
            }
        }
        fclose(f);

        std::cout<<"Initializing Per Frame Curvature"<<std::endl;
        InitPerFrameCurvature();
        std::cout<<"Done"<<std::endl;
        std::cout<<"Initializing Per Frame Jacobian"<<std::endl;
        InitPerFrameJacobian();
        std::cout<<"Done"<<std::endl;
        return true;
    }

    void UpdateProjectionBasis()
    {
        vcg::GridStaticPtr<FaceType,ScalarType> Grid;
        Grid.Set(animated_template_shape.face.begin(),animated_template_shape.face.end());

        ScalarType MaxD=animated_template_shape.bbox.Diag();

        VertFaceIdx.clear();
        VertFaceBary.clear();
        for (size_t i=0;i<target_shape.vert.size();i++)
        {
            ScalarType MinD;
            CoordType closestPt,normI,baryP;
            FaceType *f=vcg::tri::GetClosestFaceBase(animated_template_shape,Grid,target_shape.vert[i].P(),
                                                     MaxD,MinD,closestPt,normI,baryP);
            assert(f!=NULL);
            size_t IndexF=vcg::tri::Index(animated_template_shape,f);

            VertFaceIdx.push_back(IndexF);
            VertFaceBary.push_back(baryP);
        }

        FaceFaceIdx.clear();
        for (size_t i=0;i<target_shape.face.size();i++)
        {
            CoordType BaryF=(target_shape.face[i].P(0)+
                            target_shape.face[i].P(1)+
                            target_shape.face[i].P(2))/3;
            ScalarType MinD;
            CoordType closestPt,normI,baryP;
            FaceType *f=vcg::tri::GetClosestFaceBase(animated_template_shape,Grid,BaryF,
                                                     MaxD,MinD,closestPt,normI,baryP);
            assert(f!=NULL);
            size_t IndexF=vcg::tri::Index(animated_template_shape,f);

            FaceFaceIdx.push_back(IndexF);
        }
    }

    //    void UpdateCurvatureAndStretch()
    //    {
    //        ComputePerFrameCurvature();
    //    }

    void Init()
    {
        animated_template_shape.Clear();
        vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(animated_template_shape,target_shape);
        animated_template_shape.InitRPos();
        animated_template_shape.UpdateAttributes();

        VertFaceIdx.clear();
        VertFaceBary.clear();
        PerFramePos.clear();
        UpdateProjectionBasis();
    }

    size_t NumFrames()const
    {return PerFramePos.size();}

    void UpdateToFrame(size_t IndexFrame,
                       bool UpdateCurvature=false,
                       bool UpdateStretch=false)
    {
        assert(IndexFrame<PerFramePos.size());
        assert(PerFramePos[IndexFrame].size()==animated_template_shape.vert.size());

        //target_shape.RestoreRPos();
        CoordType CenterTemplate=animated_template_shape.bbox.Center();

        for (size_t i=0;i<target_shape.vert.size();i++)
            target_shape.vert[i].P()=InterpolatePos(i,IndexFrame);

        vcg::tri::UpdateBounding<TriMeshType>::Box(target_shape);
        CoordType CenterAnim=target_shape.bbox.Center();
        for (size_t i=0;i<target_shape.vert.size();i++)
            target_shape.vert[i].P()+=(CenterTemplate-CenterAnim);

        target_shape.UpdateAttributes();

        if (UpdateCurvature)
            UpdateFaceCurvatureField(IndexFrame);
        if (UpdateStretch)
            UpdateFaceStretchField(IndexFrame);
    }


    AnimationManager(TriMeshType &_target_shape):target_shape(_target_shape)
    {}
};

#endif
