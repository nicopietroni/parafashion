#include <vcg/complex/algorithms/parametrization/uv_utils.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/space/distance2.h>
#include <vcg/complex/algorithms/update/quality.h>
#include "wrap/igl/lscm_parametrization.h"

//enum LenghMode{LModeLocal,LModeGLobal};

//template <class MeshType>
//class FibreElongation
//{
//    typedef typename MeshType::FaceType FaceType;
//    typedef typename MeshType::VertexType VertexType;
//    typedef typename MeshType::CoordType CoordType;
//    typedef typename MeshType::ScalarType ScalarType;

//    typedef typename vcg::Point2<ScalarType> Point2Type;
//    typedef typename vcg::Box2<ScalarType> Box2Type;

//    MeshType &param_mesh;
//    bool WriteDebug;
//    ScalarType MaxScale;

//    struct FiberInt
//    {
//        Point2Type IntUV0,IntUV1;
//        CoordType Int3DP_0,Int3DP_1;
//        size_t FaceIdx;
//        size_t IndexE0,IndexE1;
//        size_t Direction;


//        ScalarType Lenght3D()const
//        {
//            return (Int3DP_0 - Int3DP_1).Norm();
//        }

//        ScalarType LenghtUV()const
//        {
//            return (IntUV0 - IntUV1).Norm();
//        }

//        inline bool operator <(const FiberInt &FInt)const
//        {
//            Point2Type AvgUV0=(IntUV0+IntUV1)/2;
//            Point2Type AvgUV1=(FInt.IntUV0+FInt.IntUV1)/2;
//            if (Direction==0)return (AvgUV0.Y()<AvgUV1.Y());
//            if (Direction==1)return (AvgUV0.X()<AvgUV1.X());
//        }

//        FiberInt(Point2Type &_IntUV0,
//                 Point2Type &_IntUV1,
//                 CoordType &_Int3DP_0,
//                 CoordType &_Int3DP_1,
//                 size_t _FaceIdx,
//                 size_t _IndexE0,
//                 size_t _IndexE1,
//                 size_t _Direction)
//        {
//            IntUV0=_IntUV0;
//            IntUV1=_IntUV1;
//            Int3DP_0=_Int3DP_0;
//            Int3DP_1=_Int3DP_1;
//            FaceIdx=_FaceIdx;
//            Direction=_Direction;
//            IndexE0=_IndexE0;
//            IndexE1=_IndexE1;
//            if (((Direction==0)&&(IntUV0.Y()>IntUV1.Y()))||
//                    ((Direction==1)&&(IntUV0.X()>IntUV1.X())))
//            {
//                std::swap(IntUV0,IntUV1);
//                std::swap(Int3DP_0,Int3DP_1);
//                std::swap(IndexE0,IndexE1);
//            }
//        }

//    };

//    struct Fiber
//    {
//        ScalarType LenParam;
//        ScalarType Len3D;
//        ScalarType RatioL;
//        ScalarType DiffL;
//        std::vector<FiberInt> FiberI;


//        Point2Type CenterUV()
//        {
//           Point2Type P0=FiberI[0].IntUV0;
//           Point2Type P1=FiberI.back().IntUV1;
//           return((P0+P1)/(ScalarType)2);
//        }

//        Fiber()
//        {LenParam=0;Len3D=0;RatioL=0;DiffL=0;}
//    };

//    std::vector<Fiber> UFibers,VFibers;
//    Box2Type UVBox;
//    bool perVert;
//    Point2Type fibre_space;

//    ScalarType MaxElongR;
//    ScalarType MaxCompressR;
//    ScalarType MaxElongL;
//    ScalarType MaxCompressL;

//    ScalarType TotDiffL;

//    void ComputeTargetDispl(std::vector<Fiber> &Fibres,
//                            std::vector<size_t> &IndexVert,
//                            std::vector<Point2Type> &TargetUV,
//                            std::vector<ScalarType> &TargetW)
//    {
//        //std::cout<<"BBox:"<<param_mesh.bbox.Diag()<<std::endl;
//        for (size_t i=0;i<Fibres.size();i++)
//        {
//            //std::cout<<"Size:"<<Fibres[i].FiberI.size()<<std::endl;
//            if (Fibres[i].FiberI.size()==0)continue;

//            ScalarType Len3D=Fibres[i].Len3D;
//            ScalarType LenUV=Fibres[i].LenParam;
//            ScalarType DiffL=Len3D-LenUV;
//            //std::cout<<"DiffL:"<<DiffL<<std::endl;

//            //get vertices and intersection points
//            size_t IndexFMax=Fibres[i].FiberI.back().FaceIdx;
//            size_t IndexEMax=Fibres[i].FiberI.back().IndexE1;
//            size_t IndexFMin=Fibres[i].FiberI[0].FaceIdx;
//            size_t IndexEMin=Fibres[i].FiberI[0].IndexE0;

//            //intersection points
//            Point2Type IntEMax=Fibres[i].FiberI.back().IntUV1;
//            Point2Type IntEMin=Fibres[i].FiberI[0].IntUV0;

//            assert(IndexFMax<param_mesh.face.size());
//            assert(IndexEMax<3);

//            VertexType *V0Max=param_mesh.face[IndexFMax].V0(IndexEMax);
//            VertexType *V1Max=param_mesh.face[IndexFMax].V1(IndexEMax);
//            VertexType *V0Min=param_mesh.face[IndexFMin].V0(IndexEMin);
//            VertexType *V1Min=param_mesh.face[IndexFMin].V1(IndexEMin);

//            ScalarType EMinL=(V0Min->T().P()-V1Min->T().P()).Norm();
//            ScalarType EMinV0L=(V0Min->T().P()-IntEMin).Norm();
//            ScalarType EMinV1L=(V1Min->T().P()-IntEMin).Norm();
//            assert(fabs(EMinV0L+EMinV1L-EMinL)<0.0001);

//            ScalarType EMaxL=(V0Max->T().P()-V1Max->T().P()).Norm();
//            ScalarType EMaxV0L=(V0Max->T().P()-IntEMax).Norm();
//            ScalarType EMaxV1L=(V1Max->T().P()-IntEMax).Norm();
//            assert(fabs(EMaxV0L+EMaxV1L-EMaxL)<0.0001);

//            ScalarType W0Min=1-EMinV0L/EMinL;
//            ScalarType W1Min=1-W0Min;

//            ScalarType W0Max=1-EMaxV0L/EMaxL;
//            ScalarType W1Max=1-W0Max;


//            size_t IndexV0Max=vcg::tri::Index(param_mesh,V0Max);
//            size_t IndexV1Max=vcg::tri::Index(param_mesh,V1Max);
//            size_t IndexV0Min=vcg::tri::Index(param_mesh,V0Min);
//            size_t IndexV1Min=vcg::tri::Index(param_mesh,V1Min);

//            //get the direction
//            size_t Direction0=Fibres[i].FiberI[0].Direction;
//            size_t Direction1=Fibres[i].FiberI.back().Direction;
//            assert(Direction0==Direction1);

//            Point2Type DisplUp(0,0);
//            if (Direction0==0)
//                DisplUp.Y()=(DiffL/2);
//            else
//                DisplUp.X()=(DiffL/2);

//            Point2Type DisplDown=-DisplUp;

//            IndexVert.push_back(IndexV0Max);
//            TargetUV.push_back(DisplUp);
//            TargetW.push_back(W0Max);

//            IndexVert.push_back(IndexV1Max);
//            TargetUV.push_back(DisplUp);
//            TargetW.push_back(W1Max);

//            IndexVert.push_back(IndexV0Min);
//            TargetUV.push_back(DisplDown);
//            TargetW.push_back(W0Min);

//            IndexVert.push_back(IndexV1Min);
//            TargetUV.push_back(DisplDown);
//            TargetW.push_back(W1Min);
//        }
//    }

//    void SmoothDisplacements(std::vector<size_t> &IndexVert,
//                             std::vector<Point2Type> &Displ,
//                             size_t steps=1)
//    {
//        vcg::tri::UpdateSelection<MeshType>::VertexClear(param_mesh);
//        std::vector<Point2Type> CurrDispl(param_mesh.vert.size(),Point2Type(0,0));

//        //first average displacement
//        for (size_t i=0;i<IndexVert.size();i++)
//        {
//            size_t IndexV=IndexVert[i];
//            CurrDispl[IndexV]+=Displ[i];
//            param_mesh.vert[IndexV].SetS();
//        }


//        IndexVert.clear();
//        Displ.clear();
//        for (size_t i=0;i<param_mesh.vert.size();i++)
//        {
//            if (!param_mesh.vert[i].IsS())continue;
//            IndexVert.push_back(i);
//            Displ.push_back(CurrDispl[i]);
//        }


//    }



    //    vcg::tri::UpdateSelection<MeshType>::VertexClear(param_mesh);
    //    std::vector<Point2Type> CurrDispl(param_mesh.vert.size(),Point2Type(0,0));

    //    //first average displacement
    //    for (size_t i=0;i<IndexVert.size();i++)
    //    {
    //        size_t IndexV=IndexVert[i];
    //        CurrDispl[IndexV]+=Displ[i];
    //        param_mesh.vert[IndexV].SetS();
    //    }


    //    IndexVert.clear();
    //    Displ.clear();
    //    for (size_t i=0;i<param_mesh.vert.size();i++)
    //    {
    //        if (!param_mesh.vert[i].IsS())continue;
    //        IndexVert.push_back(i);
    //        Displ.push_back(CurrDispl[i]);
    //    }


//    void CumulatePerVertDispl(const std::vector<size_t> &IndexVert,
//                              const std::vector<Point2Type> &DisplV,
//                              const std::vector<ScalarType> &DisplWeight,
//                              std::vector<Point2Type> &VertDispl,
//                              std::vector<bool> &Modified)
//    {
//        VertDispl=std::vector<Point2Type>(param_mesh.vert.size(),Point2Type(0,0));
//        Modified=std::vector<bool> (param_mesh.vert.size(),false);

//        std::vector<ScalarType> VertWeight(param_mesh.vert.size(),0);
//        for (size_t i=0;i<IndexVert.size();i++)
//        {
//            size_t currV=IndexVert[i];
//            Point2Type currDispl=DisplV[i];
//            ScalarType currW=DisplWeight[i];
//            VertDispl[currV]+=(currDispl*currW);
//            VertWeight[currV]+=DisplWeight[i];
//        }

//        for (size_t i=0;i<VertDispl.size();i++)
//        {
//            if (VertWeight[i]==0)continue;
//            Modified[i]=true;
//            VertDispl[i]/=VertWeight[i];
//        }
//    }


//    void ComputeTargetDispl(std::vector<Point2Type> &PerVertTargetDispl,
//                            std::vector<bool> &Modified)
//    {
//        std::vector<size_t> IndexVert;
//        std::vector<Point2Type> TargetUV;
//        std::vector<ScalarType> TargetW;


//        ComputeTargetDispl(UFibers,IndexVert,TargetUV,TargetW);
//        ComputeTargetDispl(VFibers,IndexVert,TargetUV,TargetW);


//        CumulatePerVertDispl(IndexVert,TargetUV,TargetW,PerVertTargetDispl,Modified);

//        //SmoothDisplacements(PerVertTargetDispl,Modified);
//    }

//    void ComputeVertexDispl(const Point2Type &uvVert,
//                            Point2Type &displ)
//    {
//        ScalarType Lenght3D_U,Lenght3D_V,LenghtUV_U,LenghtUV_V;
//        InterpolateFibreLenghts(param_mesh.vert[i].T().P(),Lenght3D_U,
//                                Lenght3D_V,LenghtUV_U,LenghtUV_V);

//        ScalarType LDiffU=Lenght3D_U-LenghtUV_U;
//        ScalarType LDiffV=Lenght3D_V-LenghtUV_V;

//        CenterUV()
//    }

//    void ComputeTarget2(std::vector<Point2Type> &PerVertTargetDispl,
//                        std::vector<bool> &Modified)
//    {
//        PerVertTargetDispl=std::vector<Point2Type>(param_mesh.vert.size(),Point2Type(0,0));
//        Modified=std::vector<bool>(param_mesh.vert.size(),false);

//        for (size_t i=0;i<param_mesh.vert.size();i++)
//        {
//            if (!param_mesh.vert[i].IsB())continue;
//            Modified[i]=true;

//            ScalarType Lenght3D_U,Lenght3D_V,LenghtUV_U,LenghtUV_V;
//            InterpolateFibreLenghts(param_mesh.vert[i].T().P(),Lenght3D_U,
//                                    Lenght3D_V,LenghtUV_U,LenghtUV_V);

//            ScalarType LDiffU=Lenght3D_U-LenghtUV_U;
//            ScalarType LDiffV=Lenght3D_V-LenghtUV_V;

//            CenterUV()
//        }

//        //SmoothDisplacements(PerVertTargetDispl,Modified);
//    }


//    void ApplyTargetPos(const std::vector<Point2Type> &PerVertTargetDispl,
//                        const std::vector<bool> &Modified,
//                        ScalarType Damp)
//    {
//        assert(PerVertTargetDispl.size()==param_mesh.vert.size());
//        assert(Modified.size()==param_mesh.vert.size());


//        for (size_t i=0;i<PerVertTargetDispl.size();i++)
//        {
//            if (!Modified[i])continue;
//            param_mesh.vert[i].C()=vcg::Color4b(255,0,0,255);
//            Point2Type FinalPos=param_mesh.vert[i].T().P()+PerVertTargetDispl[i];
//            param_mesh.vert[i].T().P()=((Damp)*param_mesh.vert[i].T().P())+
//                                        ((1-Damp)*FinalPos);
//        }
//    }


//    void OptimizationStep(ScalarType Damp)
//    {
//        std::vector<bool> Modified;
//        std::vector<Point2Type> VertDispl;
//        ComputeTarget2(VertDispl,Modified);
//        //ComputeTargetDispl(VertDispl,Modified);

//        ApplyTargetPos(VertDispl,Modified,Damp);


////        vcg::tri::UpdateSelection<MeshType>::VertexClear(param_mesh);
////        for(size_t i=0;i<param_mesh.vert.size();i++)
////        {
////            if (!Modified[i])continue;
////            param_mesh.vert[i].SetS();
////        }
////        vcg::tri::OptimizeUV_LSCM(param_mesh,MeshType::VertexType::SELECTED);
//        //vcg::tri::OptimizeUV_ARAP(param_mesh,100,MeshType::VertexType::SELECTED,true);

//    }

//    bool InstanciatedFibre(const int FiberIdx,
//                           const size_t &Dir)
//    {
//        assert(Dir>=0);
//        assert(Dir<=1);
//        assert(FiberIdx>=0);

//        if (Dir==0)
//            return(FiberIdx<(int)UFibers.size());
//        else
//            return(FiberIdx<(int)VFibers.size());
//    }

//    ScalarType FiberToCoord(const int FiberIdx,
//                            const size_t &Dir)
//    {
//        assert(Dir>=0);
//        assert(Dir<=1);
//        assert(FiberIdx>=0);

//        if (Dir==0)
//        {
//            //std::cout<<"Dir "<<Dir<<std::endl;
//            //std::cout<<"FIndex "<<FiberIdx<<std::endl;
//            //std::cout<<"UVFibre "<<UFibers.size()<<std::endl;
//            assert(FiberIdx<(int)UFibers.size());
//        }
//        else
//        {
//            //std::cout<<"Dir "<<Dir<<std::endl;
//            //std::cout<<"FIndex "<<FiberIdx<<std::endl;
//            //std::cout<<"UVFibre "<<VFibers.size()<<std::endl;
//            assert(FiberIdx<VFibers.size());
//        }
//        ScalarType valMin=UVBox.min.V(Dir);
//        ScalarType valFibre=valMin+FiberIdx*fibre_space.V(Dir);
//        return valFibre;
//    }

//    void InterpolateFibreBasis(const Point2Type &sample_uv,
//                               size_t &FiberU0,size_t &FiberU1,
//                               size_t &FiberV0,size_t &FiberV1,
//                               ScalarType &alphaU,ScalarType &alphaV)
//    {
//        //get the two values
//        ScalarType MinU=UVBox.min.X();
//        ScalarType MinV=UVBox.min.Y();

//        //get the 4 surrounding fibres
//        FiberU0=floor((sample_uv.X()-MinU)/fibre_space.X());
//        FiberV0=floor((sample_uv.Y()-MinV)/fibre_space.Y());

//        FiberU1=FiberU0+1;
//        FiberV1=FiberV0+1;

//        ScalarType U0=FiberToCoord(FiberU0,0);
//        ScalarType U1=FiberToCoord(FiberU1,0);

//        ScalarType V0=FiberToCoord(FiberV0,1);
//        ScalarType V1=FiberToCoord(FiberV1,1);

//        alphaU=1-((sample_uv.X()-U0)/(U1-U0));
//        alphaV=1-((sample_uv.Y()-V0)/(V1-V0));
//    }

//    void InterpolateFibreLenghts(const Point2Type &sample_uv,
//                                 ScalarType &Lenght3D_U,
//                                 ScalarType &Lenght3D_V,
//                                 ScalarType &LenghtUV_U,
//                                 ScalarType &LenghtUV_V)
//    {

//        size_t FiberU0,FiberU1,FiberV0,FiberV1;
//        ScalarType alphaU,alphaV;

//        InterpolateFibreBasis(sample_uv,FiberU0,FiberU1,
//                              FiberV0,FiberV1,alphaU,alphaV);

//        ScalarType Lenght3D_U0=UFibers[FiberU0].Len3D;
//        ScalarType Lenght3D_U1=UFibers[FiberU1].Len3D;
//        Lenght3D_U=Lenght3D_U0*alphaU+Lenght3D_U1*(1-alphaU);

//        ScalarType LenghtUV_U0=UFibers[FiberU0].LenParam;
//        ScalarType LenghtUV_U1=UFibers[FiberU1].LenParam;
//        LenghtUV_U=LenghtUV_U0*alphaU+LenghtUV_U1*(1-alphaU);

//        ScalarType Lenght3D_V0=VFibers[FiberV0].Len3D;
//        ScalarType Lenght3D_V1=VFibers[FiberV1].Len3D;
//        Lenght3D_V=Lenght3D_V0*alphaV+Lenght3D_V1*(1-alphaV);

//        ScalarType LenghtUV_V0=VFibers[FiberV0].LenParam;
//        ScalarType LenghtUV_V1=VFibers[FiberV1].LenParam;
//        LenghtUV_V=LenghtUV_V0*alphaV+LenghtUV_V1*(1-alphaV);

//    }

//    void WhichFibers(const ScalarType &minP,
//                     const ScalarType &maxP,
//                     const size_t &Dir,
//                     int &MinIdx,int &MaxIdx)
//    {
//        assert(minP<maxP);
//        assert(Dir>=0);
//        assert(Dir<=1);
//        //get the two values
//        ScalarType MinUV=UVBox.min.V(Dir);
//        MinIdx=floor((minP-MinUV)/fibre_space.V(Dir));
//        MaxIdx=ceil((maxP-MinUV)/fibre_space.V(Dir));
//        assert(MaxIdx>MinIdx);
//    }

//    bool IntersectFibre(const Point2Type &uv0,
//                        const Point2Type &uv1,
//                        const size_t &Dir,
//                        const int &IndexFibre,
//                        ScalarType &alpha)
//    {
//        assert(Dir>=0);
//        assert(Dir<=1);

//        //get the two values
//        ScalarType val0=uv0.V(Dir);
//        ScalarType val1=uv1.V(Dir);

//        bool swapped=false;
//        if (val0>val1)
//        {
//            std::swap(val0,val1);
//            swapped=true;
//        }
//        assert(val1>=val0);

//        ScalarType FibreVal=FiberToCoord(IndexFibre,Dir);

//        if (FibreVal<val0)return false;
//        if (FibreVal>val1)return false;

//        ScalarType dist1=fabs(val1-FibreVal);
//        ScalarType lenght=val1-val0;
//        assert(lenght>=0);
//        alpha=(dist1/lenght);
//        if (swapped)alpha=1-alpha;
//        assert((alpha>0)&&(alpha<1));
//        return true;
//    }

//    bool FindFaceIntersections(const FaceType &f,
//                               const size_t &IndexFibre,
//                               const size_t &Dir,
//                               Point2Type &IntUV0,
//                               Point2Type &IntUV1,
//                               CoordType &Int3DP_0,
//                               CoordType &Int3DP_1,
//                               int &IndexE0,
//                               int &IndexE1)
//    {
//        bool Found0=false;
//        bool Found1=false;
//        for (size_t i=0;i<f.VN();i++)
//        {
//            Point2Type uv0,uv1;
//            CoordType P0,P1;
//            P0=f.cP0(i);
//            P1=f.cP1(i);
//            if (perVert)
//            {
//                uv0=f.cV(i)->cT().P();
//                uv1=f.cV((i+1)%f.VN())->cT().P();
//            }
//            else
//            {
//                uv0=f.cWT(i).P();
//                uv1=f.cWT((i+1)%f.VN()).P();
//            }
//            ScalarType alpha;
//            if (!IntersectFibre(uv0,uv1,Dir,IndexFibre,alpha))continue;
//            if (!Found0)
//            {
//                Int3DP_0=P0*alpha+P1*(1-alpha);
//                IntUV0=uv0*alpha+uv1*(1-alpha);
//                IndexE0=i;
//                Found0=true;
//            }
//            else
//            {
//                Int3DP_1=P0*alpha+P1*(1-alpha);
//                IntUV1=uv0*alpha+uv1*(1-alpha);
//                IndexE1=i;
//                Found1=true;
//            }
//        }
//        assert(Found0==Found1);
//        return (Found0 && Found1);
//    }


//    void UpdateIntersection(const size_t IndexF,
//                            const size_t &Dir)
//    {
//        FaceType *f=&param_mesh.face[IndexF];
//        Box2Type FaceUVBox;
//        if (perVert)
//        {
//            FaceUVBox.Add((*f).V(0)->T().P());
//            FaceUVBox.Add((*f).V(1)->T().P());
//            FaceUVBox.Add((*f).V(2)->T().P());
//        }
//        else
//        {
//            FaceUVBox.Add((*f).WT(0).P());
//            FaceUVBox.Add((*f).WT(1).P());
//            FaceUVBox.Add((*f).WT(2).P());
//        }

//        Point2Type uv0=FaceUVBox.min;
//        Point2Type uv1=FaceUVBox.max;

//        int MinIdx,MaxIdx;
//        WhichFibers(uv0.V(Dir),uv1.V(Dir),Dir,MinIdx,MaxIdx);
//        for (size_t j=MinIdx;j<=MaxIdx;j++)
//        {
//            if (!InstanciatedFibre(j,Dir))
//            {
//                assert(0);
//            }
//            Point2Type IntUV0,IntUV1;
//            CoordType Int3DP_0,Int3DP_1;
//            int IndexE0,IndexE1;
//            if (!FindFaceIntersections((*f),j,Dir,IntUV0,IntUV1,Int3DP_0,Int3DP_1,IndexE0,IndexE1))continue;

//            FiberInt Fint(IntUV0,IntUV1,Int3DP_0,Int3DP_1,IndexF,IndexE0,IndexE1,Dir);
//            if (Dir==0)
//            {
//                assert(j<UFibers.size());
//                UFibers[j].FiberI.push_back(Fint);
//            }
//            else
//            {
//                assert(j<VFibers.size());
//                VFibers[j].FiberI.push_back(Fint);
//            }
//        }
//    }


//    void UpdateFibreIntersection()
//    {
//        for (size_t i=0;i<param_mesh.face.size();i++)
//        {
//            UpdateIntersection(i,0);
//            UpdateIntersection(i,1);
//        }

//        for(size_t i=0;i<UFibers.size();i++)
//            std::sort(UFibers[i].FiberI.begin(),UFibers[i].FiberI.end());

//        for(size_t i=0;i<VFibers.size();i++)
//            std::sort(VFibers[i].FiberI.begin(),VFibers[i].FiberI.end());
//    }

//    void InitFibres(int resolution)
//    {
//        if (perVert)
//            UVBox=vcg::tri::UV_Utils<MeshType>::PerVertUVBox(param_mesh);
//        else
//            UVBox=vcg::tri::UV_Utils<MeshType>::PerWedgeUVBox(param_mesh);
//        //inflate a bit for numerical issues
//        UVBox.Offset(UVBox.Diag()/100);

//        //then find number of fibres
//        if (fibre_space.X()<0)
//        {
//            assert(resolution>1);
//            assert(fibre_space.Y()<0);
//            ScalarType unit_size=UVBox.Diag()/resolution;
//            fibre_space=Point2Type(unit_size,unit_size);
//        }
//        assert(fibre_space.X()>0);
//        assert(fibre_space.Y()>0);

//        size_t NumUFibersU=ceil((UVBox.max.X()-UVBox.min.X())/fibre_space.X());
//        size_t NumUFibersV=ceil((UVBox.max.Y()-UVBox.min.Y())/fibre_space.Y());

//        //then resize
//        UFibers.clear();
//        UFibers.resize(NumUFibersU+1);
//        VFibers.clear();
//        VFibers.resize(NumUFibersV+1);
//    }

//    void UpdateElongation(Fiber &f)
//    {
//        f.Len3D=0;
//        f.LenParam=0;
//        for (size_t i=0;i<f.FiberI.size();i++)
//        {
//            CoordType Pos0=f.FiberI[i].Int3DP_0;
//            CoordType Pos1=f.FiberI[i].Int3DP_1;
//            Point2Type IntUV0=f.FiberI[i].IntUV0;
//            Point2Type IntUV1=f.FiberI[i].IntUV1;
//            f.Len3D+=(Pos1-Pos0).Norm();
//            f.LenParam+=(IntUV0-IntUV1).Norm();
//        }
//        f.RatioL=f.Len3D/f.LenParam;
//        f.DiffL=(f.Len3D-f.LenParam);
//    }

//    void UpdateElongation()
//    {
//        MaxElongR=0;
//        MaxCompressR=std::numeric_limits<ScalarType>::max();
//        MaxElongL=0;
//        MaxCompressL=std::numeric_limits<ScalarType>::max();

//        TotDiffL=0;

//        for (size_t i=0;i<UFibers.size();i++)
//        {
//            UpdateElongation(UFibers[i]);
//            MaxElongR=std::max(MaxElongR,UFibers[i].RatioL);
//            MaxCompressR=std::min(MaxCompressR,UFibers[i].RatioL);

//            MaxElongL=std::max(MaxElongL,UFibers[i].DiffL);
//            MaxCompressL=std::min(MaxCompressL,UFibers[i].DiffL);

//            TotDiffL+=fabs(VFibers[i].DiffL);
//        }
//        for (size_t i=0;i<VFibers.size();i++)
//        {
//            UpdateElongation(VFibers[i]);
//            MaxElongR=std::max(MaxElongR,VFibers[i].RatioL);
//            MaxCompressR=std::min(MaxCompressR,VFibers[i].RatioL);

//            MaxElongL=std::max(MaxElongL,VFibers[i].DiffL);
//            MaxCompressL=std::min(MaxCompressL,VFibers[i].DiffL);

//            TotDiffL+=fabs(VFibers[i].DiffL);
//        }
//        if (WriteDebug)
//        {
//            std::cout<<"Max Elongation Ratio "<<MaxElongR<<std::endl;
//            std::cout<<"Max Compression Ratio "<<MaxCompressR<<std::endl;
//            std::cout<<"Max Elongation Lenght "<<MaxElongL<<std::endl;
//            std::cout<<"Max Compression Lenght "<<MaxCompressL<<std::endl;
//        }
//    }

//    void GetFiberMesh(Fiber &f,
//                      MeshType &FibreMesh3D,
//                      MeshType &FibreMeshUV,
//                      ScalarType radius)
//    {
//        FibreMesh3D.Clear();
//        FibreMeshUV.Clear();
//        for (size_t i=0;i<f.FiberI.size();i++)
//        {
//            CoordType Pos0=f.FiberI[i].Int3DP_0;
//            CoordType Pos1=f.FiberI[i].Int3DP_1;
//            CoordType PosUV0=CoordType(f.FiberI[i].IntUV0.X(),
//                                       f.FiberI[i].IntUV0.Y(),0);
//            CoordType PosUV1=CoordType(f.FiberI[i].IntUV1.X(),
//                                       f.FiberI[i].IntUV1.Y(),0);
//            MeshType Seg3D,SegUV;
//            vcg::tri::OrientedCylinder(Seg3D,Pos0,Pos1,radius,true,8,2);
//            vcg::tri::OrientedCylinder(SegUV,PosUV0,PosUV1,radius,true,8,2);
//            vcg::tri::Append<MeshType,MeshType>::Mesh(FibreMesh3D,Seg3D);
//            vcg::tri::Append<MeshType,MeshType>::Mesh(FibreMeshUV,SegUV);
//        }
//    }

//    void ExtractVectFiberMeshes(MeshType &Mesh3D,MeshType &MeshUV,
//                                std::vector<Fiber> Fibers,
//                                ScalarType MaxElongRatio=-1,
//                                ScalarType MaxComprRatio=-1,
//                                ScalarType radius=0.5)
//    {
//        assert(MaxCompressR<=1);
//        assert(MaxElongR>=1);
//        ScalarType diff_to_1_Compr=fabs(1-MaxCompressR);
//        ScalarType diff_to_1_Elong=fabs(1-MaxElongR);
//        ScalarType diff_to_1=std::max(diff_to_1_Compr,diff_to_1_Elong);

//        if (MaxComprRatio==-1)
//            MaxComprRatio=1-diff_to_1;
//        if (MaxElongRatio==-1)
//            MaxElongRatio=1+diff_to_1;

//        assert(MaxElongRatio>=1);
//        assert(MaxComprRatio<=1);

//        //        std::cout<<"Min Interval Color "<<MaxComprRatio<<std::endl;
//        //        std::cout<<"Max Interval Color "<<MaxElongRatio<<std::endl;

//        Mesh3D.Clear();
//        MeshUV.Clear();
//        ScalarType abs_radius=fibre_space.Norm()*radius;
//        for (size_t i=0;i<Fibers.size();i++)
//        {
//            MeshType FibreMesh3D,FibreMeshUV;
//            GetFiberMesh(Fibers[i],FibreMesh3D,FibreMeshUV,abs_radius);
//            ScalarType Fibre_RatioL=Fibers[i].RatioL;

//            vcg::Color4b col;
//            col=vcg::Color4b::ColorRamp(MaxComprRatio,MaxElongRatio,Fibre_RatioL);

//            vcg::tri::UpdateColor<MeshType>::PerFaceConstant(FibreMesh3D,col);
//            vcg::tri::UpdateColor<MeshType>::PerFaceConstant(FibreMeshUV,col);

//            vcg::tri::Append<MeshType,MeshType>::Mesh(Mesh3D,FibreMesh3D);
//            vcg::tri::Append<MeshType,MeshType>::Mesh(MeshUV,FibreMeshUV);
//        }
//    }

//    ScalarType Area()
//    {
//        ScalarType A=0;
//        for (size_t i=0;i<param_mesh.face.size();i++)
//            A+=vcg::DoubleArea(param_mesh.face[i]);
//        return A/2;
//    }

//public:

//    LenghMode LMode;

//    void ExtractFiberMeshes(MeshType &Mesh3D,MeshType &MeshUV,
//                            ScalarType MaxElongRatio=-1,
//                            ScalarType MaxComprRatio=-1,
//                            ScalarType radius=0.01)
//    {
//        Mesh3D.Clear();
//        MeshUV.Clear();
//        vcg::tri::UpdateBounding<MeshType>::Box(param_mesh);
//        ScalarType abs_radius=sqrt(Area())*radius;

//        MeshType FibreMesh3D_U,FibreMesh3D_V;
//        MeshType FibreMeshUV_U,FibreMeshUV_V;

//        ExtractVectFiberMeshes(FibreMesh3D_U,FibreMeshUV_U,UFibers,
//                               MaxElongRatio,MaxComprRatio,abs_radius);

//        ExtractVectFiberMeshes(FibreMesh3D_V,FibreMeshUV_V,VFibers,
//                               MaxElongRatio,MaxComprRatio,abs_radius);

//        vcg::tri::Append<MeshType,MeshType>::Mesh(Mesh3D,FibreMesh3D_U);
//        vcg::tri::Append<MeshType,MeshType>::Mesh(Mesh3D,FibreMesh3D_V);

//        vcg::tri::Append<MeshType,MeshType>::Mesh(MeshUV,FibreMeshUV_U);
//        vcg::tri::Append<MeshType,MeshType>::Mesh(MeshUV,FibreMeshUV_V);
//    }

//    void ComputeFibreElongation(int resolution=100,
//                                Point2Type _fibre_space=Point2Type(-1,-1),
//                                bool _perVert=true)
//    {
//        perVert=_perVert;
//        fibre_space=_fibre_space;
//        InitFibres(resolution);
//        UpdateFibreIntersection();
//        UpdateElongation();
//    }


//    void FibreOptimization(int resolution=100,
//                           size_t steps=10,
//                           ScalarType Damp=0.5)
//    {

//        //        std::vector<Point2Type> TargetPos;
//        //        std::vector<bool> Optimized;
//        ComputeFibreElongation(resolution);
//        ScalarType OldEnergy=TotDiffL;//std::max(MaxElongR,((ScalarType)1)/MaxCompressR);
//        std::cout<<"Energy: "<<OldEnergy<<std::endl;
//        for (size_t s=0;s<steps;s++)
//        {
//            OptimizationStep(Damp);

//            ComputeFibreElongation(resolution);
//            ScalarType NewEnergy=TotDiffL;//std::max(MaxElongR,((ScalarType)1)/MaxCompressR);
//            std::cout<<"Step: "<<s<<std::endl;
//            std::cout<<"Energy: "<<NewEnergy<<std::endl;
//            if (NewEnergy>OldEnergy)break;
//            OldEnergy=NewEnergy;
//        }
//    }

//    FibreElongation(MeshType &_param_mesh):param_mesh(_param_mesh)
//    {
//        LMode=LModeGLobal;
//        WriteDebug=false;
//        MaxScale=2;
//    }
//};



enum LenghMode{LModeLocal,LModeGLobal};

template <class MeshType>
class FibreElongation
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;

    typedef typename vcg::Point2<ScalarType> Point2Type;
    typedef typename vcg::Box2<ScalarType> Box2Type;

    MeshType &param_mesh;
    bool WriteDebug;
    ScalarType MaxScale;

    struct FiberInt
    {
        Point2Type IntUV0,IntUV1;
        CoordType Int3DP_0,Int3DP_1;
        size_t FaceIdx;
        size_t IndexE0,IndexE1;
        size_t Direction;


        ScalarType Lenght3D()const
        {
            return (Int3DP_0 - Int3DP_1).Norm();
        }

        ScalarType LenghtUV()const
        {
            return (IntUV0 - IntUV1).Norm();
        }

        inline bool operator <(const FiberInt &FInt)const
        {
            Point2Type AvgUV0=(IntUV0+IntUV1)/2;
            Point2Type AvgUV1=(FInt.IntUV0+FInt.IntUV1)/2;
            if (Direction==0)return (AvgUV0.Y()<AvgUV1.Y());
            assert(Direction==1);
            return (AvgUV0.X()<AvgUV1.X());
        }

        FiberInt(Point2Type &_IntUV0,
                 Point2Type &_IntUV1,
                 CoordType &_Int3DP_0,
                 CoordType &_Int3DP_1,
                 size_t _FaceIdx,
                 size_t _IndexE0,
                 size_t _IndexE1,
                 size_t _Direction)
        {
            IntUV0=_IntUV0;
            IntUV1=_IntUV1;
            Int3DP_0=_Int3DP_0;
            Int3DP_1=_Int3DP_1;
            FaceIdx=_FaceIdx;
            Direction=_Direction;
            IndexE0=_IndexE0;
            IndexE1=_IndexE1;
            if (((Direction==0)&&(IntUV0.Y()>IntUV1.Y()))||
                    ((Direction==1)&&(IntUV0.X()>IntUV1.X())))
            {
                std::swap(IntUV0,IntUV1);
                std::swap(Int3DP_0,Int3DP_1);
                std::swap(IndexE0,IndexE1);
            }
        }

    };

    struct Fiber
    {
        ScalarType LenParam;
        ScalarType Len3D;
        ScalarType RatioL;
        ScalarType DiffL;
        std::vector<FiberInt> FiberI;


        bool UVExtremes(Point2Type &P0,Point2Type &P1)
        {
            if (FiberI.size()==0)return false;
            P0=FiberI[0].IntUV0;
            P1=FiberI.back().IntUV1;
        }

        Point2Type CenterFibre()
        {
            Point2Type P0,P1;
            UVExtremes(P0,P1);
            return ((P0+P1)/2);
        }

        void PrintDebug()
        {
            Point2Type P0,P1;
            bool cohincident=UVExtremes(P0,P1);
            std::cout<<"Lenght Fibre 3D:"<<Len3D<<std::endl;
            std::cout<<"Lenght Fibre UV:"<<LenParam<<std::endl;
            std::cout<<"Extremes:"<<"("<<P0.X()<<","<<P0.Y()<<")"
                    <<"("<<P1.X()<<","<<P1.Y()<<")"<<std::endl;
            std::cout<<"Lenght Computed:"<<(P0-P1).Norm()<<std::endl;
        }

        Fiber()
        {LenParam=0;Len3D=0;RatioL=0;DiffL=0;}
    };

    std::vector<Fiber> UFibers,VFibers;
    Box2Type UVBox;
    bool perVert;
    Point2Type fibre_space;
public:
    ScalarType MaxElongR;
    ScalarType MaxCompressR;
    ScalarType MaxElongL;
    ScalarType MaxCompressL;
private:
    ScalarType TotDiffL;

    ScalarType MaxDispl;

    //    std::vector<Point2Type> VertDirectionU;
    //    std::vector<Point2Type> VertDirectionV;

    void InitMaxDispl()
    {
        MaxDispl=0;
        size_t Num=0;
        for (size_t i=0;i<param_mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                Point2Type P0=param_mesh.face[i].V0(j)->T().P();
                Point2Type P1=param_mesh.face[i].V1(j)->T().P();
                MaxDispl+=(P0-P1).Norm();
                Num++;
            }
        MaxDispl/=Num;
        MaxDispl/=2;
    }

    void SmoothDisplacements(std::vector<Point2Type> &Displ,
                             size_t steps=1)
    {
        for (size_t s=0;s<steps;s++)
        {
            std::vector<Point2Type> DisplSum=Displ;
            std::vector<size_t> DisplNum(param_mesh.vert.size(),1);

            for (size_t i=0;i<param_mesh.face.size();i++)
                for (size_t j=0;j<3;j++)
                {
                    size_t IndexV0=vcg::tri::Index(param_mesh,param_mesh.face[i].V0(j));
                    size_t IndexV1=vcg::tri::Index(param_mesh,param_mesh.face[i].V1(j));
                    DisplNum[IndexV0]++;
                    DisplNum[IndexV1]++;
                    DisplSum[IndexV0]+=Displ[IndexV1];
                    DisplSum[IndexV1]+=Displ[IndexV0];
                }

            for (size_t i=0;i<DisplSum.size();i++)
            {
                DisplSum[i]/=DisplNum[i];
                Displ[i]=DisplSum[i];
            }

        }
    }
    //    void InitVectDirection()
    //    {
    //        VertDirectionV.resize(param_mesh.vert.size(),Point2Type(0,0));
    //        VertDirectionU.resize(param_mesh.vert.size(),Point2Type(0,0));
    //        for (size_t i=0;i<param_mesh.face.size();i++)
    //            for (size_t j=0;j<3;j++)
    //            {
    //                if (!vcg::face::IsBorder(param_mesh.face[i],j))continue;
    //                Point2Type P0=param_mesh.face[i].V0(j)->T().P();
    //                Point2Type P1=param_mesh.face[i].V1(j)->T().P();
    //                Point2Type Dir=(P1-P0).Normalize();
    //                Point2Type OrthoDir(-Dir.Y(),Dir.X());
    //                size_t IndexV0=vcg::tri::Index(param_mesh,param_mesh.face[i].V0(j));
    //                size_t IndexV1=vcg::tri::Index(param_mesh,param_mesh.face[i].V1(j));
    //                VertDirectionU[IndexV0]+=OrthoDir;
    //                VertDirectionV[IndexV1]+=OrthoDir;
    //            }

    //        for (size_t i=0;i<VertDirectionV.size();i++)
    //        {
    //            if (VertDirectionU[i].X()>0)
    //                VertDirectionU[i]=Point2Type(1,0);
    //            else
    //                VertDirectionU[i]=Point2Type(-1,0);

    //            if (VertDirectionV[i].Y()>0)
    //                VertDirectionV[i]=Point2Type(0,1);
    //            else
    //                VertDirectionV[i]=Point2Type(0,-1);
    //        }
    //    }
    //    std::vector<std::vector<Point2Type> > StretchDir;

    //    void CumputeStretchDir()
    //    {
    //       for (size_t i=0;i<UFibers.size();i++)
    //           for (size_t i=0;i<VFibers.size();i++)
    //           {

    //           }
    //           ScalarType FiberToCoord(const int FiberIdx,
    //                                   const size_t &Dir)
    //           {
    //       }
    //    }

    //    Point2Type GetDirection(Point2Type &pos,
    //                            )
    //    {

    //    }

    //    void ComputeTargetDispl(std::vector<Point2Type> &PerVertTargetDispl,
    //                       std::vector<bool> &Modified)
    //    {
    //        PerVertTargetDispl=std::vector<Point2Type>(param_mesh.vert.size(),Point2Type(0,0));
    //        Modified=std::vector<bool>(param_mesh.vert.size(),false);

    //        //ScalarType val=vcg::tri::UV_Utils<MeshType>::PerVertUVBox(param_mesh).Diag()*0.01;
    //        for (size_t i=0;i<param_mesh.vert.size();i++)
    //        {
    //            //std::cout<<"Test"<<std::endl;
    //            if (!param_mesh.vert[i].IsB())continue;

    //            Modified[i]=true;

    //            int FiberU0,FiberU1,FiberV0,FiberV1;
    //            Point2Type UVvert=param_mesh.vert[i].T().P();
    //            WhichFibers(UVvert,FiberU0,FiberU1,FiberV0,FiberV1);

    //            //            std::cout<<"Val U:"<<UVvert.X()<<std::endl;
    //            //            std::cout<<"Val V:"<<UVvert.Y()<<std::endl;

    //            //show the uv top e bottom
    //            Point2Type ExtremeU00,ExtremeU01;
    //            Point2Type ExtremeU10,ExtremeU11;

    //            UFibers[FiberU0].UVExtremes(ExtremeU00,ExtremeU01);
    //            UFibers[FiberU1].UVExtremes(ExtremeU10,ExtremeU11);

    //            //            std::cout<<"Fibre U0"<<std::endl;
    //            //            std::cout<<"P0:"<<ExtremeU00.X()<<","<<ExtremeU00.Y()<<std::endl;
    //            //            std::cout<<"P1:"<<ExtremeU01.X()<<","<<ExtremeU01.Y()<<std::endl;

    //            //            std::cout<<"Fibre U1"<<std::endl;
    //            //            std::cout<<"P0:"<<ExtremeU10.X()<<","<<ExtremeU10.Y()<<std::endl;
    //            //            std::cout<<"P1:"<<ExtremeU11.X()<<","<<ExtremeU11.Y()<<std::endl;


    //            //then get the distance wrt fibres
    //            ScalarType UVal0=FiberToCoord(FiberU0,0);
    //            ScalarType UVal1=FiberToCoord(FiberU1,0);

    //            ScalarType VVal0=FiberToCoord(FiberV0,1);
    //            ScalarType VVal1=FiberToCoord(FiberV1,1);

    //            //            std::cout<<"Val U0:"<<UVal0<<std::endl;
    //            //            std::cout<<"Val U1:"<<UVal1<<std::endl;

    //            ScalarType alphaU=1-(UVvert.X()-UVal0)/(UVal1-UVal0);
    //            ScalarType alphaV=1-(UVvert.Y()-VVal0)/(VVal1-VVal0);

    //            //diff must be inverted because should be lenght in UV - lenght in 3D
    //            ScalarType DiffLenU0=-UFibers[FiberU0].DiffL/2;
    //            ScalarType DiffLenU1=-UFibers[FiberU1].DiffL/2;
    //            ScalarType DiffLenInterpU=DiffLenU0*alphaU+DiffLenU1*(1-alphaU);

    //            //diff must be inverted because should be lenght in UV - lenght in 3D
    //            ScalarType DiffLenV0=-VFibers[FiberV0].DiffL/2;
    //            ScalarType DiffLenV1=-VFibers[FiberV1].DiffL/2;
    //            ScalarType DiffLenInterpV=DiffLenV0*alphaV+DiffLenV1*(1-alphaV);

    //            //notice that U and V are flipped on purposes here
    //            Point2Type DirDisplU=VertDirectionU[i];
    //            DirDisplU*=DiffLenInterpV;

    //            Point2Type DirDisplV=VertDirectionV[i];
    //            DirDisplV*=DiffLenInterpU;

    //            PerVertTargetDispl[i]=(DirDisplU+DirDisplV);
    //        }

    //        //SmoothDisplacements(PerVertTargetDispl,Modified);
    //    }

    void ComputeTargetDispl(std::vector<Point2Type> &PerVertTargetDispl,
                            std::vector<bool> &Modified)
    {
        PerVertTargetDispl=std::vector<Point2Type>(param_mesh.vert.size(),Point2Type(0,0));
        Modified=std::vector<bool>(param_mesh.vert.size(),false);

        //ScalarType val=vcg::tri::UV_Utils<MeshType>::PerVertUVBox(param_mesh).Diag()*0.01;
        for (size_t i=0;i<param_mesh.vert.size();i++)
        {
            //std::cout<<"Test"<<std::endl;
            //if (!param_mesh.vert[i].IsB())continue;

            Modified[i]=true;

            int FiberU0,FiberU1,FiberV0,FiberV1;
            Point2Type UVvert=param_mesh.vert[i].T().P();
            WhichFibers(UVvert,FiberU0,FiberU1,FiberV0,FiberV1);

            //                        std::cout<<"Val U:"<<UVvert.X()<<std::endl;
            //                        std::cout<<"Val V:"<<UVvert.Y()<<std::endl;

            //show the uv top e bottom
            Point2Type ExtremeU00,ExtremeU01;
            Point2Type ExtremeU10,ExtremeU11;

            UFibers[FiberU0].UVExtremes(ExtremeU00,ExtremeU01);
            UFibers[FiberU1].UVExtremes(ExtremeU10,ExtremeU11);

            //            std::cout<<"Fibre U0"<<std::endl;
            //            std::cout<<"P0:"<<ExtremeU00.X()<<","<<ExtremeU00.Y()<<std::endl;
            //            std::cout<<"P1:"<<ExtremeU01.X()<<","<<ExtremeU01.Y()<<std::endl;

            //            std::cout<<"Fibre U1"<<std::endl;
            //            std::cout<<"P0:"<<ExtremeU10.X()<<","<<ExtremeU10.Y()<<std::endl;
            //            std::cout<<"P1:"<<ExtremeU11.X()<<","<<ExtremeU11.Y()<<std::endl;


            //then get the distance wrt fibres
            ScalarType UVal0=FiberToCoord(FiberU0,0);
            ScalarType UVal1=FiberToCoord(FiberU1,0);

            ScalarType VVal0=FiberToCoord(FiberV0,1);
            ScalarType VVal1=FiberToCoord(FiberV1,1);

            //            std::cout<<"Val U0:"<<UVal0<<std::endl;
            //            std::cout<<"Val U1:"<<UVal1<<std::endl;

            ScalarType alphaU=1-((UVvert.X()-UVal0)/(UVal1-UVal0));
            ScalarType alphaV=1-((UVvert.Y()-VVal0)/(VVal1-VVal0));

            ScalarType DiffLenU0=UFibers[FiberU0].DiffL/2;
            ScalarType DiffLenU1=UFibers[FiberU1].DiffL/2;
            ScalarType DiffLenInterpU=DiffLenU0*alphaU+DiffLenU1*(1-alphaU);

            //then get the lenght intepolated
            ScalarType LenU0=UFibers[FiberU0].LenParam;
            ScalarType LenU1=UFibers[FiberU1].LenParam;
            ScalarType LenInterpU=LenU0*alphaU+LenU1*(1-alphaU);

            //then get the interpolated bary
            ScalarType CentreU0=UFibers[FiberU0].CenterFibre()[1];
            ScalarType CentreU1=UFibers[FiberU1].CenterFibre()[1];
            ScalarType CentreInterpU=CentreU0*alphaU+CentreU1*(1-alphaU);
            ScalarType ScaleValU=(UVvert.V(1)-CentreInterpU)/(LenInterpU/2);
            ScaleValU=std::min(ScaleValU,(ScalarType)1);
            ScaleValU=std::max(ScaleValU,(ScalarType)-1);
            DiffLenInterpU*=ScaleValU;

            ScalarType DiffLenV0=VFibers[FiberV0].DiffL/2;
            ScalarType DiffLenV1=VFibers[FiberV1].DiffL/2;
            ScalarType DiffLenInterpV=DiffLenV0*alphaV+DiffLenV1*(1-alphaV);

            //then get the lenght intepolated
            ScalarType LenV0=VFibers[FiberV0].LenParam;
            ScalarType LenV1=VFibers[FiberV1].LenParam;
            ScalarType LenInterpV=LenV0*alphaV+LenV1*(1-alphaV);

            //then get the interpolated bary
            ScalarType CentreV0=VFibers[FiberV0].CenterFibre()[0];
            ScalarType CentreV1=VFibers[FiberV1].CenterFibre()[0];
            ScalarType CentreInterpV=CentreV0*alphaV+CentreV1*(1-alphaV);
            ScalarType ScaleValV=(UVvert.V(0)-CentreInterpV)/(LenInterpV/2);
            ScaleValV=std::min(ScaleValV,(ScalarType)1);
            ScaleValV=std::max(ScaleValV,(ScalarType)-1);
            DiffLenInterpV*=ScaleValV;


            Point2Type DirDisplU=Point2Type(0,1);//VertDirectionU[i];
            DirDisplU*=DiffLenInterpU;

            Point2Type DirDisplV=Point2Type(1,0);//VertDirectionU[i];
            DirDisplV*=DiffLenInterpV;

            PerVertTargetDispl[i]=(DirDisplU+DirDisplV);
        }

        //SmoothDisplacements(PerVertTargetDispl,Modified);
    }

    void ScaleTargetDispl(std::vector<Point2Type> &PerVertTargetDispl)
    {
        //normalize the diplacement by max
        ScalarType MaxNorm=0;
        for (size_t i=0;i<PerVertTargetDispl.size();i++)
            MaxNorm=std::max(MaxNorm,PerVertTargetDispl[i].Norm());

        if (MaxNorm<MaxDispl)return;

        for (size_t i=0;i<PerVertTargetDispl.size();i++)
            PerVertTargetDispl[i]/=MaxNorm;

        for (size_t i=0;i<PerVertTargetDispl.size();i++)
            PerVertTargetDispl[i]*=MaxDispl;
    }

    void ApplyTargetDispl(const std::vector<Point2Type> &PerVertTargetDispl,
                          const std::vector<bool> &Modified,
                          ScalarType Damp)
    {
        assert(PerVertTargetDispl.size()==param_mesh.vert.size());
        assert(Modified.size()==param_mesh.vert.size());

        for (size_t i=0;i<PerVertTargetDispl.size();i++)
        {
            if (!Modified[i])continue;
            //            param_mesh.vert[i].C()=vcg::Color4b(255,0,0,255);
            Point2Type FinalPos=param_mesh.vert[i].T().P()+PerVertTargetDispl[i];
            param_mesh.vert[i].T().P()=((Damp)*param_mesh.vert[i].T().P())+
                    ((1-Damp)*FinalPos);
        }
    }


    void OptimizationStep(ScalarType Damp)
    {
        std::vector<bool> Modified;
        std::vector<Point2Type> VertDispl;

        ComputeTargetDispl(VertDispl,Modified);

        ScaleTargetDispl(VertDispl);

        SmoothDisplacements(VertDispl);

        ApplyTargetDispl(VertDispl,Modified,Damp);


        //vcg::tri::OptimizeUV_LSCM(param_mesh,100,MeshType::VertexType::SELECTED,false);
        //        vcg::tri::UpdateSelection<MeshType>::VertexClear(param_mesh);
        //        for(size_t i=0;i<param_mesh.vert.size();i++)
        //        {
        //            if (!Modified[i])continue;
        //            param_mesh.vert[i].SetS();
        //        }
        //        //vcg::tri::InitializeArapWithLSCM<MyTriMesh>(param_mesh);
        //        //vcg::tri::OptimizeUV_ARAP(tri_mesh);
        //        //vcg::tri::OptimizeUV_LSCM(param_mesh,MeshType::VertexType::SELECTED);

        //        vcg::tri::OptimizeUV_ARAP(param_mesh,100,MeshType::VertexType::SELECTED,false);

        ////        vcg::tri::UpdateSelection<MeshType>::VertexInvert(param_mesh);

        ////        vcg::tri::OptimizeUV_ARAP(param_mesh,100,MeshType::VertexType::SELECTED,false);

    }

    bool InstanciatedFibre(const int FiberIdx,
                           const size_t &Dir)
    {
        assert(Dir>=0);
        assert(Dir<=1);
        assert(FiberIdx>=0);

        if (Dir==0)
            return(FiberIdx<(int)UFibers.size());
        else
            return(FiberIdx<(int)VFibers.size());
    }

    ScalarType FiberToCoord(const int FiberIdx,
                            const size_t &Dir)
    {
        assert(Dir>=0);
        assert(Dir<=1);
        assert(FiberIdx>=0);

        if (Dir==0)
        {
            //std::cout<<"Dir "<<Dir<<std::endl;
            //std::cout<<"FIndex "<<FiberIdx<<std::endl;
            //std::cout<<"UVFibre "<<UFibers.size()<<std::endl;
            assert(FiberIdx<(int)UFibers.size());
        }
        else
        {
            //std::cout<<"Dir "<<Dir<<std::endl;
            //std::cout<<"FIndex "<<FiberIdx<<std::endl;
            //std::cout<<"UVFibre "<<VFibers.size()<<std::endl;
            assert(FiberIdx<VFibers.size());
        }
        ScalarType valMin=UVBox.min.V(Dir);
        ScalarType valFibre=valMin+FiberIdx*fibre_space.V(Dir);
        return valFibre;
    }

    //    void InterpolateFibreBasis(const Point2Type &sample_uv,
    //                               size_t &FiberU0,size_t &FiberU1,
    //                               size_t &FiberV0,size_t &FiberV1,
    //                               ScalarType &alphaU,ScalarType &alphaV)
    //    {
    //        //get the two values
    //        ScalarType MinU=UVBox.min.X();
    //        ScalarType MinV=UVBox.min.Y();

    //        //get the 4 surrounding fibres
    //        FiberU0=floor((sample_uv.X()-MinU)/fibre_space.X());
    //        FiberV0=floor((sample_uv.Y()-MinV)/fibre_space.Y());

    //        FiberU1=FiberU0+1;
    //        FiberV1=FiberV0+1;

    //        ScalarType U0=FiberToCoord(FiberU0,0);
    //        ScalarType U1=FiberToCoord(FiberU1,0);

    //        ScalarType V0=FiberToCoord(FiberV0,1);
    //        ScalarType V1=FiberToCoord(FiberV1,1);

    //        alphaU=1-((sample_uv.X()-U0)/(U1-U0));
    //        alphaV=1-((sample_uv.Y()-V0)/(V1-V0));
    //    }

    //    void InterpolateFibreLenghts(const Point2Type &sample_uv,
    //                                 ScalarType &Lenght3D_U,
    //                                 ScalarType &Lenght3D_V,
    //                                 ScalarType &LenghtUV_U,
    //                                 ScalarType &LenghtUV_V)
    //    {

    //        size_t FiberU0,FiberU1,FiberV0,FiberV1;
    //        ScalarType alphaU,alphaV;

    //        InterpolateFibreBasis(sample_uv,FiberU0,FiberU1,
    //                              FiberV0,FiberV1,alphaU,alphaV);

    //        ScalarType Lenght3D_U0=UFibers[FiberU0].Len3D;
    //        ScalarType Lenght3D_U1=UFibers[FiberU1].Len3D;
    //        Lenght3D_U=Lenght3D_U0*alphaU+Lenght3D_U1*(1-alphaU);

    //        ScalarType LenghtUV_U0=UFibers[FiberU0].LenParam;
    //        ScalarType LenghtUV_U1=UFibers[FiberU1].LenParam;
    //        LenghtUV_U=LenghtUV_U0*alphaU+LenghtUV_U1*(1-alphaU);

    //        ScalarType Lenght3D_V0=VFibers[FiberV0].Len3D;
    //        ScalarType Lenght3D_V1=VFibers[FiberV1].Len3D;
    //        Lenght3D_V=Lenght3D_V0*alphaV+Lenght3D_V1*(1-alphaV);

    //        ScalarType LenghtUV_V0=VFibers[FiberV0].LenParam;
    //        ScalarType LenghtUV_V1=VFibers[FiberV1].LenParam;
    //        LenghtUV_V=LenghtUV_V0*alphaV+LenghtUV_V1*(1-alphaV);

    //    }



    void WhichFibers(const ScalarType &minP,const ScalarType &maxP,
                     const size_t &Dir,int &MinIdx,int &MaxIdx)
    {
        assert(minP<=maxP);
        assert(Dir>=0);
        assert(Dir<=1);
        //get the two values
        ScalarType MinUV=UVBox.min.V(Dir);
        MinIdx=floor((minP-MinUV)/fibre_space.V(Dir));
        MaxIdx=ceil((maxP-MinUV)/fibre_space.V(Dir));
        assert(MaxIdx>MinIdx);
    }

    void WhichFibers(const Point2Type &sample_uv,
                     int &FiberU0,int &FiberU1,
                     int &FiberV0,int &FiberV1)
    {
        WhichFibers(sample_uv.X(),sample_uv.X(),0,FiberU0,FiberU1);
        WhichFibers(sample_uv.Y(),sample_uv.Y(),1,FiberV0,FiberV1);
    }

    bool IntersectFibre(const Point2Type &uv0,
                        const Point2Type &uv1,
                        const size_t &Dir,
                        const int &IndexFibre,
                        ScalarType &alpha)
    {
        assert(Dir>=0);
        assert(Dir<=1);

        //get the two values
        ScalarType val0=uv0.V(Dir);
        ScalarType val1=uv1.V(Dir);

        bool swapped=false;
        if (val0>val1)
        {
            std::swap(val0,val1);
            swapped=true;
        }
        assert(val1>=val0);

        ScalarType FibreVal=FiberToCoord(IndexFibre,Dir);

        if (FibreVal<val0)return false;
        if (FibreVal>val1)return false;

        ScalarType dist1=fabs(val1-FibreVal);
        ScalarType lenght=val1-val0;
        assert(lenght>=0);
        alpha=(dist1/lenght);
        if (swapped)alpha=1-alpha;
        assert((alpha>0)&&(alpha<1));
        return true;
    }

    bool FindFaceIntersections(const FaceType &f,
                               const size_t &IndexFibre,
                               const size_t &Dir,
                               Point2Type &IntUV0,
                               Point2Type &IntUV1,
                               CoordType &Int3DP_0,
                               CoordType &Int3DP_1,
                               int &IndexE0,
                               int &IndexE1)
    {
        bool Found0=false;
        bool Found1=false;
        for (size_t i=0;i<f.VN();i++)
        {
            Point2Type uv0,uv1;
            CoordType P0,P1;
            P0=f.cP0(i);
            P1=f.cP1(i);
            if (perVert)
            {
                uv0=f.cV(i)->cT().P();
                uv1=f.cV((i+1)%f.VN())->cT().P();
            }
            else
            {
                uv0=f.cWT(i).P();
                uv1=f.cWT((i+1)%f.VN()).P();
            }
            ScalarType alpha;
            if (!IntersectFibre(uv0,uv1,Dir,IndexFibre,alpha))continue;
            if (!Found0)
            {
                Int3DP_0=P0*alpha+P1*(1-alpha);
                IntUV0=uv0*alpha+uv1*(1-alpha);
                IndexE0=i;
                Found0=true;
            }
            else
            {
                Int3DP_1=P0*alpha+P1*(1-alpha);
                IntUV1=uv0*alpha+uv1*(1-alpha);
                IndexE1=i;
                Found1=true;
            }
        }
        assert(Found0==Found1);
        return (Found0 && Found1);
    }


    void UpdateIntersection(const size_t IndexF,
                            const size_t &Dir)
    {
        FaceType *f=&param_mesh.face[IndexF];
        Box2Type FaceUVBox;
        if (perVert)
        {
            FaceUVBox.Add((*f).V(0)->T().P());
            FaceUVBox.Add((*f).V(1)->T().P());
            FaceUVBox.Add((*f).V(2)->T().P());
        }
        else
        {
            FaceUVBox.Add((*f).WT(0).P());
            FaceUVBox.Add((*f).WT(1).P());
            FaceUVBox.Add((*f).WT(2).P());
        }

        Point2Type uv0=FaceUVBox.min;
        Point2Type uv1=FaceUVBox.max;

        int MinIdx,MaxIdx;
        WhichFibers(uv0.V(Dir),uv1.V(Dir),Dir,MinIdx,MaxIdx);
        for (size_t j=MinIdx;j<=MaxIdx;j++)
        {
            if (!InstanciatedFibre(j,Dir))
            {
                assert(0);
            }
            Point2Type IntUV0,IntUV1;
            CoordType Int3DP_0,Int3DP_1;
            int IndexE0,IndexE1;
            if (!FindFaceIntersections((*f),j,Dir,IntUV0,IntUV1,Int3DP_0,Int3DP_1,IndexE0,IndexE1))continue;

            FiberInt Fint(IntUV0,IntUV1,Int3DP_0,Int3DP_1,IndexF,IndexE0,IndexE1,Dir);
            if (Dir==0)
            {
                assert(j<UFibers.size());
                UFibers[j].FiberI.push_back(Fint);
            }
            else
            {
                assert(j<VFibers.size());
                VFibers[j].FiberI.push_back(Fint);
            }
        }
    }


    void UpdateFibreIntersection()
    {
        for (size_t i=0;i<param_mesh.face.size();i++)
        {
            UpdateIntersection(i,0);
            UpdateIntersection(i,1);
        }

        for(size_t i=0;i<UFibers.size();i++)
            std::sort(UFibers[i].FiberI.begin(),UFibers[i].FiberI.end());

        for(size_t i=0;i<VFibers.size();i++)
            std::sort(VFibers[i].FiberI.begin(),VFibers[i].FiberI.end());
    }

    void InitFibres(int resolution)
    {
        if (perVert)
            UVBox=vcg::tri::UV_Utils<MeshType>::PerVertUVBox(param_mesh);
        else
            UVBox=vcg::tri::UV_Utils<MeshType>::PerWedgeUVBox(param_mesh);
        //inflate a bit for numerical issues
        UVBox.Offset(UVBox.Diag()/100);

        //then find number of fibres
        if (fibre_space.X()<0)
        {
            assert(resolution>1);
            assert(fibre_space.Y()<0);
            ScalarType unit_size=UVBox.Diag()/resolution;
            fibre_space=Point2Type(unit_size,unit_size);
        }
        assert(fibre_space.X()>0);
        assert(fibre_space.Y()>0);

        size_t NumUFibersU=ceil((UVBox.max.X()-UVBox.min.X())/fibre_space.X());
        size_t NumUFibersV=ceil((UVBox.max.Y()-UVBox.min.Y())/fibre_space.Y());

        //then resize
        UFibers.clear();
        UFibers.resize(NumUFibersU+1);
        VFibers.clear();
        VFibers.resize(NumUFibersV+1);
    }

    void UpdateElongation(Fiber &f)
    {
        f.Len3D=0;
        f.LenParam=0;
        for (size_t i=0;i<f.FiberI.size();i++)
        {
            CoordType Pos0=f.FiberI[i].Int3DP_0;
            CoordType Pos1=f.FiberI[i].Int3DP_1;
            Point2Type IntUV0=f.FiberI[i].IntUV0;
            Point2Type IntUV1=f.FiberI[i].IntUV1;
            f.Len3D+=(Pos1-Pos0).Norm();
            f.LenParam+=(IntUV0-IntUV1).Norm();
        }
        f.RatioL=f.Len3D/f.LenParam;
        f.DiffL=(f.Len3D-f.LenParam);
    }

    void UpdateElongation()
    {
        MaxElongR=0;
        MaxCompressR=std::numeric_limits<ScalarType>::max();
        MaxElongL=0;
        MaxCompressL=std::numeric_limits<ScalarType>::max();

        TotDiffL=0;

        for (size_t i=0;i<UFibers.size();i++)
        {
            UpdateElongation(UFibers[i]);
            MaxElongR=std::max(MaxElongR,UFibers[i].RatioL);
            MaxCompressR=std::min(MaxCompressR,UFibers[i].RatioL);

            MaxElongL=std::max(MaxElongL,UFibers[i].DiffL);
            MaxCompressL=std::min(MaxCompressL,UFibers[i].DiffL);

            TotDiffL+=fabs(VFibers[i].DiffL);
        }
        for (size_t i=0;i<VFibers.size();i++)
        {
            UpdateElongation(VFibers[i]);
            MaxElongR=std::max(MaxElongR,VFibers[i].RatioL);
            MaxCompressR=std::min(MaxCompressR,VFibers[i].RatioL);

            MaxElongL=std::max(MaxElongL,VFibers[i].DiffL);
            MaxCompressL=std::min(MaxCompressL,VFibers[i].DiffL);

            TotDiffL+=fabs(VFibers[i].DiffL);
        }
        if (WriteDebug)
        {
            std::cout<<"Max Elongation Ratio "<<MaxElongR<<std::endl;
            std::cout<<"Max Compression Ratio "<<MaxCompressR<<std::endl;
            std::cout<<"Max Elongation Lenght "<<MaxElongL<<std::endl;
            std::cout<<"Max Compression Lenght "<<MaxCompressL<<std::endl;
        }
    }

    void GetFiberMesh(Fiber &f,
                      MeshType &FibreMesh3D,
                      MeshType &FibreMeshUV,
                      ScalarType radius)
    {
        FibreMesh3D.Clear();
        FibreMeshUV.Clear();
        for (size_t i=0;i<f.FiberI.size();i++)
        {
            CoordType Pos0=f.FiberI[i].Int3DP_0;
            CoordType Pos1=f.FiberI[i].Int3DP_1;
            CoordType PosUV0=CoordType(f.FiberI[i].IntUV0.X(),
                                       f.FiberI[i].IntUV0.Y(),0);
            CoordType PosUV1=CoordType(f.FiberI[i].IntUV1.X(),
                                       f.FiberI[i].IntUV1.Y(),0);
            MeshType Seg3D,SegUV;
            vcg::tri::OrientedCylinder(Seg3D,Pos0,Pos1,radius,true,8,2);
            vcg::tri::OrientedCylinder(SegUV,PosUV0,PosUV1,radius,true,8,2);
            vcg::tri::Append<MeshType,MeshType>::Mesh(FibreMesh3D,Seg3D);
            vcg::tri::Append<MeshType,MeshType>::Mesh(FibreMeshUV,SegUV);
        }
    }

    void ExtractVectFiberMeshes(MeshType &Mesh3D,MeshType &MeshUV,
                                std::vector<Fiber> Fibers,
                                ScalarType MaxElongRatio=-1,
                                ScalarType MaxComprRatio=-1,
                                ScalarType radius=0.5)
    {
        assert(MaxCompressR<=1);
        assert(MaxElongR>=1);
        ScalarType diff_to_1_Compr=fabs(1-MaxCompressR);
        ScalarType diff_to_1_Elong=fabs(1-MaxElongR);
        ScalarType diff_to_1=std::max(diff_to_1_Compr,diff_to_1_Elong);

        if (MaxComprRatio==-1)
            MaxComprRatio=1-diff_to_1;
        if (MaxElongRatio==-1)
            MaxElongRatio=1+diff_to_1;

        assert(MaxElongRatio>=1);
        assert(MaxComprRatio<=1);

        //        std::cout<<"Min Interval Color "<<MaxComprRatio<<std::endl;
        //        std::cout<<"Max Interval Color "<<MaxElongRatio<<std::endl;

        Mesh3D.Clear();
        MeshUV.Clear();
        ScalarType abs_radius=fibre_space.Norm()*radius;
        for (size_t i=0;i<Fibers.size();i++)
        {
            MeshType FibreMesh3D,FibreMeshUV;
            GetFiberMesh(Fibers[i],FibreMesh3D,FibreMeshUV,abs_radius);
            ScalarType Fibre_RatioL=Fibers[i].RatioL;

            vcg::Color4b col;
            col=vcg::Color4b::ColorRamp(MaxComprRatio,MaxElongRatio,Fibre_RatioL);

            vcg::tri::UpdateColor<MeshType>::PerFaceConstant(FibreMesh3D,col);
            vcg::tri::UpdateColor<MeshType>::PerFaceConstant(FibreMeshUV,col);

            vcg::tri::Append<MeshType,MeshType>::Mesh(Mesh3D,FibreMesh3D);
            vcg::tri::Append<MeshType,MeshType>::Mesh(MeshUV,FibreMeshUV);
        }
    }

    ScalarType Area()
    {
        ScalarType A=0;
        for (size_t i=0;i<param_mesh.face.size();i++)
            A+=vcg::DoubleArea(param_mesh.face[i]);
        return A/2;
    }

public:

    LenghMode LMode;

    void ExtractFiberMeshes(MeshType &Mesh3D,MeshType &MeshUV,
                            ScalarType MaxElongRatio=-1,
                            ScalarType MaxComprRatio=-1,
                            ScalarType radius=0.01)
    {
        Mesh3D.Clear();
        MeshUV.Clear();
        vcg::tri::UpdateBounding<MeshType>::Box(param_mesh);
        ScalarType abs_radius=sqrt(Area())*radius;

        MeshType FibreMesh3D_U,FibreMesh3D_V;
        MeshType FibreMeshUV_U,FibreMeshUV_V;

        ExtractVectFiberMeshes(FibreMesh3D_U,FibreMeshUV_U,UFibers,
                               MaxElongRatio,MaxComprRatio,abs_radius);

        ExtractVectFiberMeshes(FibreMesh3D_V,FibreMeshUV_V,VFibers,
                               MaxElongRatio,MaxComprRatio,abs_radius);

        vcg::tri::Append<MeshType,MeshType>::Mesh(Mesh3D,FibreMesh3D_U);
        vcg::tri::Append<MeshType,MeshType>::Mesh(Mesh3D,FibreMesh3D_V);

        vcg::tri::Append<MeshType,MeshType>::Mesh(MeshUV,FibreMeshUV_U);
        vcg::tri::Append<MeshType,MeshType>::Mesh(MeshUV,FibreMeshUV_V);
    }

    void ComputeFibreElongation(int resolution=100,
                                Point2Type _fibre_space=Point2Type(-1,-1),
                                bool _perVert=true)
    {
        perVert=_perVert;
        fibre_space=_fibre_space;
        InitFibres(resolution);
        UpdateFibreIntersection();
        UpdateElongation();

    }


    void FibreOptimization(int resolution=100,
                           size_t steps=10,
                           ScalarType Damp=0.5)
    {
        vcg::tri::UpdateTopology<MeshType>::FaceFace(param_mesh);
        vcg::tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(param_mesh);
        InitMaxDispl();
        ComputeFibreElongation(resolution);
        //InitVectDirection();
        ScalarType OldEnergy=TotDiffL;//std::max(MaxElongR,((ScalarType)1)/MaxCompressR);
        std::cout<<"Energy: "<<OldEnergy<<std::endl;
        for (size_t s=0;s<steps;s++)
        {
            //save old positions
            std::vector<Point2Type> CurrUV;
            for (size_t i=0;i<param_mesh.vert.size();i++)
                CurrUV.push_back(param_mesh.vert[i].T().P());

            OptimizationStep(Damp);
            ComputeFibreElongation(resolution);
            //InitVectDirection();
            ScalarType NewEnergy=TotDiffL;//std::max(MaxElongR,((ScalarType)1)/MaxCompressR);
            std::cout<<"Step: "<<s<<std::endl;
            std::cout<<"Energy: "<<NewEnergy<<std::endl;
            if (NewEnergy>OldEnergy)
            {
                for (size_t i=0;i<param_mesh.vert.size();i++)
                    param_mesh.vert[i].T().P()=CurrUV[i];

                break;
            }
            OldEnergy=NewEnergy;
        }
    }

    FibreElongation(MeshType &_param_mesh):param_mesh(_param_mesh)
    {
        LMode=LModeGLobal;
        WriteDebug=false;
        MaxScale=2;
    }
};
