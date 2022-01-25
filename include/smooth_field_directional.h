/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2014                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef SMOOTHER_FIELD_DIRECTIONAL_H
#define SMOOTHER_FIELD_DIRECTIONAL_H

//vcg stuff
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <vcg/complex/algorithms/mesh_to_matrix.h>

//igl related stuff
#include <igl/principal_curvature.h>

//directional stuff
#include <directional/polyvector_field.h>
#include <directional/polyvector_to_raw.h>


//namespace vcg {
//namespace tri {

template <class MeshType>
class DirectionalFieldSmoother
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType CoordType;

public:

public:

    struct SmoothParam
    {
        //the 90° rotation independence while smoothing the direction field
        int Ndir;
        //the weight of curvature if doing the smoothing keeping the field close to the original one
        ScalarType alpha_curv;
        //align the field to border or not
        bool align_borders;
        //threshold to consider some edge as sharp feature and to use as hard constraint (0, not use)
        ScalarType sharp_thr;
        //threshold to consider some edge as high curvature anisotropyand to use as hard constraint (0, not use)
        ScalarType curv_thr;
        //the number of faces of the ring used ot esteem the curvature
        int curvRing;
        //this are additional hard constraints
        std::vector<std::pair<int,CoordType> > AddConstr;
        //use the predefined field and anisotropy or not
        bool use_predefined_field;

        SmoothParam()
        {
            Ndir=4;
            curvRing=3;
            alpha_curv=1000;
            align_borders=true;
            //SmoothM=SMMiq;
            sharp_thr=0.0;
            curv_thr=0.4;
            use_predefined_field=false;
            //IteN=20;
        }

    };

private:

    static void InitQualityByAnisotropyDir(MeshType &mesh)
    {
        ScalarType minV=0.00001;
        std::vector<ScalarType> QVal;
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            ScalarType N1=fabs(mesh.vert[i].K1());
            ScalarType N2=fabs(mesh.vert[i].K2());
            QVal.push_back(N1);
            QVal.push_back(N2);
        }

        std::sort(QVal.begin(),QVal.end());
        int percUp=int(floor(QVal.size()*0.9+0.5));
        ScalarType trimUp=QVal[percUp];

        for (size_t i=0;i<mesh.vert.size();i++)
        {
            ScalarType N1=(mesh.vert[i].K1());
            ScalarType N2=(mesh.vert[i].K2());

            ScalarType NMax=std::max(N1,N2)/trimUp;
            ScalarType NMin=std::min(N1,N2)/trimUp;

            if (NMax>1)NMax=1;
            if (NMax<-1)NMax=-1;

            if (NMin>1)NMin=1;
            if (NMin<-1)NMin=-1;

            ScalarType CurvAni=(NMax-NMin)/2;
            CurvAni=std::max(CurvAni,minV);
            mesh.vert[i].Q()=CurvAni;
        }
        vcg::tri::UpdateQuality<MeshType>::FaceFromVertex(mesh);
    }

    static void SetEdgeDirection(FaceType *f,int edge)
    {
        CoordType dir=f->P0(edge)-f->P1(edge);
        dir.Normalize();
        ScalarType prod1=fabs(dir*f->PD1());
        ScalarType prod2=fabs(dir*f->PD2());
        if (prod1>prod2)
        {
            f->PD1()=dir;
            f->PD2()=f->N()^dir;
        }else
        {
            f->PD2()=dir;
            f->PD1()=f->N()^dir;
        }
    }

    static void AddSharpEdgesConstraints(MeshType & mesh,
                                         const ScalarType &thr,
                                         std::vector<size_t> &IndexF)
    {
        for (size_t i=0;i<mesh.face.size();i++)
            for (int j=0;j<mesh.face[i].VN();j++)
            {
                FaceType *f0=&mesh.face[i];
                FaceType *f1=f0->FFp(j);
                if (f0==f1)continue;
                CoordType N0=f0->N();
                CoordType N1=f1->N();
                if ((N0*N1)>thr)continue;
                SetEdgeDirection(f0,j);
                IndexF.push_back(i);
                //f0->SetS();
            }
    }

    static void AddBorderConstraints(MeshType & mesh,
                                     std::vector<size_t> &IndexF)
    {
        for (size_t i=0;i<mesh.face.size();i++)
            for (int j=0;j<mesh.face[i].VN();j++)
            {
                FaceType *f0=&mesh.face[i];
                FaceType *f1=f0->FFp(j);
                assert(f1!=NULL);
                if (f0!=f1)continue;
                SetEdgeDirection(f0,j);
                IndexF.push_back(i);
                //f0->SetS();
            }
    }

    static void AddCurvatureConstraints(MeshType & mesh,
                                        const ScalarType &thr,
                                        std::vector<size_t> &IndexF)
    {
        for (size_t i=0;i<mesh.face.size();i++)
            if (mesh.face[i].Q()>thr)IndexF.push_back(i);
    }

    static void SmoothNPoly(MeshType &mesh,
                            const Eigen::VectorXi &FaceI,   //target faces
                            const Eigen::MatrixXd &FaceD,   //target directions
                            const Eigen::VectorXd &alignWeights, //target weights (-1 -> fixed)
                            double smoothWeight,
                            int Ndir)
    {
        assert((Ndir==2)||(Ndir==4));

        Eigen::MatrixXi F;
        typename vcg::tri::MeshToMatrix<MeshType>::MatrixXm Vf;

        vcg::tri::MeshToMatrix<MeshType>::GetTriMeshData(mesh,F,Vf);
        //then cast
        Eigen::MatrixXd V = Vf.template cast<double>();
        Eigen::MatrixXcd output_field;

        double roSyWeight=-1;
        directional::polyvector_field(V, F, FaceI, FaceD, smoothWeight, roSyWeight, alignWeights, Ndir, output_field);

        Eigen::MatrixXd rawOutField;
        Eigen::MatrixXd normals, B1, B2;
        igl::local_basis(Vf, F, B1, B2, normals); // Compute a local orthogonal reference system for each triangle in the given mesh
        directional::polyvector_to_raw(B1, B2, output_field, Ndir, rawOutField);

        assert(output_field.rows()==mesh.face.size());
        //finally update the principal directions
        for (size_t i=0;i<mesh.face.size();i++)
        {
            CoordType dir1;
            dir1[0]=rawOutField(i,0);
            dir1[1]=rawOutField(i,1);
            dir1[2]=rawOutField(i,2);

            dir1.Normalize();
            CoordType dir2=mesh.face[i].N()^dir1;
            dir2.Normalize();

            ScalarType Norm1=mesh.face[i].PD1().Norm();
            ScalarType Norm2=mesh.face[i].PD2().Norm();

            mesh.face[i].PD1()=dir1*Norm1;
            mesh.face[i].PD2()=dir2*Norm2;
        }
    }

    static void PickRandomDir(const MeshType &mesh,
                              int indexF,
                              CoordType &Dir)
    {
        indexF=rand()%mesh.fn;
        const FaceType *currF=&mesh.face[indexF];
        CoordType dirN=currF->cN();
        dirN.Normalize();
        Dir=CoordType(1,0,0);
        if (fabs(Dir*dirN)>0.9)
            Dir=CoordType(0,1,0);
        if (fabs(Dir*dirN)>0.9)
            Dir=CoordType(0,0,1);

        Dir=dirN^Dir;
        Dir.Normalize();
    }

    static void CollectConstraintsData(const MeshType &mesh,
                                       const bool useCurvatureSoft,
                                       const std::vector<size_t> &FixedI,
                                       Eigen::VectorXi &FaceI,   //target faces
                                       Eigen::MatrixXd &FaceD,   //target directions
                                       Eigen::VectorXd &alignWeights)//,SmoothParam &SParam)
    {
        if (useCurvatureSoft)
        {
            //in this case add one row for each face
            int sizeV=mesh.face.size();
            FaceI=Eigen::VectorXi(sizeV);
            FaceD=Eigen::MatrixXd(sizeV,3);
            alignWeights=Eigen::VectorXd(sizeV);
            for (size_t i=0;i<mesh.face.size();i++)
            {
                FaceI(i)=i;
                alignWeights(i)=mesh.face[i].Q();
                FaceD(i,0)=mesh.face[i].PD1().X();
                FaceD(i,1)=mesh.face[i].PD1().Y();
                FaceD(i,2)=mesh.face[i].PD1().Z();
            }
            //then set the weight to -1 if constrained
            for (size_t i=0;i<FixedI.size();i++)
            {
                assert(FixedI[i]<mesh.face.size());
                alignWeights(FixedI[i])=-1;
            }
        }
        else
        {
            //otherwise just add the row for selected
            int sizeV=FixedI.size();
            FaceI=Eigen::VectorXi(sizeV);
            FaceD=Eigen::MatrixXd(sizeV,3);
            //set the weights to -1, all fixed for those
            alignWeights=Eigen::VectorXd(sizeV);
            for (size_t i=0;i<FixedI.size();i++)
            {
                size_t IndexF=FixedI[i];
                FaceI(i)=IndexF;
                alignWeights(i)=-1;
                FaceD(i,0)=mesh.face[IndexF].PD1().X();
                FaceD(i,1)=mesh.face[IndexF].PD1().Y();
                FaceD(i,2)=mesh.face[IndexF].PD1().Z();
            }
        }
    }


    static void GetConstrainedFaces(MeshType &mesh,
                                    SmoothParam &SParam,
                                    std::vector<size_t> &IndexF)
    {
        //clear all selected faces
        IndexF.clear();

        if (SParam.curv_thr>0)
            AddCurvatureConstraints(mesh,SParam.curv_thr,IndexF);///Ratio);

        //add alignment to sharp features
        if (SParam.sharp_thr>0)
            AddSharpEdgesConstraints(mesh,SParam.sharp_thr,IndexF);

        //add border constraints
        if (SParam.align_borders)
            AddBorderConstraints(mesh,IndexF);

        //additional user-defined constraints
        for (size_t i=0;i<SParam.AddConstr.size();i++)
        {
            int currI=SParam.AddConstr[i].first;
            CoordType dir=SParam.AddConstr[i].second;
            mesh.face[currI].PD1()=dir;
            mesh.face[currI].PD2()=mesh.face[currI].N()^dir;
            mesh.face[currI].PD1().Normalize();
            mesh.face[currI].PD2().Normalize();
            IndexF.push_back(currI);
        }

        std::sort(IndexF.begin(), IndexF.end());
        std::vector<size_t>::iterator last = std::unique(IndexF.begin(), IndexF.end());
        IndexF.erase(last, IndexF.end());
    }


    static void ColorByFieldGuide(MeshType &mesh,
                                  bool colorbyQ,
                                  const Eigen::VectorXi &FaceI,
                                  const Eigen::VectorXd &alignWeights)
    {
        if (colorbyQ)
            vcg::tri::UpdateColor<MeshType>::PerFaceQualityGray(mesh);
        else
            vcg::tri::UpdateColor<MeshType>::PerFaceConstant(mesh);

        for (size_t i=0;i<FaceI.rows();i++)
        {
            if (alignWeights(i)!=-1)continue;
            mesh.face[FaceI[i]].C()=vcg::Color4b::Red;
        }
    }

public:

    static void InitByCurvature(MeshType & mesh,
                                unsigned Nring,
                                bool UpdateFaces=true)
    {

        vcg::tri::RequirePerVertexCurvatureDir(mesh);

        Eigen::MatrixXi F;
        typename vcg::tri::MeshToMatrix<MeshType>::MatrixXm Vf;

        Eigen::MatrixXd PD1,PD2,PV1,PV2;
        vcg::tri::MeshToMatrix<MeshType>::GetTriMeshData(mesh,F,Vf);
        Eigen::MatrixXd V = Vf.template cast<double>();

        igl::principal_curvature(V,F,PD1,PD2,PV1,PV2,Nring,true);

        //then copy curvature per vertex
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            mesh.vert[i].PD1()=CoordType(PD1(i,0),PD1(i,1),PD1(i,2));
            mesh.vert[i].PD2()=CoordType(PD2(i,0),PD2(i,1),PD2(i,2));
            mesh.vert[i].K1()=PV1(i,0);
            mesh.vert[i].K2()=PV2(i,0);
        }
        if (!UpdateFaces)return;
        vcg::tri::CrossField<MeshType>::SetFaceCrossVectorFromVert(mesh);
        InitQualityByAnisotropyDir(mesh);
    }


    static void SmoothDirections(MeshType &mesh,
                                 SmoothParam SParam,
                                 bool color_mesh=false)
    {

//        //initialize curvature if needed
//        if ((SParam.alpha_curv>0)||
//            (SParam.curv_thr>0))
//        {
        if (!SParam.use_predefined_field)
            InitByCurvature(mesh,SParam.curvRing);
//        }

        //then select the constrained faces
        std::vector<size_t> FixedI;
        GetConstrainedFaces(mesh,SParam,FixedI);

        bool useCurvatureSoft=(SParam.alpha_curv>0);

        Eigen::VectorXi FaceI;
        Eigen::MatrixXd FaceD;
        Eigen::VectorXd alignWeights;
        CollectConstraintsData(mesh,useCurvatureSoft,FixedI,FaceI,FaceD,alignWeights);

        //check if need to add a random constraint
        if (FaceI.rows()==0)
        {
            FaceI=Eigen::VectorXi(1);
            FaceD=Eigen::MatrixXd(1,3);
            alignWeights=Eigen::VectorXd(1);

            FaceI(0)=0;
            alignWeights(0)=-1;

            CoordType Dir;
            PickRandomDir(mesh,0,Dir);
            FaceD(0,0)=Dir.X();
            FaceD(0,1)=Dir.Y();
            FaceD(0,2)=Dir.Z();
        }

        ScalarType currSm=1;
        if (SParam.alpha_curv>0)
            currSm=SParam.alpha_curv;
        SmoothNPoly(mesh,FaceI,FaceD,alignWeights,currSm, SParam.Ndir);

        if (color_mesh)
        {
            bool ColorByQ=(SParam.alpha_curv>0);
            ColorByFieldGuide(mesh,ColorByQ,FaceI,alignWeights);
        }
    }
};

//} // end namespace tri
//} // end namespace vcg
#endif // SMOOTHER_FIELD_DIRECTIONAL_H
