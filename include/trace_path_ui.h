#ifndef PATH_UI
#define PATH_UI

#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/algorithms/closest.h>
#include "vcg/complex/algorithms/geodesic.h"

template <class TriMeshType>
class PathUI
{
    typedef typename TriMeshType::FaceType FaceType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;
    typedef std::pair<size_t,size_t> EdgeFace;

    TriMeshType &mesh;

    //std::vector<std::vector<CoordType> > PickedPoints;
    //std::vector<std::vector<vcg::Point2i> > Pixels;
    vcg::GridStaticPtr<VertexType,ScalarType> Grid;
    //std::vector<EdgeFace> Snapped;
    std::vector<std::vector<size_t> > VertexPath;

    void SmoothVert(std::set<std::pair<VertexType*,VertexType*> > &VertexPathSet,
                    const std::vector<bool> &OnPath,bool onlyLine,ScalarType Damp)
    {
        std::vector<CoordType> NewPos(mesh.vert.size(),CoordType(0,0,0));
        std::vector<size_t> Sum(mesh.vert.size(),0);

        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                VertexType* v0=mesh.face[i].V0(j);
                VertexType* v1=mesh.face[i].V1(j);
                size_t IndexV0=vcg::tri::Index(mesh,v0);
                size_t IndexV1=vcg::tri::Index(mesh,v1);
                std::pair<VertexType*,VertexType*> key(std::min(v0,v1),std::max(v0,v1));
                if ((VertexPathSet.count(key)==0)&&(onlyLine))continue;

                NewPos[IndexV0]+=v1->P();
                NewPos[IndexV1]+=v0->P();
                Sum[IndexV0]++;
                Sum[IndexV1]++;
            }
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (Sum[i]==0)continue;
            NewPos[i]/=Sum[i];
            if ((onlyLine)&&(!OnPath[i]))continue;
            if ((onlyLine)&&(Sum[i]==2))continue;//extreme
            if ((!onlyLine)&&(OnPath[i]))continue;
            if (mesh.vert[i].IsB())continue;
            if (!mesh.vert[i].IsS())continue;

            mesh.vert[i].P()=mesh.vert[i].P()*Damp+NewPos[i]*(1-Damp);
        }
    }

    void Reproject(vcg::GridStaticPtr<FaceType,ScalarType> &GridF)
    {
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            ScalarType maxD=mesh.bbox.Diag();
            ScalarType minD;
            CoordType closestPt;
            if (mesh.vert[i].IsB())continue;
            vcg::tri::GetClosestFaceBase(mesh,GridF,mesh.vert[i].P(),maxD,minD,closestPt);
            mesh.vert[i].P()=closestPt;
        }
    }

    void SnapPathOnVertices(const std::vector<std::vector<CoordType> > &PickedPoints)
    {
        //VertexPath.resize(VertexPath.size()+1);
        VertexPath.clear();
        VertexPath.resize(PickedPoints.size());
        ScalarType maxD=mesh.bbox.Diag();
        for (size_t i=0;i<PickedPoints.size();i++)
        {
            for (size_t j=0;j<PickedPoints[i].size();j++)
            {
                ScalarType minD;
                VertexType *v=vcg::tri::GetClosestVertex(mesh,Grid,PickedPoints[i][j],maxD,minD);
                size_t indexV=vcg::tri::Index(mesh,v);
                assert(v!=NULL);
                if ((VertexPath[i].size()>0)
                        &&(indexV==VertexPath[i].back()))continue;
                VertexPath[i].push_back(indexV);
            }
        }

    }


    void SmoothPaths(size_t steps=3)
    {
        std::set<std::pair<VertexType*,VertexType*> > VertexPathSet;
        vcg::tri::UpdateSelection<TriMeshType>::Clear(mesh);
        std::vector<bool> OnPath(mesh.vert.size(),false);
        for (size_t i=0;i<VertexPath.size();i++)
        {
            for (size_t j=0;j<VertexPath[i].size()-1;j++)
            {
                VertexType* v0=&mesh.vert[VertexPath[i][j]];
                VertexType* v1=&mesh.vert[VertexPath[i][(j+1)]];
                std::pair<VertexType*,VertexType*> key(std::min(v0,v1),std::max(v0,v1));
                VertexPathSet.insert(key);
                size_t IndexV0=vcg::tri::Index(mesh,v0);
                size_t IndexV1=vcg::tri::Index(mesh,v1);
                OnPath[IndexV0]=true;
                OnPath[IndexV1]=true;
                mesh.vert[IndexV0].SetS();
                mesh.vert[IndexV1].SetS();
            }
        }

        vcg::GridStaticPtr<FaceType,ScalarType> GridF;
        GridF.Set(mesh.face.begin(),mesh.face.end());

        int dilate_step=1;
        for (int i=0;i<dilate_step;i++)
        {
            vcg::tri::UpdateSelection<TriMeshType>::FaceFromVertexLoose(mesh);
            vcg::tri::UpdateSelection<TriMeshType>::VertexFromFaceLoose(mesh);
        }

        //first smooth line
        for (size_t s=0;s<steps;s++)
        {
            //first smooth line
            SmoothVert(VertexPathSet,OnPath,true,0.5);
            //then the rest
            SmoothVert(VertexPathSet,OnPath,false,0.5);

            Reproject(GridF);
        }
        //UpdateMesh();
    }

    void GetPath(VertexType *v0,
                 VertexType *v1,
                 std::vector<size_t> &Path)
    {
        Path.clear();

        typename TriMeshType::template PerVertexAttributeHandle<VertexType*> father;
        father=tri::Allocator<TriMeshType>::template GetPerVertexAttribute<VertexType*> (mesh,"father");

        std::vector<VertexType*> seedVec;
        seedVec.push_back(v0);
        vcg::tri::Mark(mesh,v0);
        vcg::tri::EuclideanDistance<TriMeshType> ED;
        ScalarType maxDistanceThr=mesh.bbox.Diag();
        assert(v0!=v1);
        vcg::tri::Geodesic<TriMeshType>::PerVertexDijkstraCompute(mesh,seedVec,ED,maxDistanceThr,
                                                                  NULL,NULL,&father,false,v1);
        assert(father[v1]!=NULL);
        VertexType *currV=v1;
        do{
            size_t currI=vcg::tri::Index(mesh,currV);
            Path.push_back(currI);
            currV=father[currV];
        }while (currV!=v0);
        Path.push_back(vcg::tri::Index(mesh,v0));
        std::reverse(Path.begin(),Path.end());
    }

    void CompletePath()
    {

        //create a map of edges
        std::set<std::pair<size_t,size_t> > EdgeMap;
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                assert(IndexV0!=IndexV1);
                EdgeMap.insert(std::pair<size_t,size_t>(std::min(IndexV0,IndexV1),
                                                        std::max(IndexV0,IndexV1)));

            }

        for (size_t i=0;i<VertexPath.size();i++)
        {
            std::vector<size_t> NewSeq;
            for (size_t j=0;j<VertexPath[i].size()-1;j++)
            {
                NewSeq.push_back(VertexPath[i][j]);
                size_t IndexV0=VertexPath[i][j];
                size_t IndexV1=VertexPath[i][j+1];
                std::pair<size_t,size_t> Key(std::min(IndexV0,IndexV1),
                                             std::max(IndexV0,IndexV1));

                if (EdgeMap.count(Key)>0)continue;

                std::vector<size_t> Path;
                VertexType *v0=&mesh.vert[IndexV0];
                VertexType *v1=&mesh.vert[IndexV1];
                GetPath(v0,v1,Path);

                for (size_t k=1;k<Path.size()-1;k++)
                    NewSeq.push_back(Path[k]);

            }
            NewSeq.push_back(VertexPath[i].back());
            VertexPath[i]=NewSeq;
        }
    }



public:


    void AddSharpConstraints(const std::vector<std::vector<CoordType> > &PickedPoints)
    {
        Grid.Set(mesh.vert.begin(),mesh.vert.end());

        mesh.SharpFeatures.clear();
        vcg::tri::UpdateSelection<TriMeshType>::Clear(mesh);
        SnapPathOnVertices(PickedPoints);

        CompletePath();
        SmoothPaths();

        std::set<std::pair<size_t,size_t> > EdgeSet;
        for (size_t i=0;i<VertexPath.size();i++)
            for (size_t j=0;j<VertexPath[i].size()-1;j++)
            {
                size_t IndexV0=VertexPath[i][j];
                size_t IndexV1=VertexPath[i][j+1];
                std::pair<size_t,size_t> Key(std::min(IndexV0,IndexV1),
                                             std::max(IndexV0,IndexV1));
                EdgeSet.insert(Key);
            }
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                VertexType* v0=mesh.face[i].V0(j);
                VertexType* v1=mesh.face[i].V1(j);
                size_t IndexV0=vcg::tri::Index(mesh,v0);
                size_t IndexV1=vcg::tri::Index(mesh,v1);
                std::pair<size_t,size_t> Key(std::min(IndexV0,IndexV1),
                                             std::max(IndexV0,IndexV1));
                if (EdgeSet.count(Key)==0)continue;
                mesh.face[i].SetFaceEdgeS(j);
                mesh.SharpFeatures.push_back(std::pair<size_t,size_t>(i,j));
            }
        mesh.InitRPos();
        //SmoothPaths();
    }

    PathUI(TriMeshType &_mesh):mesh(_mesh){}

};

#endif
