#ifndef MY_TRI_MESH_TYPE
#define MY_TRI_MESH_TYPE

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/simplex/face/topology.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>
#include <wrap/gl/trimesh.h>

class MyTriFace;
class MyTriVertex;


struct TriUsedTypes: public vcg::UsedTypes<vcg::Use<MyTriVertex>::AsVertexType,
        vcg::Use<MyTriFace>::AsFaceType>{};

class MyTriVertex:public vcg::Vertex<TriUsedTypes,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d,
        vcg::vertex::VFAdj,
        vcg::vertex::BitFlags,
        vcg::vertex::TexCoord2d,
        //vcg::vertex::Curvatured,
        //vcg::vertex::CurvatureDird,
        vcg::vertex::Qualityd,
        vcg::vertex::Color4b
        //vcg::vertex::Mark
        >
{
public:
    CoordType RestPos;
};

class MyTriFace:public vcg::Face<TriUsedTypes,
        vcg::face::VertexRef,
        vcg::face::VFAdj,
        vcg::face::FFAdj,
        vcg::face::BitFlags,
        vcg::face::Normal3d,
        vcg::face::Color4b,
        vcg::face::Qualityd,
        //vcg::face::CurvatureDird,
        //
        //vcg::face::WedgeTexCoord2d,
        vcg::face::Mark
        >
{

};



class MyTriMesh: public vcg::tri::TriMesh< std::vector<MyTriVertex>,
        std::vector<MyTriFace > >
{

public:
    //VCG UPDATING STRUCTURES
    void UpdateDataStructures()
    {
        vcg::tri::UpdateBounding<MyTriMesh>::Box(*this);
        vcg::tri::UpdateNormal<MyTriMesh>::PerVertexNormalizedPerFace(*this);
        vcg::tri::UpdateNormal<MyTriMesh>::PerFaceNormalized(*this);
        vcg::tri::UpdateTopology<MyTriMesh>::FaceFace(*this);
        vcg::tri::UpdateTopology<MyTriMesh>::VertexFace(*this);
    }

    void InitRestPos()
    {
        for (size_t i=0;i<vert.size();i++)
            vert[i].RestPos=vert[i].P();
    }

    void RestoreRestPos()
    {
        for (size_t i=0;i<vert.size();i++)
            vert[i].P()=vert[i].RestPos;
    }

    void SetUVPosition()
    {
        for (size_t i=0;i<vert.size();i++)
            vert[i].P()=CoordType(vert[i].T().P().X(),vert[i].T().P().Y(),0);

        vcg::tri::UpdateBounding<MyTriMesh>::Box(*this);
        vcg::tri::UpdateNormal<MyTriMesh>::PerFaceNormalized(*this);
        vcg::tri::UpdateNormal<MyTriMesh>::PerVertexNormalizedPerFace(*this);
    }

    bool LoadTriMesh(const std::string &filename)
    {
        Clear();
        if(filename.empty()) return false;
        int position0=filename.find(".ply");
        int position1=filename.find(".obj");
        int position2=filename.find(".off");

        if (position0!=-1)
        {
            int err=vcg::tri::io::ImporterPLY<MyTriMesh>::Open(*this,filename.c_str());
            if (err!=vcg::ply::E_NOERROR)return false;
            InitRestPos();
            return true;
        }
        if (position1!=-1)
        {
            int mask;
            vcg::tri::io::ImporterOBJ<MyTriMesh>::LoadMask(filename.c_str(),mask);
            int err=vcg::tri::io::ImporterOBJ<MyTriMesh>::Open(*this,filename.c_str(),mask);
            if ((err!=0)&&(err!=5))return false;
            InitRestPos();
            return true;
        }
        if (position2!=-1)
        {
            int err=vcg::tri::io::ImporterOFF<MyTriMesh>::Open(*this,filename.c_str());
            if (err!=0)return false;
            InitRestPos();
            return true;
        }
        return false;
    }



};


#endif
