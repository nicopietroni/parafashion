#ifndef PARAFASHION_INTERFACE
#define PARAFASHION_INTERFACE

#include "parafashion.h"

template <class MeshType>
void MeshToVectors(const MeshType &mesh,
                   std::vector<std::vector<size_t> > &Faces,
                   std::vector<std::vector<double> > &Coord,
                   std::vector<int> &Partition,
                   std::vector<std::vector<std::vector<double> > > &FaceVertUV)
{
    Faces.clear();
    Coord.clear();
    Partition.clear();
    FaceVertUV.clear();

    Coord.resize(mesh.vert.size(),std::vector<double>(3,0));
    FaceVertUV.resize(mesh.face.size());
    for (size_t i=0;i<mesh.vert.size();i++)
    {
        Coord[i][0]=mesh.vert[i].cP().X();
        Coord[i][1]=mesh.vert[i].cP().Y();
        Coord[i][2]=mesh.vert[i].cP().Z();
    }

    for (size_t i=0;i<mesh.face.size();i++)
    {
        FaceVertUV[i].clear();
        FaceVertUV[i].resize(3);
        for (size_t j=0;j<3;j++)
        {
            FaceVertUV[i][j].push_back(mesh.face[i].WT(j).P().X());
            FaceVertUV[i][j].push_back(mesh.face[i].WT(j).P().Y());
        }
    }

    Faces.resize(mesh.face.size());
    Partition.resize(mesh.face.size(),-1);
    for (size_t i=0;i<mesh.face.size();i++)
    {
        Partition[i]=mesh.face[i].Q();
        for (int j=0;j<mesh.face[i].VN();j++)
        {
            size_t IndexV=vcg::tri::Index(mesh,mesh.face[i].cV(j));
            Faces[i].push_back(IndexV);
        }
    }
}

template <class MeshType>
void VectorsToTriMesh(const std::vector<std::vector<size_t> > &Faces,
                      const std::vector<std::vector<double> > &Coord,
                      MeshType &mesh)
{
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;

    mesh.Clear();
    for (size_t i=0;i<Coord.size();i++)
    {
        CoordType pos(Coord[i][0],Coord[i][1],Coord[i][2]);
        vcg::tri::Allocator<MeshType>::AddVertex(mesh,pos);
    }

    vcg::tri::Allocator<MeshType>::AddFaces(mesh,Faces.size());
    for (size_t i=0;i<Faces.size();i++)
    {
        size_t NumV=Faces[i].size();
        assert(NumV==3);
        VertexType *v0=&mesh.vert[Faces[i][0]];
        VertexType *v1=&mesh.vert[Faces[i][1]];
        VertexType *v2=&mesh.vert[Faces[i][2]];
        assert(v0!=v1);
        assert(v1!=v2);
        assert(v2!=v0);
        mesh.face[i].V(0)=v0;
        mesh.face[i].V(1)=v1;
        mesh.face[i].V(2)=v2;
    }
    mesh.UpdateAttributes();
}

template <class MeshType>
void VectorsToTriMesh(const std::vector<std::vector<size_t> > &Faces,
                      const std::vector<std::vector<double> > &Coord,
                      std::vector<int> &Partition,
                      std::vector<std::vector<std::vector<double> > > &FaceVertUV,
                      MeshType &mesh)
{
    VectorsToTriMesh(Faces,Coord,mesh);
    assert(mesh.face.size()==FaceVertUV.size());
    for (size_t i=0;i<mesh.face.size();i++)
    {
        assert(FaceVertUV[i].size()==3);
        for (size_t j=0;j<3;j++)
        {
            assert(FaceVertUV[i][j].size()==2);
            mesh.face[i].WT(j).P().X()=FaceVertUV[i][j][0];
            mesh.face[i].WT(j).P().Y()=FaceVertUV[i][j][1];
        }
    }

    for (size_t i=0;i<mesh.face.size();i++)
    {
        mesh.face[i].Q()=Partition[i];
        mesh.face[i].C()=vcg::Color4b::Scatter(Partition.size(),Partition[i]);
    }

}

void DerivePatchLayout(const std::vector<std::vector<size_t> > &FacesReference,
                       const std::vector<std::vector<double> > &CoordReference,
                       const std::vector<std::vector<size_t> > &FacesDeformed,
                       const std::vector<std::vector<double> > &CoordDeformed,
                       const std::vector<std::vector<std::vector<double> > > &PickedPoints,
                       std::vector<std::vector<size_t> > &FacesOutput,
                       std::vector<std::vector<double> > &CoordOutput,
                       std::vector<int> &Partition,
                       std::vector<std::vector<std::vector<double> > > &FaceVertUV)
{
    typedef typename CMesh::CoordType CoordType;

    //TRANSFORM THE CONSTRAINTS
    std::vector<std::vector<CoordType> > ConstrPoints;
    ConstrPoints.resize(PickedPoints.size());
    for (size_t i=0;i<PickedPoints.size();i++)
    {
        ConstrPoints[i].resize(PickedPoints[i].size(),CoordType(0,0,0));
        for (size_t j=0;j<PickedPoints[i].size();j++)
        {
            assert(PickedPoints[i][j].size()==3);
            ConstrPoints[i][j].X()=PickedPoints[i][j][0];
            ConstrPoints[i][j].Y()=PickedPoints[i][j][1];
            ConstrPoints[i][j].Z()=PickedPoints[i][j][2];
        }
    }

    //TRANSFORM THE MESH
    assert(FacesReference.size()==FacesDeformed.size());
    assert(CoordReference.size()==CoordDeformed.size());

    CMesh deformed_mesh,reference_mesh;
    VectorsToTriMesh(FacesReference,CoordReference,reference_mesh);
    VectorsToTriMesh(FacesDeformed,CoordDeformed,deformed_mesh);

    //PROCESS THE DECOMPOSITION
    Parafashion<CMesh> PFashion(deformed_mesh,reference_mesh);
    PFashion.Init();
    PFashion.BatchProcess(ConstrPoints,false,false);


    Partition.clear();
    FaceVertUV.clear();

    //OUTPUT THE RESULTS
    MeshToVectors(deformed_mesh,FacesOutput,CoordOutput,Partition,FaceVertUV);
}

#endif
