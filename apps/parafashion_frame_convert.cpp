#include <iostream>
#include <string>
#include <QDir>
#include <QDirIterator>

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/simplex/face/topology.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/algorithms/closest.h>

std::string pathFrames,pathRefMesh,pathRefGarnment;
std::vector<std::string> PathFrames;

class MyTriFace;
//class MyTriEdge;
class MyTriVertex;


struct TriUsedTypes: public vcg::UsedTypes<vcg::Use<MyTriVertex>::AsVertexType,
                                            //vcg::Use<MyTriEdge>::AsEdgeType,
                                            vcg::Use<MyTriFace>::AsFaceType>{};

class MyTriVertex:public vcg::Vertex<TriUsedTypes,
                                    vcg::vertex::Coord3f,
                                    vcg::vertex::Normal3f,
                                    vcg::vertex::VFAdj,
                                    vcg::vertex::BitFlags,
                                    vcg::vertex::Qualityf,
                                    vcg::vertex::Mark>
{};


class MyTriFace:public vcg::Face<TriUsedTypes,
                                vcg::face::VertexRef,
                                vcg::face::VFAdj,
                                vcg::face::FFAdj,
                                vcg::face::BitFlags,
                                vcg::face::Normal3f,
                                vcg::face::CurvatureDirf,
                                vcg::face::Qualityf,
                                vcg::face::Mark>
{};


class MyTriMesh: public vcg::tri::TriMesh< std::vector<MyTriVertex>,
                                           std::vector<MyTriFace > >
{

public:

    void UpdateDataStructures()
    {
        vcg::tri::UpdateNormal<MyTriMesh>::PerFaceNormalized(*this);
        vcg::tri::UpdateNormal<MyTriMesh>::PerVertexNormalized(*this);
        vcg::tri::UpdateTopology<MyTriMesh>::FaceFace(*this);
        vcg::tri::UpdateTopology<MeshType>::VertexFace(*this);
        vcg::tri::UpdateBounding<MyTriMesh>::Box(*this);
    }

    bool LoadTriMesh(const std::string &filename)
    {
        Clear();
        if(filename.empty()) return false;
        int position0=filename.find(".ply");
        int position1=filename.find(".obj");
        int position2=filename.find(".off");
        std::cout<<"Loading:"<<filename.c_str()<<std::endl;

        if (position0!=-1)
        {
            int err=vcg::tri::io::ImporterPLY<MyTriMesh>::Open(*this,filename.c_str());
            if (err!=vcg::ply::E_NOERROR)return false;
            UpdateDataStructures();
            return true;
        }
        if (position1!=-1)
        {
            int Mask=0;
            vcg::tri::io::ImporterOBJ<MyTriMesh>::LoadMask(filename.c_str(), Mask);
            int err=vcg::tri::io::ImporterOBJ<MyTriMesh>::Open(*this,filename.c_str(),Mask);
            std::cout<<vcg::tri::io::ImporterOBJ<MyTriMesh>::ErrorMsg(err)<<std::endl;
            if ((err!=0)&&(err!=5))return false;
            UpdateDataStructures();
            return true;
        }
        if (position2!=-1)
        {
            int err=vcg::tri::io::ImporterOFF<MyTriMesh>::Open(*this,filename.c_str());
            if (err!=0)return false;
            UpdateDataStructures();
            return true;
        }
        return false;
    }
};

typedef typename MyTriMesh::CoordType CoordType;
typedef typename MyTriMesh::ScalarType ScalarType;
std::vector<std::vector<CoordType> > PerFrameTemplatePos;
std::vector<std::vector<CoordType> > PerFrameTemplateNorm;


void GetAllObjInDir(std::string dir_path)
{
    QString QPath=(dir_path.c_str());
    QDir dir(QPath);
    if (!dir.exists()){
        std::cout<<"Wrong Directory"<<std::endl;
        exit(0);
    }
    QDirIterator it(QPath, QStringList() << "*.obj", QDir::Files);//, QDirIterator::NoIteratorFlags);
    while (it.hasNext())
    {
        std::string NameF;
        it.next();
        NameF=it.fileName().toStdString();
        PathFrames.push_back(NameF);
    }

    std::sort(PathFrames.begin(),PathFrames.end());

    //then load all frames
    PerFrameTemplatePos.clear();
    PerFrameTemplateNorm.clear();
    for (size_t i=0;i<PathFrames.size();i++)
    {
        PerFrameTemplatePos.resize(PerFrameTemplatePos.size()+1);
        PerFrameTemplateNorm.resize(PerFrameTemplateNorm.size()+1);

        MyTriMesh m;
        std::string FinalPath=dir_path+PathFrames[i];
        bool Loaded=m.LoadTriMesh(FinalPath.c_str());
        assert(Loaded);

        for (size_t i=0;i<m.vert.size();i++)
        {
            CoordType Pos=m.vert[i].P();
            CoordType Norm=m.vert[i].N();
            PerFrameTemplatePos.back().push_back(Pos);
            PerFrameTemplateNorm.back().push_back(Norm);
        }
    }
}

std::vector<size_t> VertToFaceMap;
std::vector<CoordType> VertToFaceBary;
std::vector<CoordType> VertInterpN;
std::vector<CoordType> VertCloseP;

MyTriMesh m_template,m_garnment;

void InitMapping(std::string dir_template_path,
                 std::string dir_garnment_path)
{
    m_template.Clear();
    m_garnment.Clear();

    bool loaded_template=m_template.LoadTriMesh(dir_template_path.c_str());
    if (!loaded_template)
    {
        std::cout<<"Error Loading Template Mesh"<<std::endl;
        exit(0);
    }
    bool loaded_garnment=m_garnment.LoadTriMesh(dir_garnment_path);
    if (!loaded_garnment)
    {
        std::cout<<"Error Loading Garnment Mesh"<<std::endl;
        exit(0);
    }

    vcg::GridStaticPtr<MyTriFace,ScalarType> Grid;
    Grid.Set(m_template.face.begin(),m_template.face.end());
    VertToFaceMap.clear();
    VertToFaceBary.clear();
    VertInterpN.clear();
    VertCloseP.clear();

    ScalarType MaxD=m_garnment.bbox.Diag();
    for (size_t i=0;i<m_garnment.vert.size();i++)
    {
        CoordType Pos=m_garnment.vert[i].P();
        CoordType closestPt,normF,Bary;

        ScalarType MinD;
        MyTriFace *f=vcg::tri::GetClosestFaceBase(m_template,Grid,Pos,MaxD,MinD,closestPt,normF,Bary);
        assert(f!=NULL);

        size_t IndexF=vcg::tri::Index(m_template,f);
        VertToFaceMap.push_back(IndexF);
        VertToFaceBary.push_back(Bary);
        VertInterpN.push_back(normF);
        VertCloseP.push_back(Pos);
    }
}


std::vector<std::vector<CoordType> > PerFrameGarnmentPos;

void ComputeGarnmentPosition()
{
    assert(VertToFaceMap.size()==m_garnment.vert.size());
    assert(VertToFaceBary.size()==m_garnment.vert.size());

    //for each frame
    PerFrameGarnmentPos.clear();
    for (size_t i=0;i<PerFrameTemplatePos.size();i++)
    {
        //for each vertex
        PerFrameGarnmentPos.resize(PerFrameGarnmentPos.size()+1);
        for (size_t j=0;j<m_garnment.vert.size();j++)
        {
            int IndexF=VertToFaceMap[j];
            CoordType Bary=VertToFaceBary[j];
            CoordType Nrest=VertInterpN[j];
            CoordType Clos=VertCloseP[j];
            CoordType DiffVect=m_garnment.vert[j].P()-Clos;

            size_t IndexV0Templ=vcg::tri::Index(m_template,m_template.face[IndexF].V(0));
            size_t IndexV1Templ=vcg::tri::Index(m_template,m_template.face[IndexF].V(1));
            size_t IndexV2Templ=vcg::tri::Index(m_template,m_template.face[IndexF].V(2));

            CoordType P0=PerFrameTemplatePos[i][IndexV0Templ];
            CoordType P1=PerFrameTemplatePos[i][IndexV1Templ];
            CoordType P2=PerFrameTemplatePos[i][IndexV2Templ];

            CoordType N0=PerFrameTemplateNorm[i][IndexV0Templ];
            CoordType N1=PerFrameTemplateNorm[i][IndexV1Templ];
            CoordType N2=PerFrameTemplateNorm[i][IndexV2Templ];
            CoordType InterpNorm=N0*Bary.X()+N1*Bary.Y()+N2*Bary.Z();
            vcg::Matrix33<ScalarType> RotM=vcg::RotationMatrix<ScalarType>(Nrest,InterpNorm);

            DiffVect=RotM*DiffVect;
            CoordType InterpPos=P0*Bary.X()+P1*Bary.Y()+P2*Bary.Z();
            //InterpPos+=DiffVect;
            PerFrameGarnmentPos.back().push_back(InterpPos);
        }
    }
}

void SetGarnmentFrame(size_t IndexFrame)
{
    assert(IndexFrame>=0);
    assert(IndexFrame<PerFrameGarnmentPos.size());
    assert(PerFrameGarnmentPos[IndexFrame].size()==m_garnment.vert.size());
    for (size_t i=0;i<m_garnment.vert.size();i++)
        m_garnment.vert[i].P()=PerFrameGarnmentPos[IndexFrame][i];
}

void SaveFrames()
{
    for (size_t i=0;i<PerFrameGarnmentPos.size();i++)
    {
        char buffer [50];
        sprintf (buffer, "%d.ply",(int)i);
        SetGarnmentFrame(i);
        vcg::tri::io::ExporterPLY<MyTriMesh>::Save(m_garnment,buffer);
    }
}

void SaveGarnmentPosFrames(const std::string &path)
{
    FILE *f=fopen(path.c_str(),"wt");
    fprintf(f,"%d\n",(int)PerFrameGarnmentPos.size());
    for (size_t i=0;i<PerFrameGarnmentPos.size();i++)
        for (size_t j=0;j<PerFrameGarnmentPos[i].size();j++)
        {
            fprintf(f,"%f,%f,%f\n",PerFrameGarnmentPos[i][j].X(),
                                   PerFrameGarnmentPos[i][j].Y(),
                                   PerFrameGarnmentPos[i][j].Z());

        }

    fclose(f);
}

int main(int argc, char *argv[])
{
    if (argc<4)
    {
        std::cout<<"Arguments: Template Mesh, Template Garnment, Animated Template"<<std::endl;
        exit(0);
    }

    std::string path=std::string(std::string(argv[3]));

    std::cout<<"Loading Frames"<<std::endl;
    GetAllObjInDir(path);

    std::cout<<"Init Mapping"<<std::endl;
    InitMapping(std::string(argv[1]),std::string(argv[2]));

    std::cout<<"Transferring Garnment Pos"<<std::endl;
    ComputeGarnmentPosition();

    std::cout<<"Saving output pos"<<std::endl;
    SaveGarnmentPosFrames("out_frames.txt");

//    std::cout<<"Saving Debug garnment"<<std::endl;
//    SaveFrames();
    std::cout<<"Done"<<std::endl;

}
