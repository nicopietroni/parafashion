#ifndef SVG_EXPORTER
#define SVG_EXPORTER

#include <vcg/complex/complex.h>
#include <wrap/qt/Outline2ToQImage.h>
#include <vcg/simplex/face/pos.h>
#include <set>
#include <vcg/complex/algorithms/parametrization/uv_utils.h>
//#include "marching_squares.h"
#include <clipper.hpp>
#include <vcg/complex/algorithms/point_sampling.h>
#include <wrap/io_trimesh/export.h>
//#include "./lib/CavalierContours/include/cavc/polylineoffset.hpp"
//#include <vcg/complex/algorithms/update/flag.h>
#include <QSvgRenderer>
#include <QIcon>

template <class TriMeshType>
class SvgExporter
{
    typedef typename TriMeshType::FaceType FaceType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;

    static void GetOutLines(TriMeshType &mesh,
                            std::vector< std::vector<vcg::Point2f> > &outline2Vec,
                            std::vector< std::vector<CoordType> > &outline3Vec,
                            float scaleVal)
    {
        mesh.UpdateAttributes();
        //make the selection coherent
        for (size_t i=0;i<mesh.face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                if (vcg::face::IsBorder(mesh.face[i],j))
                {
                    mesh.face[i].SetFaceEdgeS(j);
                    continue;
                }
                if (!mesh.face[i].IsFaceEdgeS(j))continue;
                FaceType *FOpp=mesh.face[i].FFp(j);
                int EOpp=mesh.face[i].FFi(j);
                FOpp->SetFaceEdgeS(EOpp);
            }
        }
        //vcg::tri::io::ExporterPLY<TriMeshType>::Save(mesh,"test.ply");

        outline2Vec.clear();
        std::set<std::pair<size_t,size_t> > FaceEdge;

        for (size_t i=0;i<mesh.face.size();i++)
        {
            //if (mesh.face[i].IsS())continue;
            for (size_t j=0;j<3;j++)
            {
                if (!mesh.face[i].IsFaceEdgeS(j))continue;
                if (FaceEdge.count(std::pair<size_t,size_t> (i,j))>0)continue;
                outline2Vec.resize(outline2Vec.size()+1);
                outline3Vec.resize(outline3Vec.size()+1);
                vcg::face::Pos<FaceType> PosF(&mesh.face[i],j);
                vcg::face::Pos<FaceType> PosInit=PosF;
                std::cout<<"Starting new polyline"<<std::endl;
                do{
                    size_t IndexF=vcg::tri::Index(mesh,PosF.F());
                    size_t IndexE=PosF.E();
                    std::cout<<"IndexF"<<IndexF<<std::endl;
                    std::cout<<"IndexE"<<IndexE<<std::endl;
                    std::pair<size_t,size_t> Key(IndexF,IndexE);
                    assert(FaceEdge.count(Key)==0);
                    FaceEdge.insert(std::pair<size_t,size_t>(IndexF,IndexE));
                    vcg::Point2f UV;
                    UV.Import(PosF.F()->WT(PosF.E()).P()*scaleVal);
                    outline2Vec.back().push_back(UV);
                    outline3Vec.back().push_back(PosF.V()->P());
                    std::cout<<"Step"<<std::endl;
                    PosF.NextEdgeS();
                }while (PosF!=PosInit);
            }
        }
        std::cout<<"Number Polylines:"<<outline2Vec.size()<<std::endl;
    }

    //    static vcg::Box2<ScalarType> GetBox(TriMeshType &mesh,
    //                                        float scaleVal,
    //                                        float boundSize)
    //    {
    //        vcg::Box2<ScalarType> UVBBox=vcg::tri::UV_Utils<TriMeshType>::PerWedgeUVBox(mesh);
    //        UVBBox.min*=scaleVal;
    //        UVBBox.max*=scaleVal;
    //        UVBBox.Offset(boundSize*2);
    //        return UVBBox;
    //    }

    static vcg::Box2<ScalarType> GetBox(const std::vector< std::vector<vcg::Point2f> > &Outline)
    {
        vcg::Box2<ScalarType> UVBBox;
        for (size_t i=0;i<Outline.size();i++)
            for (size_t j=0;j<Outline[i].size();j++)
                UVBBox.Add(Outline[i][j]);
        return UVBBox;
    }


    static void CreateOffsetBoundary(const std::vector< std::vector<vcg::Point2f> > &Outline,
                                     std::vector< std::vector<vcg::Point2f> > &OffSetted,
                                     float boundSize)
    {
        ClipperLib::ClipperOffset co;

        for (size_t i=0;i<Outline.size();i++)
        {
            ClipperLib::Path clipOutline;
            for (size_t j=0;j<Outline[i].size();j++)
                clipOutline<<ClipperLib::IntPoint(Outline[i][j].X(),Outline[i][j].Y());
            co.AddPath(clipOutline, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
        }

        ClipperLib::Paths solution;
        co.Execute(solution, boundSize);

        for (size_t i=0;i<solution.size();i++)
        {
            OffSetted.resize(OffSetted.size()+1);
            for (size_t j=0;j<solution[i].size();j++)
            {
                ClipperLib::IntPoint clipPos=solution[i][j];
                vcg::Point2f pos(clipPos.X,clipPos.Y);
                OffSetted.back().push_back(pos);
            }
        }
    }

    static vcg::Point2f GetClosestPoint(const vcg::Point2f &Pos,const std::vector< std::vector<vcg::Point2f> > &Outline)
    {
        float minD=std::numeric_limits<float>::max();
        vcg::Point2f closest=Outline[0][0];
        for (size_t i=0;i<Outline.size();i++)
            for (size_t j=0;j<Outline[i].size();j++)
            {
                ScalarType dist_test=(Pos-Outline[i][j]).Norm();
                if (dist_test>minD)continue;
                minD=dist_test;
                closest=Outline[i][j];
            }
        return closest;
    }

    static void SawingLabels(const std::vector< std::vector<CoordType> > &outline3Vec,
                             const std::vector< std::vector<vcg::Point2f> > &Outline,
                             const std::vector< std::vector<vcg::Point2f> > TextLine,
                             std::vector<std::vector<vcg::Point2f> > &Pos2D,
                             std::vector<std::vector<std::string> > &label,
                             float space)
    {
        //        //merge because there might be duplicate on symmetry line
        //        TriMeshType merged_mesh;
        //        vcg::tri::Append<TriMeshType,TriMeshType>::Mesh(merged_mesh,mesh);
        //        vcg::tri::Clean<TriMeshType>::RemoveDuplicateVertex(merged_mesh);
        //        merged_mesh.UpdateAttributes();
        //        std::set<CoordType> BorderPos;
        //        for (size_t i=0;i<merged_mesh.vert.size();i++)
        //            if (merged_mesh.vert[i].IsB())
        //                BorderPos.insert(merged_mesh.vert[i].P());

        //vcg::tri::io::ExporterPLY<TriMeshType>::Save(merged_mesh,"test.ply");

        std::set<CoordType> VisitedPos;
        std::map<CoordType,size_t> PosLabel;
        size_t currLabel=0;

        Pos2D.resize(Outline.size());
        label.resize(Outline.size());

        for (size_t i=0;i<Outline.size();i++)
        {
            float curr_sp=0;
            size_t sizeOut=Outline[i].size();
            for (size_t j=0;j<Outline[i].size();j++)
            {
                CoordType currPos=outline3Vec[i][j];
                vcg::Point2f posUV=Outline[i][j];
                //if already added by somebody else
                if (PosLabel.count(currPos)>0)
                {
                    //retrieve the label
                    size_t foundLabel=PosLabel[currPos];

                    //add the marker
                    label[i].push_back("X");
                    Pos2D[i].push_back(posUV);

                    label[i].push_back(std::to_string((int)foundLabel));
                    posUV=GetClosestPoint(posUV,TextLine);
                    Pos2D[i].push_back(posUV);
                    //reset lenght to zero
                    curr_sp=0;
                }
                else
                {
                    //not already added the position
                    if (VisitedPos.count(currPos)==0)
                    {
                        float dist=(Outline[i][j]-Outline[i][(j+1)%sizeOut]).Norm();
                        curr_sp+=dist;
                        if (curr_sp>space)
                        {
                            label[i].push_back("X");
                            Pos2D[i].push_back(posUV);

                            label[i].push_back(std::to_string((int)currLabel));
                            posUV=GetClosestPoint(posUV,TextLine);
                            Pos2D[i].push_back(posUV);
                            PosLabel[currPos]=currLabel;
                            currLabel++;
                            curr_sp=0;
                        }
                    }
                }
                //add the position
                VisitedPos.insert(outline3Vec[i][j]);
            }
        }
    }

public:

    static void ExportUVPolyline(TriMeshType &mesh,
                                 const char *pathSVG,
                                 const char *pathPNG,
                                 QImage &SVGTxt,
                                 float scaleVal=1000,
                                 float penwidth=2,
                                 float boundSize=15,
                                 float fontsize=7)
    {
        //vcg::Box2<ScalarType> uv_box=vcg::tri::UV_Utils<CMesh>::PerWedgeUVBox(mesh);

        std::vector< std::vector<vcg::Point2f> > outline2Vec;
        std::vector< std::vector<CoordType> > outline3Vec;
        GetOutLines(mesh,outline2Vec,outline3Vec,scaleVal);


        Outline2Dumper::Param pp;
        pp.penWidth=penwidth;
        pp.fontSize=fontsize;

        std::vector<std::vector<std::string> > Label;
        std::vector<std::vector<float> > LabelRad;
        std::vector<std::vector<vcg::Similarity2f> > trText;

        trText.resize(outline2Vec.size());
        LabelRad.resize(outline2Vec.size());
        Label.resize(outline2Vec.size());

        if (fontsize!=0)
        {
            //create the line to position the text
            std::vector< std::vector<vcg::Point2f> > TextLine;
            CreateOffsetBoundary(outline2Vec,TextLine,boundSize/2);

            std::vector<std::vector<vcg::Point2f> > Pos2D;
            std::vector<std::vector<std::string> > labelInt;
            //SawingLabels(mesh,outline3Vec,outline2Vec,TextLine,Pos2D,labelInt,fontsize*4);
            SawingLabels(outline3Vec,outline2Vec,TextLine,Pos2D,labelInt,fontsize*4);
            for (size_t i=0;i<labelInt.size();i++)
                for (size_t j=0;j<labelInt[i].size();j++)
                {
                    vcg::Similarity2f Sim;
                    Sim.tra=Pos2D[i][j];

                    trText[i].push_back(Sim);
                    LabelRad[i].push_back(fontsize*2);
                    //Label[i].push_back(std::to_string((int)labelInt[i][j]));
                    Label[i].push_back(labelInt[i][j]);
                }
        }


        //create the Offset
        if (boundSize!=0)
        {
            std::vector< std::vector<vcg::Point2f> > OffSetted;
            CreateOffsetBoundary(outline2Vec,OffSetted,boundSize);
            outline2Vec.insert(outline2Vec.end(),OffSetted.begin(),OffSetted.end());
        }

        std::vector<vcg::Similarity2f> trVec;
        trVec.resize(outline2Vec.size());
        trText.resize(outline2Vec.size());
        LabelRad.resize(outline2Vec.size());
        Label.resize(outline2Vec.size());

        //get the bounding box
        vcg::Box2<ScalarType> UVBBox=GetBox(outline2Vec);
        pp.height=UVBBox.DimY();
        pp.width=UVBBox.DimX();

        std::cout<<"UV Box Dim Y:"<<UVBBox.DimY()<<std::endl;
        std::cout<<"UV Box Dim X:"<<UVBBox.DimX()<<std::endl;


        std::vector< std::vector< std::vector<Point2f> > > outline2VecVec(outline2Vec.size());
        for(size_t i=0;i<outline2Vec.size();++i)
        {
            outline2VecVec[i].resize(1);
            outline2VecVec[i][0]=outline2Vec[i];
        }

        Outline2Dumper::dumpOutline2VecSVG(pathSVG,outline2VecVec,trVec,Label,trText,LabelRad,pp);


        SVGTxt = QIcon(pathSVG).pixmap(QSize(4000,4000)).toImage();
        for (size_t x=0;x<SVGTxt.width()*SVGTxt.height();x++)
            {
                size_t baseI=x*4;
                size_t alpha=baseI+3;
                if (SVGTxt.bits()[alpha]==0)
                {
                    SVGTxt.bits()[baseI]=255;
                    SVGTxt.bits()[baseI+1]=255;
                    SVGTxt.bits()[baseI+2]=255;
                    SVGTxt.bits()[baseI+3]=255;
                }
            }

        // Save, image format based on file extension
        SVGTxt.save(pathPNG);
    }
};

#endif
