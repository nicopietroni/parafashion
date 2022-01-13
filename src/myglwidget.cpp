/****************************************************************************
**
** Copyright (C) 2011 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
**     the names of its contributors may be used to endorse or promote
**     products derived from this software without specific prior written
**     permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
** $QT_END_LICENSE$
**
****************************************************************************/

#include <GL/glew.h>
#include <QMouseEvent>

#include <math.h>
#include <tracing/GL_mesh_drawing.h>
#include <tracing/mesh_type.h>
#include "myglwidget.h"
#include <wrap/qt/trackball.h>
#include <wrap/qt/anttweakbarMapper.h>
#include <wrap/gl/gl_field.h>
#include <qtimer.h>
#include <qfont.h>
#include <QGLWidget>
#include <QGuiApplication>
#include <QOpenGLTexture>

#include "trace_path_GL.h"
#include "wrap/qt/Outline2ToQImage.h"
#include <svg_exporter.h>
#include "parafashion.h"
#include "parafashion_interface.h"

std::string pathRef="";
std::string pathDef="";
std::string pathFrames="";

vcg::Trackball track;//the active manipulator

bool drawfield=false;


typename TraceMesh::CoordType CenterDef,CenterRef;
TraceMesh deformed_mesh;
TraceMesh reference_mesh;
TraceMesh half_def_mesh;

TwBar *barFashion;

vcg::GlTrimesh<TraceMesh> glWrap;

vcg::GLField<TraceMesh> glField;

typedef typename TraceMesh::ScalarType ScalarType;
typedef typename TraceMesh::CoordType CoordType;
typedef typename TraceMesh::FaceType TraceFaceType;
typedef typename vcg::face::Pos<TraceFaceType> TracePosType;

int Iterations;
ScalarType EdgeStep;
bool drawRefMesh=true;
bool drawDefMesh=true;
bool drawSymmetryPlane=true;
bool parametrized=false;
bool colored_distortion=false;
bool draw3D=true;
bool textured=false;
bool drawParam=false;
bool matchcurvature=false;
bool hasFrames=false;

int xMouse,yMouse;

//typedef FieldSmoother<TraceMesh> FieldSmootherType;
//FieldSmootherType::SmoothParam FieldParam;
vcg::GLW::DrawMode drawmode=vcg::GLW::DMSmooth;     /// the current drawmode

bool do_rotate=false;
bool do_anim=false;
int curr_frame=0;

QTimer *timerRot;
QTimer *timerAnim;

bool HasTxt=false;

PathGL<TraceMesh> GPath;

//#define MAX_DIST 0.05

Parafashion<TraceMesh> PFashion(deformed_mesh,reference_mesh);
AnimationManager<TraceMesh> AManager(deformed_mesh);

//std::vector<std::vector<TracePosType > > TestPosSeq;

QImage SVGTxt;
GLuint layoutTxtIdx=0;
bool HasLayoutTxt=false;
bool has_to_update_layout=false;

#define COLOR_LIMITS 1.2

void GLDrawPoints(const std::vector<CoordType> &DrawPos,
                  const ScalarType &GLSize,const vcg::Color4b &Col)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glDepthRange(0,0.9995);
    glPointSize(GLSize);

    for (size_t i=0;i<DrawPos.size();i++)
    {

        vcg::glColor(Col);
        glBegin(GL_POINTS);
        vcg::glVertex(DrawPos[i]);
        glEnd();
    }
    glPopAttrib();
}


//void GlDrawPosSeq()
//{
//    glPushAttrib(GL_ALL_ATTRIB_BITS);
//    glDisable(GL_LIGHTING);
//    glDepthRange(0,0.9995);
//    glLineWidth(20);
//    glPointSize(40);

//    for (size_t i=0;i<TestPosSeq.size();i++)
//    {
//        vcg::Color4b col;
//        col=vcg::Color4b::Scatter(TestPosSeq.size(),i);
//        vcg::glColor(col);
//        glBegin(GL_LINES);
//        for (size_t j=0;j<TestPosSeq[i].size();j++)
//        {
//           CoordType P0=TestPosSeq[i][j].V()->P();
//           CoordType P1=TestPosSeq[i][j].VFlip()->P();
//           vcg::glVertex(P0);
//           vcg::glVertex(P1);
//        }
//        glEnd();
//    }

//    vcg::Color4b col=vcg::Color4b::Red;
//    vcg::glColor(col);
//    glBegin(GL_POINTS);
//    for (size_t i=0;i<deformed_mesh.vert.size();i++)
//    {
//        if (!deformed_mesh.vert[i].IsS())continue;
//        CoordType P=deformed_mesh.vert[i].P();
//        vcg::glVertex(P);
//    }
//    glEnd();

//    glPopAttrib();
//}

template <class ScalarType>
void GlDrawPlane(const vcg::Plane3<ScalarType> &Pl,
                 const ScalarType &size)
{
    typedef typename vcg::Point3<ScalarType> CoordType;
    CoordType p[4];
    p[0]=CoordType(-size,-size,0);
    p[1]=CoordType(size,-size,0);
    p[2]=CoordType(size,size,0);
    p[3]=CoordType(-size,size,0);
    CoordType N0=CoordType(0,0,1);
    CoordType N1=Pl.Direction();

    ///then get rotation matrix
    vcg::Matrix33<ScalarType> RotPl=vcg::RotationMatrix(N0,N1);
    ///then rotate-translate the points to math the right direction
    for (int i=0;i<4;i++)
    {
        p[i]=RotPl*p[i];
        CoordType Transl(0,0,0);
        Transl+=Pl.Direction()*Pl.Offset();
        p[i]+=Transl;
        //p[i]+=Center;
    }

    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glDisable(GL_LIGHTING);
    vcg::glColor(vcg::Color4b(255,0,0,255));
    glLineWidth(3);
    glBegin(GL_LINE_LOOP);
    for (int i=0;i<4;i++)
        vcg::glVertex(p[i]);
    glEnd();

    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glEnable(GL_BLEND);                //activate blending mode
    //glBlendFunc(GL_SRC_ALPHA,GL_ONE);  //define blending factors
    glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);
    vcg::glColor(vcg::Color4b(255,0,0,100));

    glBegin(GL_QUADS);
    for (int i=0;i<4;i++)
        vcg::glVertex(p[i]);
    glEnd();

    glPopAttrib();
}

void GLDrawPatchEdges(vcg::Color4b col=vcg::Color4b(0,0,255,255),
                      ScalarType GLSize=5)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glDepthRange(0,0.9999);
    glLineWidth(GLSize);
    vcg::glColor(col);
    glBegin(GL_LINES);
    for (size_t i=0;i<deformed_mesh.face.size();i++)
    {
        //size_t IndexP0=deformed_mesh.face[i].Q();
        for (size_t j=0;j<3;j++)
        {
            bool drawEdge=false;
            drawEdge|=vcg::face::IsBorder(deformed_mesh.face[i],j);
            drawEdge|=deformed_mesh.face[i].IsFaceEdgeS(j);
            //            size_t IndexP1=deformed_mesh.face[i].FFp(j)->Q();
            //            drawEdge|=(IndexP1!=IndexP0);
            if (!drawEdge)continue;
            CoordType Pos0=deformed_mesh.face[i].P0(j);
            CoordType Pos1=deformed_mesh.face[i].P1(j);
            vcg::glVertex(Pos0);
            vcg::glVertex(Pos1);
        }
    }
    glEnd();
    glPopAttrib();
}

void MyGLWidget::GLDrawLegenda()
{
    vcg::Color4b col0=vcg::Color4b::ColorRamp(-COLOR_LIMITS,COLOR_LIMITS,-1);
    vcg::Color4b colM=vcg::Color4b::ColorRamp(-COLOR_LIMITS,COLOR_LIMITS,0);
    vcg::Color4b col1=vcg::Color4b::ColorRamp(-COLOR_LIMITS,COLOR_LIMITS,1);

    ScalarType sizeX=0.5;
    ScalarType sizeY=0.05;
    glPushMatrix();
    vcg::glTranslate(CoordType(0,-1,0));
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);

    //renderText(floor(-2*sizeX), floor(-2*sizeY), "<5% Compression", QFont("Arial", 16, QFont::Bold,false) );

    glBegin(GL_QUADS);
    //glColor(vcg::Color4b(0,0,255,255));
    glColor(col1);
    glVertex(CoordType(-sizeX,-sizeY,0));
    //glColor(vcg::Color4b(0,255,0,255));
    glColor(colM);
    glVertex(CoordType(0,-sizeY,0));
    glVertex(CoordType(0,sizeY,0));
    //glColor(vcg::Color4b(0,0,255,255));
    glColor(col1);
    glVertex(CoordType(-sizeX,sizeY,0));

    //glColor(vcg::Color4b(0,255,0,255));
    glColor(colM);
    glVertex(CoordType(0,-sizeY,0));
    //glColor(vcg::Color4b(255,0,0,255));
    glColor(col0);
    glVertex(CoordType(sizeX,-sizeY,0));
    glVertex(CoordType(sizeX,sizeY,0));
    //glColor(vcg::Color4b(0,255,0,255));
    glColor(colM);
    glVertex(CoordType(0,sizeY,0));
    glEnd();

    glPopAttrib();
    glPopMatrix();

    vcg::glColor(vcg::Color4b(0,0,0,255));

    std::string ComprString=std::to_string(((int)fabs(PFashion.max_compression*100)));
    ComprString+="% Compression";
    int Width=width();
    int Height=height();
    double posX=(double)Width*0.25;
    double posY=(double)Height*0.9;
    renderText(floor(posX), floor(posY), ComprString.c_str(), QFont("Arial", 16, QFont::Bold,false) );

    std::string TensionString=std::to_string(((int)fabs(PFashion.max_compression*100)));
    TensionString+="% Tension";
    posX=(double)Width*0.65;
    posY=(double)Height*0.9;
    renderText(floor(posX), floor(posY), TensionString.c_str(), QFont("Arial", 16, QFont::Bold,false) );
}

void DoColorByDistortion()
{
    if (!parametrized)return;
    //save old partitioning first
    std::vector<ScalarType> FaceQ;
    for (size_t i=0;i<deformed_mesh.face.size();i++)
        FaceQ.push_back(deformed_mesh.face[i].Q());

    vcg::tri::Distortion<TraceMesh,true>::SetQasDistorsion(deformed_mesh,vcg::tri::Distortion<TraceMesh,true>::EdgeComprStretch);

    //copy per vertex
    vcg::tri::UpdateQuality<TraceMesh>::VertexFromFace(deformed_mesh);

    //    //clamp
    //    for (size_t i=0;i<deformed_mesh.vert.size();i++)
    //    {
    //        deformed_mesh.vert[i].Q()=std::min(deformed_mesh.vert[i].Q(),PFashion.max_tension);
    //        deformed_mesh.vert[i].Q()=std::max(deformed_mesh.vert[i].Q(),PFashion.max_);
    //    }

    colored_distortion=true;

    for (size_t i=0;i<deformed_mesh.vert.size();i++)
    {
        ScalarType Val=deformed_mesh.vert[i].Q();
        if (deformed_mesh.vert[i].Q()<=PFashion.max_compression)
        {
            deformed_mesh.vert[i].C()=vcg::Color4b::Blue;
            continue;
        }
        if (deformed_mesh.vert[i].Q()>=PFashion.max_tension)
        {
            deformed_mesh.vert[i].C()=vcg::Color4b::Red;
            continue;
        }
        deformed_mesh.vert[i].C()=vcg::Color4b::ColorRamp(PFashion.max_compression*COLOR_LIMITS,
                                                         PFashion.max_tension*COLOR_LIMITS,
                                                         Val);
    }

    //vcg::tri::UpdateColor<TraceMesh>::PerVertexQualityRamp(deformed_mesh,MAX_DIST,-MAX_DIST);

    //cut everything bigger than a threshold
    for (size_t i=0;i<deformed_mesh.face.size();i++)
    {
        deformed_mesh.face[i].Q()=FaceQ[i];
    }
}

void DoColorByPatch()
{
    //SelectPatchBorderEdges(deformed_mesh);
    MakePartitionOnQConsistent(deformed_mesh);
    deformed_mesh.ScatterColorByQualityFace();
}

void TW_CALL ColorByNone(void *)
{
    vcg::tri::UpdateColor<TraceMesh>::PerFaceConstant(deformed_mesh,vcg::Color4b(255,255,255,255));
}

void TW_CALL ColorByArapDist(void *)
{
    DoColorByDistortion();
}

void TW_CALL ColorByPatch(void *)
{
    DoColorByPatch();
    colored_distortion=false;
}

void DoSymmetrize()
{
    PFashion.MakeMeshSymmetric(GPath.PickedPoints);
}

void TW_CALL SymmetrizeDeformed(void *)
{
    DoSymmetrize();
}

void DoParametrize()
{

    PFashion.DoParametrize();
}

void TW_CALL Parametrize(void *)
{
    DoParametrize();
    drawParam=true;
}

void DoTracePath()
{
    PFashion.TracePatch();
    if (hasFrames)
        AManager.UpdateProjectionBasis();
    //TestPosSeq.clear();
}

void TW_CALL TracePath(void *)
{
    DoTracePath();
    //TestPosSeq.clear();
}

void DoGenerateSVG(std::string ProjM)
{
    //THEN SAVE THE PATCH DATA
    std::string pathPartitions=ProjM;
    pathPartitions=ProjM+"_patch.txt";
    FILE *F=fopen(pathPartitions.c_str(),"wt");
    assert(F!=NULL);
    fprintf(F,"%d\n",(int)deformed_mesh.face.size());
    for (size_t i=0;i<deformed_mesh.face.size();i++)
        fprintf(F,"%d\n",(int)deformed_mesh.face[i].Q());
    fclose(F);

    //SAVE THE SVG
    std::string pathPatch=ProjM;
    pathPatch=ProjM+"_patch.svg";
    std::string pathPatchPNG=ProjM;
    pathPatchPNG=ProjM+"_patch.png";
    float scale=1000;
    SvgExporter<TraceMesh>::ExportUVPolyline(deformed_mesh,
                                             pathPatch.c_str(),
                                             pathPatchPNG.c_str(),
                                             SVGTxt,
                                             scale,4,
                                             PFashion.param_boundary*scale,
                                             PFashion.param_boundary*scale*0.75);
    //prepare the texure
    HasLayoutTxt=true;
    has_to_update_layout=true;
}

void TW_CALL GenerateSVG(void *)
{
    std::string ProjM=pathDef;
    size_t indexExt=ProjM.find_last_of(".");
    ProjM=ProjM.substr(0,indexExt);
    DoGenerateSVG(ProjM);
}

void TW_CALL SaveData(void *)
{

    PFashion.SaveDebugPatches("../data/");

    //get the path
    std::string ProjM=pathDef;
    size_t indexExt=ProjM.find_last_of(".");
    ProjM=ProjM.substr(0,indexExt);
    std::string saveMeshName=ProjM+std::string("_patch.obj");

    //SAVE THE MESH
    TraceMesh saveM;
    vcg::tri::Append<TraceMesh,TraceMesh>::Mesh(saveM,deformed_mesh);


    vcg::tri::Clean<TraceMesh>::RemoveDuplicateVertex(saveM);

    for (size_t i=0;i<saveM.vert.size();i++)
        saveM.vert[i].P()+=CenterDef;

    saveM.UpdateAttributes();
    vcg::tri::io::ExporterOBJ<TraceMesh>::Save(saveM,saveMeshName.c_str(),
                                               vcg::tri::io::Mask::IOM_WEDGTEXCOORD|
                                               vcg::tri::io::Mask::IOM_FACECOLOR);

    //THEN SAVE THE PATCH DATA
    std::string pathPartitions=ProjM;
    pathPartitions=ProjM+"_patch.txt";
    FILE *F=fopen(pathPartitions.c_str(),"wt");
    assert(F!=NULL);
    fprintf(F,"%d\n",(int)deformed_mesh.face.size());
    for (size_t i=0;i<deformed_mesh.face.size();i++)
        fprintf(F,"%d\n",(int)deformed_mesh.face[i].Q());
    fclose(F);

    //    //SAVE THE SVG
    //    std::string pathPatch=ProjM;
    //    pathPatch=ProjM+"_patch.svg";
    //    std::string pathPatchPNG=ProjM;
    //    pathPatchPNG=ProjM+"_patch.png";
    //    float scale=1000;
    //    SvgExporter<TraceMesh>::ExportUVPolyline(deformed_mesh,
    //                                             pathPatch.c_str(),
    //                                             pathPatchPNG.c_str(),
    //                                             SVGTxt,
    //                                             scale,4,
    //                                             PFashion.param_boundary*scale,
    //                                             PFashion.param_boundary*scale*0.75);
    //    //prepare the texure
    //    HasLayoutTxt=true;
    //    has_to_update_layout=true;

    DoGenerateSVG(ProjM);

    //THEN save the mesh in UV
    std::string pathUV=ProjM;
    pathUV=ProjM+"_UV.txt";
    F=fopen(pathUV.c_str(),"wt");
    assert(F!=NULL);
    fprintf(F,"%d\n",(int)deformed_mesh.face.size());
    for (size_t i=0;i<deformed_mesh.face.size();i++)
    {
        vcg::Point2<ScalarType> T0=deformed_mesh.face[i].WT(0).P();
        vcg::Point2<ScalarType> T1=deformed_mesh.face[i].WT(1).P();
        vcg::Point2<ScalarType> T2=deformed_mesh.face[i].WT(2).P();
        fprintf(F,"%f,%f;%f,%f;%f,%f\n",
                T0.X(),T0.Y(),
                T1.X(),T1.Y(),
                T2.X(),T2.Y());
    }
    fclose(F);

}

void DoSmoothField()
{
    PFashion.ComputeField();

}

void TW_CALL RemoveAlongSymmetryLine(void *)
{
    //    for (size_t i=0;i<deformed_mesh.face.size();i++)
    //        for (size_t j=0;j<3;j++)
    //        {
    //            if (!vcg::face::IsBorder(deformed_mesh.face[i],j))continue;
    //            deformed_mesh.face[i].SetFaceEdgeS(j);
    //        }
    //    vcg::tri::Clean<TraceMesh>::RemoveDuplicateVertex(deformed_mesh);
    //    deformed_mesh.UpdateAttributes();
    //    RetrievePosSeqFromSelEdges(deformed_mesh,TestPosSeq);
    //vcg::tri::io::ExporterPLY<TraceMesh>::Save(deformed_mesh,"test_mesh.ply");

    PFashion.RemoveOnSymmetryPathIfPossible();
    DoParametrize();
    DoColorByPatch();

}


void TW_CALL SmoothField(void *)
{
    DoSmoothField();
    drawfield=true;
}

void UpdateBaseColorMesh()
{
    if (colored_distortion)
        DoColorByDistortion();
    else
        DoColorByPatch();
}

void DoBatchProcess()
{

    PFashion.BatchProcess(GPath.PickedPoints);
    if (hasFrames)
    {
        AManager.UpdateProjectionBasis();
        //AManager.UpdateCurvatureAndStretch();
    }

    drawDefMesh=true;
    drawRefMesh=false;
    drawParam=true;
    parametrized=true;

    UpdateBaseColorMesh();

    //TestPosSeq.clear();
}

void TW_CALL BatchProcess(void *)
{
    DoBatchProcess();
}


void SetFieldBarSizePosition(QWidget *w)
{
    int params[2];
    params[0] = QTDeviceWidth(w) / 2.0; // AntTweakBar Menu size here
    params[1] = QTDeviceHeight(w)*1.5;
    TwSetParam(barFashion, NULL, "size", TW_PARAM_INT32, 2, params);
    params[0] = QTLogicalToDevice(w, 10);
    params[1] = 30;//QTDeviceHeight(w) - params[1] - QTLogicalToDevice(w, 10);
    TwSetParam(barFashion, NULL, "position", TW_PARAM_INT32, 2, params);
}

enum FieldAnimMode{FANone,FACurvature,FAStretchCompress};
FieldAnimMode FAnimMode=FANone;

void InitBar(QWidget *w) // AntTweakBar menu
{
    barFashion = TwNewBar("parafashion menu");

    SetFieldBarSizePosition(w);

    TwEnumVal drawmodes[4] = { {vcg::GLW::DMSmooth, "Smooth"},
                               {vcg::GLW::DMPoints, "Per Points"},
                               {vcg::GLW::DMFlatWire, "FlatWire"},
                               {vcg::GLW::DMFlat, "Flat"}};
    // Create a type for the enum shapeEV
    TwType drawMode = TwDefineEnum("DrawMode", drawmodes, 4);
    TwAddVarRW(barFashion, "Draw Mode", drawMode, &drawmode, " keyIncr='<' keyDecr='>' help='Change draw mode.' ");

    TwAddVarRW(barFashion,"draw3D",TW_TYPE_BOOLCPP, &draw3D," label='Draw 3D Mesh'");
    TwAddVarRW(barFashion,"drawUV",TW_TYPE_BOOLCPP, &drawParam," label='Draw UV Mesh'");
    TwAddVarRW(barFashion,"textured",TW_TYPE_BOOLCPP, &textured," label='Textured'");

    TwAddVarRW(barFashion,"doRotate",TW_TYPE_BOOLCPP, &do_rotate," label='Rotate'");

    TwEnumVal field_anim_mode[3] = { {FANone, "Field Anim None"},
                                     {FACurvature, "Field Anim Curv"},
                                     {FAStretchCompress, "Field Anim Stretch/Compress"}
                                   };
    TwType fieldAnimMode = TwDefineEnum("FieldAnimMode", field_anim_mode, 3);
    TwAddVarRW(barFashion, "Field Anim Mode", fieldAnimMode, &FAnimMode, " keyIncr='<' keyDecr='>' help='Change Field Anim mode.' ");

    TwAddVarRW(barFashion,"doAnim",TW_TYPE_BOOLCPP, &do_anim," label='Animate'");



    TwAddButton(barFashion,"ColorDist",ColorByArapDist,0,"label='Color Distortion'");
    TwAddButton(barFashion,"ColorPatch",ColorByPatch,0,"label='Color Patch'");
    TwAddButton(barFashion,"ColorNone",ColorByNone,0,"label='Color None'");


    TwAddVarRW(barFashion,"drawDef",TW_TYPE_BOOLCPP, &drawDefMesh," label='Draw Deformed'");
    TwAddVarRW(barFashion,"drawRef",TW_TYPE_BOOLCPP, &drawRefMesh," label='Draw Reference'");
    TwAddVarRW(barFashion,"drawField",TW_TYPE_BOOLCPP, &drawfield," label='Draw Field'");
    TwAddVarRW(barFashion,"drawPlane",TW_TYPE_BOOLCPP, &drawSymmetryPlane," label='Draw Symm Plane'");

    TwAddSeparator(barFashion,NULL,NULL);

    TwAddButton(barFashion,"SymmetrizeDef",SymmetrizeDeformed,0,"label='Symmetrize Deformed'");
    //    TwEnumVal fieldmodes[2] = { {FMBoundary, "Boundary"},
    //                                {FMCurvature, "Curvature"}
    //                              };
    //TwType fieldMode = TwDefineEnum("FieldMode", fieldmodes, 2);
    //TwAddVarRW(barFashion, "Field Mode", fieldMode, &FMode, " keyIncr='<' keyDecr='>' help='Change field mode.' ");

    TwAddButton(barFashion,"ComputeField",SmoothField,0,"label='Compute Field'");
    TwAddVarRW(barFashion,"matchCurv",TW_TYPE_BOOLCPP,
               &PFashion.match_valence," label='Match Valence'");

    TwEnumVal patchmode[3] = { {PMMinTJuncions, "Min T-Junctions"},
                               {PMAvgTJuncions, "Avg T-Junctions"},
                               {PMAllTJuncions, "All T-Junctions"}
                             };
    TwType patchMode = TwDefineEnum("PatchMode", patchmode, 3);
    TwAddVarRW(barFashion, "Patch Mode", patchMode, &PFashion.PMode, " keyIncr='<' keyDecr='>' help='Change patch mode.' ");


    TwAddVarRW(barFashion,"MaxCorners",TW_TYPE_INT32,&PFashion.max_corners," label='Max Corners'");
    TwAddVarRW(barFashion,"SelfGlue",TW_TYPE_BOOLCPP,&PFashion.allow_self_glue," label='Allow SelfGlue'");
    TwAddVarRW(barFashion,"Darts",TW_TYPE_BOOLCPP,&PFashion.use_darts," label='Allow Darts'");


    TwAddButton(barFashion,"TracePaths",TracePath,0,"label='Trace Paths'");
    TwAddVarRW(barFashion,"ParamBound",TW_TYPE_DOUBLE,
               &PFashion.param_boundary," label='Boundary Tolerance'");

    TwAddButton(barFashion,"Parametrize",Parametrize,0,"label='Parametrize Deformed'");

    TwAddSeparator(barFashion,NULL,NULL);
    TwAddVarRW(barFashion,"RemoveSym",TW_TYPE_BOOLCPP,&PFashion.remove_along_symmetry," label='Rem Symmetry'");
    TwAddVarRW(barFashion,"checkStress",TW_TYPE_BOOLCPP,&PFashion.check_stress," label='Check Stress'");
    TwAddVarRW(barFashion,"MaxCompr",TW_TYPE_DOUBLE,&PFashion.max_compression," label='Max Compression'");
    TwAddVarRW(barFashion,"MaxTens",TW_TYPE_DOUBLE,&PFashion.max_tension," label='Max Tension'");


    TwAddButton(barFashion,"BatchProcess",BatchProcess,0,"label='Batch Process'");
    //TwAddButton(barFashion,"BatchProcess 2",BatchProcess2,0,"label='Batch Process 2'");
    TwAddButton(barFashion,"RemoveAlongSym",RemoveAlongSymmetryLine,0,"label='Remove Along Symmetry'");
    TwAddButton(barFashion,"GenerateSVG",GenerateSVG,0,"label='Generate SVG'");


    TwAddButton(barFashion,"SaveData",SaveData,0,"label='Save Data'");


}

ScalarType angleR=0;
ScalarType alpha=0.4;

void MyGLWidget::UpdateAnim()
{
    if (!hasFrames)return;

    curr_frame=(curr_frame+1)%AManager.NumFrames();

    if (FAnimMode==FANone)
    {
        AManager.UpdateToFrame(curr_frame,false,false);
        UpdateBaseColorMesh();
    }
    if (FAnimMode==FACurvature)
    {
        AManager.UpdateToFrame(curr_frame,true,false);
        AManager.ColorByAnisotropy();
    }
    if (FAnimMode==FAStretchCompress)
    {
        AManager.UpdateToFrame(curr_frame,false,true);
        AManager.ColorByStretch();
    }

    update();
}

void MyGLWidget::UpdateRot()
{
    angleR+=alpha;
    update();
}

QImage txt_image;
//GLuint texture_dress;

//QOpenGLTexture * texture = nullptr;

MyGLWidget::MyGLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
    timerRot= new QTimer(this);
    connect(timerRot, SIGNAL(timeout()), this, SLOT(UpdateRot()));

    timerAnim=new QTimer(this);
    connect(timerAnim, SIGNAL(timeout()), this, SLOT(UpdateAnim()));


    bool Loaded=deformed_mesh.LoadMesh(pathDef.c_str());
    //bool Loaded=LoadTriMesh(deformed_mesh,pathDef.c_str());
    if (!Loaded)
    {
        std::cout<<"Error Loading Mesh"<<std::endl;
        exit(0);
    }
    std::cout<<"Loaded "<<deformed_mesh.face.size()<<" faces of deformed mesh "<<std::endl;
    std::cout<<"Loaded "<<deformed_mesh.vert.size()<<" vertices of deformed mesh  "<<std::endl;
    deformed_mesh.UpdateAttributes();

    Loaded=reference_mesh.LoadMesh(pathRef.c_str());
    //Loaded=LoadTriMesh(reference_mesh,pathRef.c_str());
    if (!Loaded)
    {
        std::cout<<"Error Loading Mesh"<<std::endl;
        exit(0);
    }
    std::cout<<"Loaded "<<reference_mesh.face.size()<<" faces of reference mesh "<<std::endl;
    std::cout<<"Loaded "<<reference_mesh.vert.size()<<" vertices of reference mesh  "<<std::endl;
    reference_mesh.UpdateAttributes();

    CenterRef=reference_mesh.MoveCenterOnZero();
    CenterDef=deformed_mesh.MoveCenterOnZero();

    vcg::tri::UpdateColor<TraceMesh>::PerFaceConstant(deformed_mesh,vcg::Color4b(220,220,220,255));
    vcg::tri::UpdateColor<TraceMesh>::PerFaceConstant(reference_mesh,vcg::Color4b(255,185,15,255));

    //see if there is the texture
    std::string pathTxt=pathDef;
    pathTxt.erase(pathTxt.find_last_of("/"));
    pathTxt.append("/txt.png");
    //    texture = LoadDDS(textureFile);
    //	Q_ASSERT(texture != nullptr);
    txt_image=QImage(pathTxt.c_str());
    HasTxt=(!txt_image.isNull());
    if (HasTxt)
        std::cout<<"Texture Loaded "<<std::endl;

    hasFrames=false;
    hasDoubleClick=false;

    if (pathFrames!=std::string(""))
    {
        AManager.Init();
        std::cout<<"Loading Frames"<<std::endl;
        hasFrames=AManager.LoadPosFrames(pathFrames.c_str());
        if (hasFrames)
            std::cout<<"Frames Loaded Successfully"<<std::endl;
        else
        {
            std::cout<<"Error Loading Frames"<<std::endl;
            exit(0);
        }
        std::cout<<"Update X Frames curvature and stretch"<<std::endl;
        //AManager.UpdateCurvatureAndStretch();
    }

    PFashion.Init();
}


GLuint texture_dress;

void MyGLWidget::initializeGL ()
{
    //initialize Glew
    glewInit();
    //CaptInt.GLInit( MyGLWidget::width(),MyGLWidget::height());
    glClearColor(0, 0, 0, 0);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

    if (HasTxt)
    {
        //        glGenTextures(1, &texture_dress);
        //        glBindTexture(GL_TEXTURE_2D, texture_dress);

        //        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, txt_image.width(), txt_image.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, txt_image.constBits());
        //            // glGenerateMipmap(GL_TEXTURE_2D);

        //        glActiveTexture(GL_TEXTURE0);
        //        //glBindTexture(GL_TEXTURE_2D, texture_dress);

        glGenTextures( 1, & texture_dress );
        glEnable(GL_TEXTURE_2D);
        glBindTexture( GL_TEXTURE_2D, texture_dress );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
        glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        //load texture at level i
        //txt_image.
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, txt_image.width(), txt_image.height(), 0, GL_BGRA_EXT, GL_UNSIGNED_BYTE, txt_image.constBits());
        glDisable(GL_TEXTURE_2D);
        glBindTexture( GL_TEXTURE_2D, 0 );
    }
}


void MyGLWidget::resizeGL (int w, int h)
{
    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    TwWindowSize(w, h);
    InitBar(this);
    initializeGL();
}


static void GLDrawSVGLayout(size_t sizeU,
                            size_t sizeV,
                            GLuint TxtIndex)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glPushMatrix();

    ScalarType MaxSize=std::max(sizeU,sizeV);
    ScalarType DimX=sizeU/MaxSize;
    ScalarType DimY=sizeV/MaxSize;

    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);

    if (TxtIndex>=0)
    {
        glActiveTexture(GL_TEXTURE0);
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, TxtIndex);
        //std::cout<<"USE TXT"<<std::endl;
    }
    vcg::glColor(vcg::Color4b(255,255,255,255));
    glBegin(GL_QUADS);
    vcg::glTexCoord(vcg::Point2<ScalarType>(0,0));
    vcg::glVertex(CoordType(0,0,0));
    vcg::glTexCoord(vcg::Point2<ScalarType>(1,0));
    vcg::glVertex(CoordType(DimX,0,0));
    vcg::glTexCoord(vcg::Point2<ScalarType>(1,1));
    vcg::glVertex(CoordType(DimX,DimY,0));
    vcg::glTexCoord(vcg::Point2<ScalarType>(0,1));
    vcg::glVertex(CoordType(0,DimY,0));
    glEnd();

    glPopMatrix();
    glPopAttrib();
}

void MyGLWidget::paintGL ()
{
    if (has_to_update_layout)
    {
        glGenTextures( 1, & layoutTxtIdx );
        glEnable(GL_TEXTURE_2D);
        glBindTexture( GL_TEXTURE_2D, layoutTxtIdx );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
        glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        //load texture at level i
        //txt_image.

        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, SVGTxt.width(), SVGTxt.height(), 0, GL_BGRA_EXT, GL_UNSIGNED_BYTE, SVGTxt.constBits());
        glDisable(GL_TEXTURE_2D);
        glBindTexture( GL_TEXTURE_2D, 0 );
        has_to_update_layout=false;
        //std::cout<<"LOADED TXT"<<std::endl;
    }

    if (do_rotate)
        timerRot->start(0);
    else
        timerRot->stop();

    if (do_anim)
        timerAnim->start(0);
    else
        timerAnim->stop();

    glClearColor(255,255,255,255);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, MyGLWidget::width()/(float)MyGLWidget::height(), 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,3.5f,   0,0,0,   0,1,0);
    track.center=vcg::Point3f(0, 0, 0);
    track.radius= 1;
    track.GetView();

    glPushMatrix();

    bool draw_uv_and_3D=(draw3D && drawParam && parametrized);

    if (colored_distortion)
        GLDrawLegenda();

    if ((parametrized)&&(drawParam))
    {
        if (draw_uv_and_3D)
        {
            glPushMatrix();
            glTranslate(CoordType(1.4, 0.6,0)); // TODO find formula to fit exactly in corner
            //vcg::Box2<ScalarType> uv_box=vcg::tri::UV_Utils<TraceMesh>::PerWedgeUVBox(deformed_mesh);
            //glTranslate(CoordType(uv_box.DimX()/3,0,0));
            vcg::glScale(0.5);
        }
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDepthRange(0.0001,1);
        MeshDrawing<TraceMesh>::GLDrawUV(deformed_mesh,textured,colored_distortion);
        //deformed_mesh.GLDrawUV(textured,colored_distortion);
        glEnable( GL_LINE_SMOOTH );
        glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
        glDepthRange(0,0.999);
        MeshDrawing<TraceMesh>::GLDrawEdgeUV(deformed_mesh);
        //deformed_mesh.GLDrawEdgeUV();

        if (draw_uv_and_3D)
        {
            glPopMatrix();
            //            glMatrixMode(GL_MODELVIEW);
            //            glLoadIdentity();
            if (HasLayoutTxt)
            {
                glPushMatrix();
                glTranslate(CoordType(0.8,-1,0));
                //vcg::glScale(0.5);
                GLDrawSVGLayout(SVGTxt.width(),SVGTxt.height(),layoutTxtIdx);
                glPopMatrix();
            }
        }
        //        if (draw_uv_and_3D)
        //        {
        //            glPopMatrix();
        //        }
        glPopAttrib();
    }


    if (draw_uv_and_3D)
    {
        glPushMatrix();
        glTranslate(CoordType(0,0,0)); // Offset for main 3D mesh when UV is also displayed here
        //vcg::glScale(2.0);
    }

    track.Apply();
    glPushMatrix();

    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    glEnable( GL_POLYGON_SMOOTH );
    glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST );

    glRotated(angleR,0,1,0);
    vcg::glScale(2.0f/reference_mesh.bbox.Diag());
    glTranslate(-reference_mesh.bbox.Center());
    //glTranslate(CoordType(-reference_mesh.bbox.Diag()/2,0,0));


    if (do_anim && hasFrames)
    {
        if (FAnimMode==FACurvature)
            vcg::GLField<TraceMesh>::GLDrawFaceField(deformed_mesh,false,false,0.007,
                                                     AManager.MaxAnisotropy(),0,false);

        if (FAnimMode==FAStretchCompress)
            vcg::GLField<TraceMesh>::GLDrawFaceField(deformed_mesh,false,false,0.007,
                                                     AManager.MaxStretchCompress(),
                                                     -AManager.MaxStretchCompress(),true);
    }

    if ((drawDefMesh)&&(draw3D))
    {

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        //glPushMatrix();

        if ((HasTxt)&&(textured))
        {
            glActiveTexture(GL_TEXTURE0);
            glEnable(GL_TEXTURE_2D);
            glBindTexture(GL_TEXTURE_2D, texture_dress);
        }
        else
            glDisable(GL_TEXTURE_2D);

        vcg::GLW::ColorMode CM=vcg::GLW::CMPerFace;
        if (colored_distortion)
            CM=vcg::GLW::CMPerVert;

        glWrap.m=&deformed_mesh;
        if ((textured)&&(HasTxt))
            glWrap.Draw(drawmode,CM,vcg::GLW::TMPerWedge);
        else
            glWrap.Draw(drawmode,CM,vcg::GLW::TMNone);

        glDisable(GL_TEXTURE_2D);

        //half_def_mesh.GLDrawSharpEdges(vcg::Color4b(255,0,0,255),10);
        MeshDrawing<TraceMesh>::GLDrawSharpEdges(deformed_mesh,vcg::Color4b(255,0,0,255),10);
        //deformed_mesh.GLDrawSharpEdges(vcg::Color4b(255,0,0,255),10);

        GLDrawPatchEdges();

        glPopAttrib();
    }

    if ((drawRefMesh)&&(draw3D))
    {
        glWrap.m=&reference_mesh;
        glWrap.Draw(drawmode,vcg::GLW::CMPerFace,vcg::GLW::TMNone);
    }

    if ((drawfield)&&(draw3D))
    {
        vcg::GLField<TraceMesh>::GLDrawFaceField(deformed_mesh,false,false,0.007);
        vcg::GLField<TraceMesh>::GLDrawSingularity(deformed_mesh);
    }

    if ((drawSymmetryPlane)&&(draw3D))
    {
        GlDrawPlane(Symmetrizer<TraceMesh>::SymmPlane(),deformed_mesh.bbox.Diag()/2);//,deformed_mesh.bbox.Center());
    }

    if  (user_is_picking)
    {
        GPath.GLAddPoint(vcg::Point2i(PickX,PickY));
        GPath.GlDrawLastPath();
        //GPath.GlDrawLastPath();
        //TPath.GLAddPoint(vcg::Point2i(PickX,PickY));
        //TPath.GlDrawLastPath();
    }
    if (hasDoubleClick)
    {
        std::cout<<"test"<<std::endl;
        bool has_removed=GPath.GLRemovePathFromPoint(vcg::Point2i(PickX,PickY),
                                                     half_def_mesh.bbox.Diag()/10);

        if (has_removed)
        {
            std::cout<<"removed"<<std::endl;
            DoBatchProcess();
        }

        hasDoubleClick=false;
    }
    //TPath.GLDrawSnapped();
    //TPath.GlDrawPath();

    //GlDrawPosSeq();

    if (draw_uv_and_3D)
    {
        glPopMatrix();
    }


    glPopMatrix();
    glPopMatrix();
    //    glPopMatrix();

    //    if (draw_uv_and_3D)
    //    {
    //        glPopMatrix();
    //    }

    TwDraw();

    //GLenum GlErr=glGetError();

    switch(glGetError()) {

    case GL_INVALID_ENUM :
        std::cout<<"Invalid Enum"<<std::endl;
        assert(0);
        break;
    case GL_INVALID_VALUE :
        std::cout<<"Invalid Value"<<std::endl;
        assert(0);
        break;
    case GL_INVALID_OPERATION :
        std::cout<<"Invalid Operation"<<std::endl;
        assert(0);
        break;
    case GL_INVALID_FRAMEBUFFER_OPERATION :
        std::cout<<"Invalid FB Operation"<<std::endl;
        assert(0);
        break;
    case GL_OUT_OF_MEMORY :
        std::cout<<"Out of Memory"<<std::endl;
        assert(0);
        break;
    case GL_STACK_UNDERFLOW :
        std::cout<<"Stuck Underflow"<<std::endl;
        assert(0);
        break;
    case GL_STACK_OVERFLOW :
        std::cout<<"Stuck Overflow"<<std::endl;
        assert(0);
        break;
        //default :
    }

    //assert(glGetError()==GL_NO_ERROR);
}

void MyGLWidget::keyReleaseEvent (QKeyEvent * e)
{
    e->ignore ();
    if (e->key () == Qt::Key_Control)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ControlModifier));
    if (e->key () == Qt::Key_Shift)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
    if (e->key () == Qt::Key_Alt) track.ButtonUp (QT2VCG (Qt::NoButton, Qt::AltModifier));
    if (e->key () == Qt::Key_Space)
    {
        spacebar_being_pressed = false;
        return; // This gets called often when spacebar is pressed... don't updateGL
    }
    updateGL ();
}


void MyGLWidget::keyPressEvent (QKeyEvent * e)
{
    e->ignore ();
    if (e->key () == Qt::Key_Control) track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ControlModifier));
    if (e->key () == Qt::Key_Shift)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
    if (e->key () == Qt::Key_Alt)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::AltModifier));
    if (e->key () == Qt::Key_Space) {
        spacebar_being_pressed = true;
    }

    TwKeyPressQt(e);
    updateGL ();
}

void MyGLWidget::mousePressEvent (QMouseEvent * e)
{
    if(!TwMousePressQt(this,e))
    {
        e->accept ();
        setFocus ();
        track.MouseDown(QT2VCG_X(this, e), QT2VCG_Y(this, e), QT2VCG (e->button (), e->modifiers ()));

        //if(e->button() == Qt::RightButton)
        //        if (QGuiApplication::keyboardModifiers().testFlag(Qt::ShiftModifier))
        //        {
        if (spacebar_being_pressed)
        {
            user_is_picking = true;
            xMouse=QT2VCG_X(this, e);
            yMouse=QT2VCG_Y(this, e);
            PickX=xMouse;
            PickY=yMouse;
            GPath.AddNewPath();
            //pointToPick=Point2i(xMouse,yMouse);
        }
    }
    updateGL ();
}

void MyGLWidget::mouseMoveEvent (QMouseEvent * e)
{
    if (e->buttons ()) {
        if (user_is_picking)
        {
            xMouse=QT2VCG_X(this, e);
            yMouse=QT2VCG_Y(this, e);
            PickX=xMouse;
            PickY=yMouse;
        }
        else
            track.MouseMove(QT2VCG_X(this, e), QT2VCG_Y(this, e));

        updateGL ();
    }
    TwMouseMotion(QTLogicalToDevice(this, e->x()), QTLogicalToDevice(this, e->y()));
}

void MyGLWidget::mouseDoubleClickEvent (QMouseEvent * e)
{
    if (e->buttons ())
    {
        xMouse=QT2VCG_X(this, e);
        yMouse=QT2VCG_Y(this, e);
        if(e->button() == Qt::LeftButton)
        {
            hasDoubleClick=true;
            PickX=xMouse;
            PickY=yMouse;
        }
        updateGL ();
    }
    updateGL();
}


void MyGLWidget::mouseReleaseEvent (QMouseEvent * e)
{
    track.MouseUp(QT2VCG_X(this, e), QT2VCG_Y(this, e), QT2VCG(e->button (), e->modifiers ()));
    TwMouseReleaseQt(this,e);
    if (user_is_picking)
    {
        user_is_picking = false;
        //TPath.AddSharpConstraints(GPath.PickedPoints);
        //TPath.AddSharpConstraints();
        if (GPath.PickedPoints.back().size() < 1){
            GPath.PickedPoints.pop_back();
            return;
        }
        DoBatchProcess();
    }
    updateGL ();
}

void MyGLWidget::wheelEvent (QWheelEvent * e)
{
    const int WHEEL_STEP = 120;
    track.MouseWheel (e->delta () / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
    updateGL ();
}
