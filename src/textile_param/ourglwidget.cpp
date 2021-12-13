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
#include "ourglwidget.h" // must be included before vcg/wrap stuff ?

#include "wrap/igl/arap_parametrization.h"
#include "fibre_elongation.h"
#include "triangle_mesh_type.h"

#include <wrap/qt/trackball.h>
#include <wrap/qt/anttweakbarMapper.h>
#include <wrap/io_trimesh/import_field.h>
#include <wrap/io_trimesh/export_field.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/gl/trimesh.h>

std::string pathM="";
//std::string pathH="";

vcg::Trackball track;//the active manipulator

MyTriMesh tri_mesh;

TwBar *barRem;

vcg::GlTrimesh<MyTriMesh> glWrap;
vcg::GlTrimesh<MyTriMesh> glWrapFibre;

typedef typename MyTriMesh::ScalarType ScalarType;
typedef typename MyTriMesh::CoordType CoordType;

//bool drawTris=true;

enum DrawMode{Draw3D,DrawUV};
DrawMode DMode=Draw3D;
DrawMode OldDMode=Draw3D;

FibreElongation <MyTriMesh> FElong(tri_mesh);
MyTriMesh Fibre3D,FibreUV;
bool drawFibres=false;
size_t Fibre_res=50;

MyTriMesh::ScalarType maxCompr,maxElong;

void UpdateCoords()
{
    if (DMode==Draw3D)
        tri_mesh.RestoreRestPos();
    else
        tri_mesh.SetUVPosition();

    tri_mesh.UpdateDataStructures();
}

void UpdateFibreMesh()
{
    tri_mesh.RestoreRestPos();
    FElong.ExtractFiberMeshes(Fibre3D,FibreUV,maxElong,maxCompr);
    Fibre3D.UpdateDataStructures();
    FibreUV.UpdateDataStructures();
    UpdateCoords();
}

void TW_CALL InitialParametrize(void *)
{
    //HMesh.ReprojectOnMesh(tri_mesh);
    tri_mesh.RestoreRestPos();
    vcg::tri::InitializeArapWithLSCM<MyTriMesh>(tri_mesh);
    vcg::tri::OptimizeUV_ARAP(tri_mesh);
    //bool Loaded=tri_mesh.LoadTriMesh("../../quad_dist.obj");
//    vcg::tri::OptimizeUV_LSCM(tri_mesh);
    //vcg::tri::OptimizeUV_ARAP(tri_mesh);
    UpdateCoords();
}

void TW_CALL InitFibres(void *)
{
    tri_mesh.RestoreRestPos();
    FElong.ComputeFibreElongation(Fibre_res);
    maxCompr=FElong.MaxCompressR;
    maxElong=FElong.MaxElongR;
    UpdateCoords();
    UpdateFibreMesh();
}

void TW_CALL FibreParametrize(void *)
{
    tri_mesh.RestoreRestPos();
    FElong.FibreOptimization(Fibre_res,1,0.2);
    UpdateFibreMesh();
    UpdateCoords();
}

void SetFieldBarSizePosition(QWidget *w)
{
    int params[2];
    params[0] = QTDeviceWidth(w) / 3;
    params[1] = QTDeviceHeight(w) / 1.8;
    TwSetParam(barRem, NULL, "size", TW_PARAM_INT32, 2, params);
    params[0] = QTLogicalToDevice(w, 10);
    params[1] = 30;//QTDeviceHeight(w) - params[1] - QTLogicalToDevice(w, 10);
    TwSetParam(barRem, NULL, "position", TW_PARAM_INT32, 2, params);
}

void InitBar(QWidget *w)
{
    barRem = TwNewBar("HexMeshSimple");

    SetFieldBarSizePosition(w);


    TwEnumVal drawmode[2] = { {Draw3D, "Draw 3D"},
                              {DrawUV, "Draw UV"}
                            };

    TwType drawMode = TwDefineEnum("DrawMode", drawmode, 2);
    TwAddVarRW(barRem, "Draw Mode", drawMode, &DMode, " help='Change draw mode.' ");

    TwAddVarRW(barRem,"DrawFibres",TW_TYPE_BOOL8, &drawFibres," label='Draw Fibres'");
    //    TwAddVarRW(barRem,"DrawHexa",TW_TYPE_BOOL8, &drawHex," label='Draw Hexa'");

    //    TwAddButton(barRem,"Reload",Reload,0,"label='Reload'");

    //    TwAddVarRW(barRem,"Iterations",TW_TYPE_INT32, &Iterations," label='Iterations'");
    //    TwAddVarRW(barRem,"Edge Size",TW_TYPE_DOUBLE, &EdgeStep," label='Edge Size'");

    //    TwAddVarRW(barRem,"Multiplier",TW_TYPE_DOUBLE, &Multiplier," label='Multiplier'");
    //    TwAddButton(barRem,"Average",AverageEdge,0,"label='Set Average Edge'");
    //    //TwAddButton(barQuad,"Optimal",OptimalEdge,0,"label='Set Optimal Edge'");

    //    TwAddVarRW(barRem,"Do Collapse",TW_TYPE_BOOL8, &do_collapse," label='Do Collapse'");
    //    TwAddVarRW(barRem,"Do Smooth",TW_TYPE_BOOL8, &do_smooth," label='Do Smooth'");
    //    TwAddVarRW(barRem,"Do Project",TW_TYPE_BOOL8, &do_project," label='Do Project'");

    //    TwAddVarRW(barRem,"SharpDegreeRem",TW_TYPE_DOUBLE, &SharpDegreeRemesh," label='Sharp Degree'");
    //    TwAddButton(barRem,"Remesh",Remesh,0,"label='Remesh'");

    //    TwAddVarRW(barRem,"SharpDegreeVis",TW_TYPE_DOUBLE, &SharpDegreeVisual," label='Sharp Degree Smooth'");
    //    TwAddButton(barRem,"SetSharp",InitSharpFeatures,0,"label='InitSharp'");
    //    TwAddVarRW(barRem,"Smooth Steps",TW_TYPE_INT32, &SmoothSteps," label='Smooth Steps'");
    //    TwAddButton(barRem,"Smooth",SmoothReproject,0,"label='Smooth'");


    //    TwAddVarRW(barRem,"Flip",TW_TYPE_BOOL8, &flip," label='Flip'");
    TwAddButton(barRem,"Initial Parametrize",InitialParametrize,0,"label='Init Parametrization'");
    TwAddVarRW(barRem,"FibreRes",TW_TYPE_INT32, &Fibre_res," label='Fibre Resolution'");
    TwAddButton(barRem,"Init Fibres",InitFibres,0,"label='Init Fibres'");

    TwEnumVal lenghtmode[2] = { {LModeLocal, "LMode Local"},
                                {LModeGLobal, "LMode Global"}
                                };

    TwType lenghtMode = TwDefineEnum("LenghtMode", lenghtmode, 2);
    TwAddVarRW(barRem, "Lenght Mode", lenghtMode, &FElong.LMode, " help='Change length mode.' ");

//    TwAddVarRW(barRem,"DrawFibres",TW_TYPE_BOOL8, &drawFibres," label='Draw Fibres'");

    TwAddButton(barRem,"Fibre Param",FibreParametrize,0,"label='Fibre Parametrize'");


    //TwAddButton(barRem,"Save",SaveHexMesh,0,"label='Save Hex'");
}


OurGLWidget::OurGLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{

    bool Loaded=tri_mesh.LoadTriMesh(pathM);
    if (!Loaded)
    {
        std::cout<<"Error Loading Mesh"<<std::endl;
        exit(0);
    }
    std::cout<<"Loaded "<<tri_mesh.face.size()<<" faces "<<std::endl;
    std::cout<<"Loaded "<<tri_mesh.vert.size()<<" vertices "<<std::endl;

    glWrap.m=&tri_mesh;

    tri_mesh.UpdateDataStructures();


    //HMesh.Load(pathH.c_str());
}



void OurGLWidget::initializeGL ()
{
    //initialize Glew
    glewInit();
    glClearColor(255, 0, 0, 0);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
}


void OurGLWidget::resizeGL (int w, int h)
{
    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    TwWindowSize(w, h);
    InitBar(this);
    initializeGL();
}

void OurGLWidget::paintGL ()
{
    if (OldDMode!=DMode)
    {
        UpdateCoords();
        OldDMode=DMode;
    }
    glClearColor(0,255,255,255);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, OurGLWidget::width()/(float)OurGLWidget::height(), 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,3.5f,   0,0,0,   0,1,0);
    track.center=vcg::Point3f(0, 0, 0);
    track.radius= 1;
    track.GetView();

    glPushMatrix();
    track.Apply();
    glPushMatrix();

    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    vcg::glScale(2.0f/tri_mesh.bbox.Diag());
    glTranslate(-tri_mesh.bbox.Center());

    //    if (drawTris)
    //    {
    //tri_mesh.GLDrawSharpEdges();
    vcg::glColor(vcg::Color4b(230,230,230,255));
    glWrap.Draw(vcg::GLW::DMSmooth,vcg::GLW::CMPerVert,vcg::GLW::TMNone);

    if (drawFibres)
    {
        if (DMode==Draw3D)
            glWrapFibre.m= &Fibre3D;
        else
            glWrapFibre.m= &FibreUV;

        glWrapFibre.Draw(vcg::GLW::DMSmooth,vcg::GLW::CMPerFace,vcg::GLW::TMNone);
    }
    // }
    /*
    if (drawHex)
        HMesh.GLDraw();*/

    glPopMatrix();
    glPopMatrix();


    TwDraw();

}

void OurGLWidget::keyReleaseEvent (QKeyEvent * e)
{
    e->ignore ();
    if (e->key () == Qt::Key_Control)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ControlModifier));
    if (e->key () == Qt::Key_Shift)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
    if (e->key () == Qt::Key_Alt) track.ButtonUp (QT2VCG (Qt::NoButton, Qt::AltModifier));
    updateGL ();
}


void OurGLWidget::keyPressEvent (QKeyEvent * e)
{
    e->ignore ();
    if (e->key () == Qt::Key_Control) track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ControlModifier));
    if (e->key () == Qt::Key_Shift)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
    if (e->key () == Qt::Key_Alt)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::AltModifier));

    TwKeyPressQt(e);
    updateGL ();
}

void OurGLWidget::mousePressEvent (QMouseEvent * e)
{
    if(!TwMousePressQt(this,e))
    {
        e->accept ();
        setFocus ();
        track.MouseDown(QT2VCG_X(this, e), QT2VCG_Y(this, e), QT2VCG (e->button (), e->modifiers ()));
    }
    updateGL ();
}

void OurGLWidget::mouseMoveEvent (QMouseEvent * e)
{
    if (e->buttons ()) {
        track.MouseMove(QT2VCG_X(this, e), QT2VCG_Y(this, e));
        updateGL ();
    }
    TwMouseMotion(QTLogicalToDevice(this, e->x()), QTLogicalToDevice(this, e->y()));
}

void OurGLWidget::mouseDoubleClickEvent (QMouseEvent * e)
{
    //    if (e->buttons ())
    //    {
    //        xMouse=QT2VCG_X(this, e);
    //        yMouse=QT2VCG_Y(this, e);
    //        //pointToPick=Point2i(e->x(),height()-e->y());
    //        //pointToPick=Point2i(xMouse,yMouse);
    //        hasToPick=true;
    //        updateGL ();
    //    }
}


void OurGLWidget::mouseReleaseEvent (QMouseEvent * e)
{
    track.MouseUp(QT2VCG_X(this, e), QT2VCG_Y(this, e), QT2VCG(e->button (), e->modifiers ()));
    TwMouseReleaseQt(this,e);
    updateGL ();
}

void OurGLWidget::wheelEvent (QWheelEvent * e)
{
    const int WHEEL_STEP = 120;
    track.MouseWheel (e->delta () / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
    updateGL ();
}
