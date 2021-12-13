#include <QApplication>
#include <QDesktopWidget>
#include <GL/glew.h>
#include <wrap/qt/anttweakbarMapper.h>
#include <QWindow>
#include <QTextStream>
#include<vcg/complex/algorithms/hole.h>

#include "myglwidget.h"

extern std::string pathDef,pathRef,pathFrames;

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QWindow dummy;
    QString def_string = QString("GLOBAL fontscaling=%1").arg((int)dummy.devicePixelRatio());
    TwDefine(def_string.toStdString().c_str());
    printf("%s\n",qPrintable(def_string));
    fflush(stdout);

    // Set functions to handle string copy
    TwCopyCDStringToClientFunc(CopyCDStringToClient);
    TwCopyStdStringToClientFunc(CopyStdStringToClient);

    if( !TwInit(TW_OPENGL, NULL) )
    {
        fprintf(stderr, "AntTweakBar initialization failed: %s\n", TwGetLastError());
        return 1;
    }

    assert(argc>2);

    pathDef=std::string(argv[1]);
    std::cout <<"Loading Deformed Mesh "<< pathDef << std::endl;

    pathRef=std::string(argv[2]);
    std::cout <<"Loading Reference Mesh "<< pathRef << std::endl;

    if (argc>3)
    {
        pathFrames=std::string(argv[3]);
        std::cout <<"Loading Frames "<< pathFrames << std::endl;
    }

    MyGLWidget window;

    window.show();
    return app.exec();
}
