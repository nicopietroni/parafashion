############################ PROJECT FILES ############################

include(libs.pri)

SOURCES = \
    src/myglwidget.cpp \
    apps/parafashion.cpp

DEFINES += GLEW_STATIC


############################ TARGET ############################

#App config
TARGET = parafashion

TEMPLATE = app
CONFIG += qt
CONFIG += c++11
CONFIG -= app_bundle

QT += core gui opengl xml widgets
QT += svg

#Debug/release optimization flags
CONFIG(debug, debug|release){
    DEFINES += DEBUG
}
CONFIG(release, debug|release){
    DEFINES -= DEBUG
    #just uncomment next line if you want to ignore asserts and got a more optimized binary
    CONFIG += FINAL_RELEASE
}

#Final release optimization flag
FINAL_RELEASE {
    unix:!macx{
        QMAKE_CXXFLAGS_RELEASE -= -g -O2
        QMAKE_CXXFLAGS += -O3 -DNDEBUG
    }
}

macx {
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.13
    QMAKE_MAC_SDK = macosx10.13
}


############################ INCLUDES ############################

#vcglib
INCLUDEPATH += $$VCGLIBPATH
SOURCES += $$VCGLIBPATH/wrap/ply/plylib.cpp
SOURCES += $$VCGLIBPATH/wrap/gui/trackball.cpp
SOURCES += $$VCGLIBPATH/wrap/gui/trackmode.cpp
SOURCES += $$VCGLIBPATH/wrap/qt/anttweakbarMapperNew.cpp
SOURCES += lib/Clipper/clipper.cpp

#eigen
INCLUDEPATH += $$EIGENPATH
INCLUDEPATH += $$TRACINGPATH
INCLUDEPATH += $$GLEWPATH/include
INCLUDEPATH += $$CLIPPERPATH

#libigl
INCLUDEPATH += $$LIBIGLPATH/include
HEADERS += include/symmetrizer.h \
    include/parametrizer.h \
    include/svg_exporter.h \
    include/field_computation.h \
    include/parafashion_interface.h \
    include/animation_manager.h \
    include/trace_path_ui.h \
    include/myglwidget.h \
    include/parafashion.h \
    include/smooth_field_directional.h \
    $$LIBIGLPATH/include/igl/principal_curvature.h \
    $$LIBIGLPATH/include/igl/copyleft/comiso/nrosy.h \
    $$VCGLIBPATH/wrap/qt/Outline2ToQImage.h \

SOURCES += \
    $$LIBIGLPATH/include/igl/principal_curvature.cpp \
    $$LIBIGLPATH/include/igl/copyleft/comiso/nrosy.cpp \
    $$VCGLIBPATH/wrap/qt/Outline2ToQImage.cpp \
    src/animation_manager.cpp \
    src/field_computation.cpp \

#AntTweakBar
INCLUDEPATH += $$ANTTWEAKBARPATH/include
LIBS += -L$$ANTTWEAKBARPATH/lib -lAntTweakBar
win32{ # Awful problem with windows..
    DEFINES += NOMINMAX
}

#glew
#LIBS += -lGLU
INCLUDEPATH += ./include/
INCLUDEPATH += $$GLEWPATH/include
INCLUDEPATH += $$DIRECTIONALPATH
SOURCES += $$GLEWPATH/src/glew.c

# Mac specific Config required to avoid to make application bundles
macx{
    QMAKE_POST_LINK +="cp -P ../../../code/lib/AntTweakBar1.16/lib/libAntTweakBar.dylib . ; "
    QMAKE_POST_LINK +="install_name_tool -change ../lib/libAntTweakBar.dylib ./libAntTweakBar.dylib $$TARGET ; "
    DEPENDPATH += .
}
