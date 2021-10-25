#ifndef PATH_GL
#define PATH_GL

#include <vcg/space/distance3.h>
#include <vcg/space/segment3.h>
#include <vcg/space/point2.h>
#include <wrap/gl/picking.h>
#include <wrap/gl/trimesh.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/algorithms/closest.h>

template <class TriMeshType>
class PathGL
{
    typedef typename TriMeshType::FaceType FaceType;
    typedef typename TriMeshType::VertexType VertexType;
    typedef typename TriMeshType::CoordType CoordType;
    typedef typename TriMeshType::ScalarType ScalarType;

    int getClosestPolyline(const CoordType &pos,ScalarType maxD)
    {
        int closestI=-1;
        ScalarType currD=std::numeric_limits<ScalarType>::max();
        for (size_t i=0;i<PickedPoints.size();i++)
            for (size_t j=0;j<PickedPoints[i].size()-1;j++)
            {
                CoordType P0=PickedPoints[i][j];
                CoordType P1=PickedPoints[i][j+1];
                vcg::Segment3<ScalarType> S(P0,P1);
                CoordType clos;
                ScalarType distTest;
                vcg::SegmentPointDistance(S,pos,clos,distTest);
                if (distTest>maxD)continue;
                if (distTest>currD)continue;
                currD=distTest;
                closestI=i;
            }
        //assert(closestI>=0);
        return (closestI);
    }

public:
    std::vector<std::vector<CoordType> > PickedPoints;


    void GlDrawPath(size_t IndexPath)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDisable(GL_LIGHT0);
        glLineWidth(10);
        glDepthRange(0,0.9999);
        vcg::glColor(vcg::Color4b(255,0,0,255));
        //glBegin(GL_LINE_STRIP);
        glBegin(GL_LINE_STRIP);
        for (size_t j=0;j<PickedPoints[IndexPath].size();j++)
            vcg::glVertex(PickedPoints[IndexPath][j]);
        glEnd();

        glPopAttrib();
    }

    void GlDrawPaths()
    {
        if (PickedPoints.size()==0)return;

        for (size_t i=0;i<PickedPoints.size();i++)
            GlDrawPath(i);
    }

    void GlDrawLastPath()
    {
        if (PickedPoints.size()==0)return;
        GlDrawPath(PickedPoints.size()-1);
    }

    void ClearPath()
    {
        PickedPoints.resize(PickedPoints.size()+1);
        //Pixels.resize(Pixels.size()+1);
    }

    void AddNewPath()
    {
        PickedPoints.resize(PickedPoints.size()+1);
    }

    void GLAddPoint(const vcg::Point2i &pixel)
    {
        CoordType pp;
        if(vcg::Pick<CoordType>(pixel.X(),pixel.Y(),pp))
            PickedPoints.back().push_back(pp);
    }

    bool GLRemovePathFromPoint(const vcg::Point2i &pixel,ScalarType maxD)
    {
        CoordType pp;
        if(!vcg::Pick<CoordType>(pixel.X(),pixel.Y(),pp))return false;
        if (PickedPoints.size()==0)return false;
        int closest=getClosestPolyline(pp,maxD);
        if (closest==-1)return false;
        PickedPoints.erase(PickedPoints.begin() + closest );
        return true;
    }

    PathGL(){}//TriMeshType &_mesh):mesh(_mesh){}

};

#endif
