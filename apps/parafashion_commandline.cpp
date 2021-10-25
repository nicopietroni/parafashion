#include "parafashion_interface.h"

extern std::string pathDef,pathRef;

int main(int argc, char *argv[])
{
    if (argc==1)
    {
        std::cout<<"Pass a Mesh (OBJ or PLY) as an argument"<<std::endl;
        exit(0);
    }
    CMesh mesh;
    std::cout<<"Loading:"<<argv[1]<<std::endl;
    bool loaded=mesh.LoadMesh(std::string(argv[1]));
    assert(loaded);

    std::vector<std::vector<size_t> > FacesReference;
    std::vector<std::vector<double> > CoordReference;
    std::vector<int> Partition;
    std::vector<std::vector<std::vector<double> > > VertUV;
    MeshToVectors<CMesh>(mesh,FacesReference,CoordReference,Partition,VertUV);

    //FOR NOW THE DEFORMED AND THE REFERENCE ARE THE SAME
    std::vector<std::vector<size_t> > FacesDeformed=FacesReference;
    std::vector<std::vector<double> > CoordDeformed=CoordReference;

    //EMPTY PICKED POINTS
    std::vector<std::vector<std::vector<double> > > PickedPoints;

    //OUTPUT
    std::vector<std::vector<size_t> > FacesOutput;
    std::vector<std::vector<double> > CoordOutput;
    std::vector<int> PartitionOutput;
    std::vector<std::vector<std::vector<double> > > VertUVOutput;

    //EXECUTE THE COMPUTATION
    std::cout<<"Deriving Patch Layout"<<std::endl;
    DerivePatchLayout(FacesReference,CoordReference,FacesDeformed,CoordDeformed,
                      PickedPoints,FacesOutput,CoordOutput,PartitionOutput,VertUVOutput);

    std::cout<<"Done"<<std::endl;

    //EXPORTING THE RESULT TO A MESH
    CMesh SaveM;
    VectorsToTriMesh(FacesOutput,CoordOutput,PartitionOutput,VertUVOutput,SaveM);
    SaveM.SaveTriMesh("./test_out.ply");
    return 0;
}
