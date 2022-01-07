
#include <igl/readOBJ.h>

//#define ENABLE_PARAM_VIEWER
#ifdef ENABLE_PARAM_VIEWER
#include <igl/opengl/glfw/Viewer.h>
#endif

#include "param/cloth_param.h"

int main(int argc, char *argv[]){

    // example input: ./param ../data/mark_skirt/mark_skirt_back_left_cut.obj

    std::string input_path = "../data/dress_front.obj";
    if (argc > 1) input_path = std::string(argv[1]);

    Eigen::MatrixXd V_3d, V_2d;
    Eigen::MatrixXi F;
    igl::readOBJ(input_path, V_3d, F);


    ClothParam cloth(V_3d, F, 0.05);
    bool success = cloth.paramAttempt();

    std::cout << "Success: " << success << std::endl;
    cloth.printStretchStats();
    V_2d = cloth.getV2d();


    // ------------- VIEWER ------------- //
    #ifdef ENABLE_PARAM_VIEWER
    bool show_uv = false;

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V_3d, F);
    viewer.data().set_uv(10 * V_2d / V_2d.maxCoeff());
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
        if (key == '1')
            show_uv = false;
        else if (key == '2')
            show_uv = true;

        if (show_uv){
            viewer.data().set_mesh(V_2d, F);
            viewer.core().align_camera_center(V_2d, F);
        }
        else
        {
            viewer.data().set_mesh(V_3d, F);
            viewer.core().align_camera_center(V_3d, F);
        }

        viewer.data().compute_normals();

        return false;
    };

    viewer.data().show_lines = false;
    viewer.data().show_texture = true;
    viewer.launch();
    #endif
}