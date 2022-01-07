
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>

#include "param/cloth_param.h"

int main(int argc, char *argv[]){
    Eigen::MatrixXd V_3d, V_2d;
    Eigen::MatrixXi F;
    igl::readOBJ("../data/dress_front.obj", V_3d, F);


    ClothParam cloth(V_3d, F, 0.05);
    bool success = cloth.paramAttempt();


    V_2d = cloth.getV2d();

    // ------------- VIEWER ------------- //

    bool show_uv = false;

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V_3d, F);
    viewer.data().set_uv(10 * V_2d);
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
}