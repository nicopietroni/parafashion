
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>

#include "param/bary_optimizer.h"
#include "param/param_utils.h"

int main(int argc, char *argv[]){

    Eigen::MatrixXd V_3d, V_2d;
    Eigen::MatrixXi F;
    igl::readOBJ("../data/dress_front.obj", V_3d, F);


    V_2d = paramARAP(V_3d, F);

    BaryOptimizer bo(F.rows(), V_3d.rows());

    auto printStretchStats = [&](){
        Eigen::VectorXd stretch_u, stretch_v; // per triangle stretch
        bo.measureScore(V_2d, V_3d, F, stretch_u, stretch_v);
        stretch_u = stretch_u.array() - 1.0;
        stretch_v = stretch_v.array() - 1.0;
        std::cout << "Stretch min -> max (avg): " << std::endl;
        printf("U: %f -> %f (%f)\n", stretch_u.minCoeff(), stretch_u.maxCoeff(), stretch_u.mean());
        printf("V: %f -> %f (%f)\n", stretch_v.minCoeff(), stretch_v.maxCoeff(), stretch_v.mean());
    };


    std::cout << "Before: " << std::endl;
    printStretchStats();

    for (int i=0; i<10; i++){
        V_2d = bo.localGlobal(V_2d, V_3d, F);
    }

    std::cout << "After: " << std::endl;
    printStretchStats();



    // --- VIEWER --- //

    bool show_uv = false;

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V_2d, F);
    viewer.data().set_uv(V_2d);
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

    // Disable wireframe
    viewer.data().show_lines = false;

    // Draw checkerboard texture
    viewer.data().show_texture = true;

    // Launch the viewer
    viewer.launch();
}