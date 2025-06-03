//
// Created by neo on 2024/7/10.
//

#ifndef SRC_GL_GUI_H
#define SRC_GL_GUI_H

#include <iostream>
#include <stdio.h>
#include <Eigen/Core>
#include <csignal>
#include <fstream>
#include <math.h>
#include <mutex>
#include <omp.h>
#include <so3_math.h>
#include <unistd.h>
#include <chrono>
#include <thread>

#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <opencv2/opencv.hpp>

#include "tools/tools_color_printf.hpp"
#include "tools/tools_data_io.hpp"
#include "tools/tools_logger.hpp"
#include "tools/tools_color_printf.hpp"
#include "tools/tools_eigen.hpp"
#include "tools/tools_random.hpp"
#include "tools/lib_sophus/so3.hpp"
#include "tools/lib_sophus/se3.hpp"
#include "tools/tools_serialization.hpp"
#include "tools/tools_graphics.hpp"
#include "tools/openGL_libs/glad.h"
#include "tools/openGL_libs/openGL_camera.hpp"
#include "tools/openGL_libs/gl_draw_founction.hpp"

#include "shader/tools_my_texture_triangle_shader.h"
#include "shader/tools_my_camera_pose_shader.h"
#include "shader/tools_my_dbg_utility.h"
#include "shader/tools_my_point_shader.h"

#include "tools_timer.hpp"
#include "tinycolormap.hpp"
#include "triangle.hpp"
#include "ikd_Tree.h"
//#include "mesh_rec_display.hpp"
#include "region_tri_shader.hpp"
//#include "MeshFragment.h"
//#include "NonPlaneMesh.h"


#define NUMBER_OF_POSE_SIZE -1
typedef std::vector< std::pair< std::vector< vec_4 >, Eigen::Matrix< double, NUMBER_OF_POSE_SIZE, 1 > > > LiDAR_frame_pts_and_pose_vec;

class GL_gui
{
public:
//    double g_maximum_pe_error = 40;
//    double g_initial_camera_exp_tim = 1.0;

//    double              g_max_incidence_angle = 90;
//    Common_tools::Timer g_cost_time_logger;

//    Eigen::Matrix3d g_camera_K;

    std::string data_path_file = std::string( Common_tools::get_home_folder() ).append( "/ImMesh_output/" );

    int    appending_pts_frame = ( int ) 5e3;
    double threshold_scale = 1.0; // normal
    double region_size = 10.0;    // normal

    double minimum_pts = 0.1 * threshold_scale;
    double g_meshing_voxel_size = 0.4 * threshold_scale;

    GL_camera g_gl_camera;

//    bool                        g_if_automatic_save_mesh = false;
//    float                       g_alpha = 1.0;
//    std::vector< vec_3 > pts_of_maps;

//    bool             show_background_color = false;
    float            camera_focus = 2000;

    Triangle_manager g_triangles_manager;


//    std::vector< GLuint >      g_texture_id_vec;
    long                       img_width = 0, img_heigh = 0;

    LiDAR_frame_pts_and_pose_vec g_eigen_vec_vec;
    bool                         g_flag_pause = false;
//    int                          if_first_call = 1;

//    std::string          g_debug_string;
    int                  g_current_frame = -1;

// GUI settting
    bool   g_display_mesh = true;
    int    g_enable_mesh_rec = true;
    int    g_save_to_offline_bin = false;
    int    g_display_face = 1;
    bool   g_draw_LiDAR_point = true;
    float  g_draw_path_size = 2.0;
    float  g_display_camera_size = 1.0;
    float  g_ply_smooth_factor = 1.0;
    int    g_ply_smooth_k = 20.0;
    bool   g_display_main_window = false;
    bool   g_display_camera_pose_window = false;
    bool   g_display_help_win = false;
    bool   g_follow_cam = true;
    bool   g_mesh_if_color = false;
    bool   g_if_draw_z_plane = true;
    bool   g_if_draw_wireframe = false;
    bool   g_if_draw_depth = false;
    bool   g_if_depth_bind_cam = true;

    bool   g_force_refresh_triangle = false;

    ImVec4 g_mesh_color = ImVec4( 1.0, 1.0, 1.0, 1.0 );

    // mesh_rec_display
    Common_tools::Point_cloud_shader    g_path_shader;
    Common_tools::Point_cloud_shader    g_LiDAR_point_shader;
    Common_tools::Triangle_facet_shader g_triangle_facet_shader;
    Common_tools::Axis_shader           g_axis_shader;
    Common_tools::Ground_plane_shader   g_ground_plane_shader;
    Common_tools::Camera_pose_shader    g_camera_pose_shader;

//    vec_3f        g_axis_min_max[ 2 ];
    std::mutex                                                                  g_region_triangle_shader_mutex;
    std::vector< std::shared_ptr< Region_triangles_shader > >                   g_region_triangles_shader_vec;
    std::map< Sync_triangle_set *, std::shared_ptr< Region_triangles_shader > > g_map_region_triangles_shader;
    std::vector< vec_3 > pt_camera_traj;

    // 三角形idx对应的指针
    std::unordered_map<int, Triangle_ptr>               m_triangle_hash_np;
    std::vector<std::unordered_map<int, Triangle_ptr>>  m_triangle_hash_p;
    std::unordered_map<int, int>                        m_tri_set_np;
    std::vector<std::unordered_map<int, int>>           m_tri_set_p;

    GL_gui()
    {
        g_eigen_vec_vec.resize( 1e6 );
        m_triangle_hash_p.resize( 1e6 );
        m_tri_set_p.resize( 1e6 );
    };
    ~GL_gui() = default;

    void print_help_window( bool *if_display_help_win );
    void get_last_avr_pose( int current_frame_idx, Eigen::Quaterniond &q_avr, vec_3 &t_vec );

    void display_current_LiDAR_pts( int current_frame_idx, double pts_size, vec_4f color = vec_4f( 1.0, 1.0, 1.0, 0.5 ) );
    void display_reinforced_LiDAR_pts(  std::vector< vec_3f > & pt_vec, double pts_size, vec_3f color = vec_3f (1.0, 0.0, 1.0) );
    void init_openGL_shader();
    void draw_triangle( const Cam_view &  gl_cam);
    void display_camera_traj( float display_size );
    void draw_camera_trajectory(int current_frame_idx, float pt_disp_size = 5);
    void draw_camera_pose( int current_frame_idx, float pt_disp_size = 5, float display_cam_size = 1.0);

    void display();

    void synchronize_triangle_list_for_disp();
    void service_refresh_and_synchronize_triangle();

    Eigen::Vector3d get_tri_center(int A_idx, int B_idx, int C_idx, int idx);
    int if_index_flip(int A_idx, int B_idx, int C_idx);
    void triangle_manager_update(const std::vector<int>& update_plane_idx_list);







};


#endif //SRC_GL_GUI_H
