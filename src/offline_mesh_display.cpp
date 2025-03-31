//
// Created by neo on 2024/9/27.
//
#include <mutex>
#include <thread>
#include <cmath>
#include <math.h>
#include <unistd.h>
#include <fstream>
#include <stdio.h>
#include <csignal>
//#include <Eigen/Core>
//#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <condition_variable>
// opencv
#include <opencv2/opencv.hpp>
#include <cv_bridge/cv_bridge.h>
// ros
#include <ros/ros.h>
// pcl
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>

#include "GL_gui.h"


int main(int argc, char **argv) {
    ros::init(argc, argv, "offline_mesh_display");
    ros::NodeHandle m_node;

    // 参数
    GL_gui m_gui;
    std::string ply_file_path;
    m_node.param<string>("input_path", ply_file_path, "xxxxx");

    // 读取离线mesh
    pcl::PolygonMesh poly_mesh;
    if (pcl::io::loadPLYFile(ply_file_path, poly_mesh) == -1) {
        std::cout << "Couldn't read file: " << ply_file_path << std::endl;
        return -2;
    }
    PointCloudXYZI::Ptr point_cloud = boost::make_shared<PointCloudXYZI>();
    pcl::fromPCLPointCloud2(poly_mesh.cloud, *point_cloud);

    // 可视化
    for (auto it : poly_mesh.polygons)
    {
        // 新建Triangle对象
        int temp_tri_pts_id[ 3 ];
        temp_tri_pts_id[0] = it.vertices[0];
        temp_tri_pts_id[1] = it.vertices[1];
        temp_tri_pts_id[2] = it.vertices[2];
        Triangle_ptr triangle_ptr;
        triangle_ptr = std::make_shared< Triangle >( temp_tri_pts_id[0], temp_tri_pts_id[1], temp_tri_pts_id[2] );
        triangle_ptr->m_index_flip = m_gui.if_index_flip( temp_tri_pts_id[0], temp_tri_pts_id[1], temp_tri_pts_id[2] );
        triangle_ptr->rec_source = 0;

        // 插入到体素中
        Eigen::Vector3d triangle_pos = m_gui.get_tri_center(temp_tri_pts_id[0], temp_tri_pts_id[1], temp_tri_pts_id[2], 0);
        int          hash_3d_x = std::round( triangle_pos( 0 ) / m_gui.g_triangles_manager.m_region_size );     // 根据三角形重心点确定所属体素
        int          hash_3d_y = std::round( triangle_pos( 1 ) / m_gui.g_triangles_manager.m_region_size );
        int          hash_3d_z = std::round( triangle_pos( 2 ) / m_gui.g_triangles_manager.m_region_size );
        Sync_triangle_set* sync_triangle_set_ptr = g_triangles_manager.m_triangle_set_in_region.get_data( hash_3d_x, hash_3d_y, hash_3d_z );
        int temp_set_idx;
        if ( sync_triangle_set_ptr == nullptr )
        {
            sync_triangle_set_ptr = new Sync_triangle_set();
            sync_triangle_set_ptr->insert( triangle_ptr );
            sync_triangle_set_ptr->vec_idx = g_triangles_manager.m_triangle_set_vector.size();

            m_gui.g_triangles_manager.m_triangle_set_in_region.insert( hash_3d_x, hash_3d_y, hash_3d_z, *sync_triangle_set_ptr );
            m_gui.g_triangles_manager.m_triangle_set_vector.push_back( m_gui.g_triangles_manager.m_triangle_set_in_region.get_data( hash_3d_x, hash_3d_y, hash_3d_z ) );
        }
        else
        {
            sync_triangle_set_ptr->insert(triangle_ptr);
        }
    }

    m_gui.service_refresh_and_synchronize_triangle();



























}
