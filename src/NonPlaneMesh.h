//
// Created by neo on 2024/4/7.
//

#ifndef SRC_NONPLANEMESHFRAGMENT_H
#define SRC_NONPLANEMESHFRAGMENT_H

// Eigen
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>
// pcl
#include <pcl/point_types.h>
#include <pcl/common/distances.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/io/ply_io.h>
#include <pcl/PolygonMesh.h>
// CGAL
#include <CGAL/Point_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Polygon_2.h>
// #include <CGAL/Boolean_set_operations_2.h>
// https://doc.cgal.org/latest/Kernel_23/index.html#Kernel_23Kernel
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
// 自定义
#include "DataType.h"
#include "hash_idx.hpp"
//#include "iostream"
#include "MeshFragment.h"

class NonPlaneMesh {
public:

    PointCloudXYZI::Ptr                     pts_list;               // 点列表，带法向量
    std::vector<int>                        pts_state;              // 点的状态 0-新加入的点 1-已完成np三角形顶点 2-np未完成/边缘顶点 3-要删除的点 4-平面内部点 5-平面边缘点 6-连接完成的平面边缘点
    std::vector<int>                        pts_process_num;        // 点更新次数
    std::vector<std::vector<int>>           pt2face_idx;            // 每个点对应的三角形列表
    std::vector<std::array<int, 3>>         faces_list_all;         // 所有的三角形列表
    std::vector<int>                        faces_to_delete;        // 最后要删除的三角形列表
    std::vector<int>                        plane_edge_v_list_all;  // 所有平面边缘点的索引，包括已删除的

    Hash_map_3d_base< long, std::shared_ptr<m_voxel> >      voxel_pts;  // 点的体素hash索引

    double                  voxel_resolution;           // 体素边长
    double                  r_dis;                      // 最近邻搜索距离
    double                  minimum_pt_dis;             // 点之间的最小距离，用以处理点云不均匀分布的情况
    int                     N_pts;                      // 拟合新三角形的最少自由点数量

    std::mutex              mtx_edge_v_extract;
    std::mutex              mtx_pts_state, mtx_pts_list, mtx_pts_process_num;
    std::mutex              mtx_pt2face_idx, mtx_faces_list_all, mtx_faces_to_delete;
    std::mutex              mtx_plane_edge_v_list_all, mtx_tri_edge_v_list_all;
    std::mutex              mtx_voxel;

    // mesh_update 三角形质量阈值
    double                  N_r_dis = 1.0;              // A-B-C形式添加三角形时的边长阈值
    double                  ang_thre = 0.866;           // 三角形角度阈值 30du

    // npp_connect相关
    double                  N_r_dis_for_edge_v = 1.0;   // 查找平面边缘点最近的非平面边缘点时的距离，N * r_dis

    // 补洞
    std::vector<int>        hole_v_process;         // 补洞处理流程中边缘点是否处理 0-未处理 1-已处理
    std::vector<int>        hole_process_num;       // 判断边缘点是否是洞的次数
    std::vector<int>        tri_edge_v_list_all;    // 所有的边缘顶点索引列表

    // gui相关
    int                     last_faces_list_all_num = 0;
    int                     last_faces_to_delete_num = 0;


    NonPlaneMesh()
    {
        pts_list = boost::make_shared<PointCloudXYZI>();
    };

    ~NonPlaneMesh() {};


    // 搜索最近邻点
    void search_nearest_vertex_under_distance_r(const int& idx, std::vector<int>& pt_vec, std::vector<double>& dis_vec);
    void search_edge_vertex_under_distance_r(const int& idx, std::vector<int>& vec);
    void search_edge_vertex_under_distance_kr(const int& idx, double max_dis, std::vector<int>& vec);
    void search_nearest_point_under_distance_r(const int& idx, std::vector<int>& vec);
    void search_nearest_valid_ptv(const Eigen::Vector3d& center_pt, double max_dis, std::vector<int>& pt_vec);
    // 判断lidar point沿法向投影在三角形内部还是外部
    bool if_projected_inner(const Eigen::Vector3d& Q, const Eigen::Vector3d& A, const Eigen::Vector3d& B,
                            const Eigen::Vector3d& C,
                            Eigen::Vector3d P = Eigen::Vector3d::Zero() );
    // 判断有无三角形顶点位于三角形内部，处理vertex
    bool if_p_in_tri(const Eigen::Vector3d& Q, const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C);
    // 判断是否是边缘vertex
    bool if_edge_vertex(const int& pt_idx);
    // 找到相邻的边缘vertex
    void find_edge_vertex(const int& pt_idx, std::vector<int> &adj_edge_v);
    // 找到相邻的边缘vertex和错误连接的点
    void find_edge_and_abnormal_vertex(const int& pt_idx, std::vector<int> &adj_edge_v, std::vector<int> &abnormal_v);
    // 判断三角形法向投影是否在边的两侧
    bool if_different_side(const Eigen::Vector3d& edge_A, const Eigen::Vector3d& edge_B, const Eigen::Vector3d& P,
                           const Eigen::Vector3d& Q);
    // 没用到
    bool if_reasonable_normal(const int& A_idx, const int& B_idx, const int& C_idx);
    // 拟合两个三角形
    void fit_2_tri(std::vector<int> &pts_idx_list, Eigen::Vector3d &avg_normal,
                   std::vector<std::array<int, 3>> &tri_pt_indices);
    // 逆时针排序
    void anticlockwise_sort(std::vector<int> &pts_idx_list, int ref_pt_idx, std::vector<std::vector<int>> &tri_v_idx,
                            std::vector<int> &pts_idx_list_sorted, std::vector<int> &pts_axis_idx);
    // 判断是否洞
    bool if_lay_in_hole(int pt_idx, std::vector<int>& adj_v_list, std::vector<int>& all_edge_v);
    // 补洞
    void fill_hole(std::vector<int>& hole_pts_list);
    // 主函数调用的补洞函数
    void mesh_hole_process(FILE *f);
    // 主函数中调用的mesh更新函数
    void mesh_update(const std::vector<std::array<long, 3>>& voxel_to_update,
                     Eigen::Vector3d& body_pos/*, FILE *f*/);
    // 判断点是否投影在平面内部
    bool if_project_inner_plane(const Eigen::Vector3d& pt, const std::shared_ptr<MeshFragment>& plane);
    // 根据平面边缘更新
    void mesh_update_from_plane_map(const std::shared_ptr<MeshFragment>& plane);
    // 根据删除的平面更新
    void mesh_delete_plane_map(const std::shared_ptr<MeshFragment>& plane);
    // 主函数中调用的连接函数
    void plane_and_noplane_connect();
    // 保存
    void save_to_ply(const std::string &fileName);
    // 保存时删掉重复的vertex
    bool save_to_ply_without_redundancy(const std::string &fileName);
    // 四叉树简化后调整平面与非平面中重复的点（连接平面、非平面的点）
    void plane_vertex_adjust(const std::shared_ptr<MeshFragment>& plane);
};


#endif //SRC_NONPLANEMESHFRAGMENT_H
