//
// Created by neo on 2024/2/29.
// 定义mesh片段类
//

#ifndef SRC_MESHFRAGMENT_H
#define SRC_MESHFRAGMENT_H

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
//
#include "DataType.h"
#include "unordered_map"
#include "QuadTree.hpp"

// 与坐标轴平行的包围框
struct AABB {
    Eigen::Vector3d center;   // 中心点
    Eigen::Vector3d max_pt;   // xyz轴坐标系下包围盒最大坐标点（对应特征值从小到大顺序）
    Eigen::Vector3d min_pt;   // 最小坐标点
};


class MeshFragment
{
public:
    PointCloudXYZI::Ptr                         ptcl_all;           // 所有点云-->世界坐标系下
    PointCloudXYZI::Ptr                         ptcl_grid;          // grid_decimate后的网格顶点
    std::vector<std::array<int, 3>>             faces_list;         // 所有三角形列表，包含三个点的索引
    std::vector<std::vector<int>>               pt2face_idx;        // 每个点对应的三角形列表
    std::vector<int>                            faces_to_delete;    // 最后要删除的三角形列表
    std::vector<int>                            pts_state;          // 点的状态 1-内部 2-边缘 3-删除

    Eigen::Vector3d                             center_point;       // 点云簇中心点
    Eigen::Vector3d                             normal_vector;      // 平面法向
    AABB                                        rec_index;          // 有向包围盒

    double                                      plane_quality;                   // 最小特征值，筛选平面用
    double                                      ptsNum_thres = 1000;             // 点数量超过这个阈值就降采样
    bool                                        if_first_decimate = true;        // 是否第一次抽稀
    Eigen::Matrix3d                             quadtree_axes;                   // 投影旋转矩阵
    Eigen::Vector3d                             quadtree_center;                 // 投影中心点
    int                                         plane_update_num = 0;            // 平面更新次数
    int                                         last_grid_decimate_num = 0;      // 上次降采样平面点的数量，用于确定下次降采样的起点
    int                                         grid_ptsNum_threshold = 4;       // 网格内点数量少于该阈值则判断该区域内没有点
    double                                      yMax, yMin, zMax, zMin;          // 投影点的范围
    double                                      y0, z0;                          // grid中心原点坐标
    double                                      y_start, z_start;                // grid左上角位置
    int                                         topH, bottomH, leftW, rightW;    // 以grid中心为原点分块
    double                                      dis_resolution;                  // cell边长
    std::deque<std::deque<grid>>                uniform_grid;                    // grid
    std::deque<std::deque<int>>                 tri_vertex_idx_grid;             // 记录grid顶点索引
    std::vector<std::array<int, 2>>             grid_need_update_list;           // 需要更新mesh的cell列表

    // 非平面相关
    std::vector<int>                            faces_to_delete_frame;           // 删除的三角形索引
    std::vector<int>                            pts_edge_delete;                 // 从删除的三角形中提取的不重复的边缘点
    std::vector<int>                            pts_edge_add;                    // 从添加的边缘三角形中提取的不重复点，包括边缘点和内部点
    std::vector<int>                            pts_edge_update;                 // cell状态2 -> 2，点坐标需更新的点索引
    std::vector<int>                            pts_inner_add;                   // 新增的内部点，不重复
    std::vector<int>                            faces_with_edge_vertex;          // 新添加的包含边缘点的三角形
    std::vector<int>                            faces_with_2_edge_vertex;
    std::unordered_map<std::string, int>        tri_delete_string_idx;           // 删除的三角形索引，例：三个顶点索引 16 21 15 -》 "15 16 21"
    std::unordered_map<int, int>                plane_vidx_to_noplnae_vidx;      // 点在平面map中的索引 -》 点在非平面map 中的索引
    int                                         plane_idx = 0;                   // 记录平面索引

    // gui相关
    int                                                     last_faces_list_num = 0;         // 上次拷贝到shader时总的三角形数量
    int                                                     last_faces_to_delete_num = 0;    // 上次拷贝到shader时总的要删除的三角形数量
    bool                                                    if_delete = false;               // 该平面是否要删除
    std::mutex                                              mtx_ptcl_grid;



    MeshFragment()
    {
        // 初始化
        ptcl_all = boost::make_shared<PointCloudXYZI>();
        ptcl_grid = boost::make_shared<PointCloudXYZI>();
    };
    ~MeshFragment() {};


    // 赋值：平面点，pcl/eigen
    void give_plane_point(const Eigen::Tensor<double, 3>& project_image, const std::vector<Eigen::Vector2i>& point_uv,
                            Eigen::Matrix3d rot_mat, Eigen::Vector3d pos_vec);
    // 赋值：平面属性，中心、法向、最小包围框（带方向的）
    void compute_plane_parameter(const Eigen::Vector3d& frame_pos);
    // 赋值：平面属性，中心、法向、最小包围框（带方向的）
    void update_plane(const std::shared_ptr<MeshFragment> new_plane);
    // 四叉树原理降采样
    void quadtree_decimate();
    // 初始化一个cell
    grid build_new_grid(double y_start, double z_start, int i, int j);
    // 均匀网格降采样
    void grid_decimate(double dis_resolution, FILE *fpTime_mesh);
    // 确定每个cell的vertices
    void MC_mesh();
    void MC_mesh_fast();
    // 根据每个cell的vertices提取vertex和facet
    void vertex_and_face_list_extract();
    void vertex_and_face_list_update();
    void vertex_and_face_list_update_fast();
    // 两个平面的vertex和facet合为一个
    void vertex_and_face_list_merge(const std::shared_ptr<MeshFragment>& new_plane);
    // 保存
    bool save_to_ply(const std::string &fileName);
    // 自定义保存，用于debug
    void mtri_save_to_ply(const std::string &fileName);
    // 保存时删掉重复的vertex
    bool save_to_ply_without_redundancy(const std::string &fileName);
};


#endif //SRC_MESHFRAGMENT_H
