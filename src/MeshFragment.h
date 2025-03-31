//
// Created by neo on 2024/2/29.
// 定义mesh片段类
//

/*
 * 对于每一帧lidar点云，将其分为平面点和非平面点，接下来就是如何构建mesh
 * 定义一个mesh片段的类，包含平面片段和非平面片段（定义一个flag区分），包含原始点云：数据（先不用索引），包含处理后的点云（拟合、抽稀）
 * 对于平面点云簇，确定法向、中心点、矩形包围框，判断是否与已有平面mesh片段实际上属于同一平面，根据法向和中心点判断，
   如果是则将新的点云簇加入到已有的，然后更新已有的（重新构建），如果不是就新建mesh片段，包含构建好的mesh数据
   按照法向投影到二维，在二维构建delaunay三角剖分，并且确定mesh片段的边缘点（注意中心空洞）
 * 除此之外，判断平面点云簇是否是真平面，保真，对于真平面，执行拟合、抽稀（一段时间未更新则如此操作？？？）
 * 对于非平面点云簇，一开始保存在一个vector中，依次处理每一个点，首先根据距离阈值查找最近的非平面mesh片段，如果有的话，
   进一步判断距离边缘点的最近距离，如果在mesh投影外则xxxx，如果在mesh投影中则xxxx（同一表面/正反面），如果大于距离阈值则++++，
   如果没有的话，在空间（二维投影平面）查找未处理的最近邻点，若成功则新建mesh片段（投影重建？），若不成功则暂不处理并保留到下一帧，
   一段时间一直未处理则丢掉（视为离群点）
   无法解决的一个问题是：如何检测错误连接的face？？？？
 * 如何提取二维点云簇的边界/边缘点？
 * 未处理的点投影到当前帧的平面查找最近邻？
 * 在二维平面查找最近邻，三维距离，当前帧投影平面构建三角剖分/计算最小特征值对应的平面的三角剖分/沿激光线投影
 * 利用矩形框索查找与目标点最近的mesh片段，然后再在mesh片段中查找最近邻vertex，借鉴IDTMM的方法扩展mesh片段
 * 如何更新，先全部保留，一般法向投影在某个三角形中，一个三角形变三个，然后新点对应的菱形调整对角线避免狭长三角形
 *
 *
 *
*/

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

// 有向包围框
struct OBB {
    Eigen::Vector3d center;   // 中心点
    Eigen::Matrix3d axes;     // 轴向/变换矩阵/特征分解的特征向量按列向量排列
    Eigen::Vector3d max_pt;   // 特征向量从小到大对应xyz轴坐标系下包围盒最大坐标点
    Eigen::Vector3d min_pt;   // 最小坐标点
};

// 与坐标轴平行的包围框
struct AABB {
    Eigen::Vector3d center;   // 中心点
    Eigen::Vector3d max_pt;   // 特征向量从小到大对应xyz轴坐标系下包围盒最大坐标点
    Eigen::Vector3d min_pt;   // 最小坐标点
};


class MeshFragment
{
public:
    // 包含的点云、中心、法向量、矩形包围框
    double plane_quality;                   // 最小特征值
    PointCloudXYZI::Ptr ptcl_all;           // 所有点云-->世界坐标系下
    PointCloudXYZI::Ptr ptcl_edge;          // 抽稀后的边界点云
    PointCloudXYZI::Ptr ptcl_grid;          // 抽稀后的四叉树节点点云

    Eigen::Vector3d center_point;           // 点云簇中心点
    Eigen::Vector3d normal_vector;          // 平面法向
    AABB rec_index;                         // 有向包围盒
    double ptsNum_thres = 1000;             // 点数量超过这个阈值就降采样
    bool if_first_decimate = true;          // 是否第一次抽稀
    Eigen::Matrix3d quadtree_axes;          // 投影旋转矩阵
    Eigen::Vector3d quadtree_center;        // 四叉树根节点/投影中心点
    int plane_update_num = 0;               // 平面更新次数
    int last_grid_decimate_num = 0;             // 上次降采样平面点的数量，用于确定下次降采样的起点
    int grid_ptsNum_threshold = 4;         // 网格内点数量少于该阈值则判断该区域内没有点
//    std::vector<std::vector<std::vector<grid>>> all_grid;
//    QuadtreeTreeNode *rootNode;
    double yMax, yMin, zMax, zMin;
    double y0, z0;
    int topH, bottomH, leftW, rightW;
    std::deque<std::deque<grid>> uniform_grid;
    std::deque<std::deque<int>> tri_vertex_idx_grid;
    double dis_resolution;
    std::vector<std::array<int, 2>> grid_need_update_list;


    // 三角面片列表，包含三个点的索引
    // 点列表，包含与点链接的三角面片的索引
    std::vector<std::array<int, 3>> faces_list;           // 所有三角形列表
    std::vector<std::vector<int>> pt2face_idx;            // 每个点对应的三角形列表
    std::vector<int> faces_to_delete;                     // 最后要删除的三角形列表
    std::vector<int> pts_state;               // 点的状态 1-内部 2-边缘 3-删除

    std::vector<int> faces_to_delete_frame;                 // 删除的三角形索引
    std::vector<int> pts_edge_delete;                       // 从删除的三角形中提取的不重复的边缘点
    std::vector<int> pts_edge_add;                          // 从添加的边缘三角形中提取的不重复点，包括边缘点和内部点
    std::vector<int> pts_edge_update;                       // 2 -》2，点坐标需更新的点索引
    std::vector<int> pts_inner_add;                         // 新增的内部点，不重复
    std::vector<int> faces_with_edge_vertex;                // 新添加的包含边缘点的三角形
    std::vector<int> faces_with_2_edge_vertex;
    std::unordered_map<std::string, int> tri_delete_string_idx; // 删除的三角形索引，例：三个顶点索引 16 21 15 -》 "15 16 21"
    std::unordered_map<int, int> plane_vidx_to_noplnae_vidx;// 点在平面map中的索引 -》 点在非平面map 中的索引
    std::set<int> pro_tri_list;
    std::set<int> pro_tri_list_1;
    int plane_np = 0;
    int plane_idx = 0;

    // gui相关
    int last_faces_list_num = 0;
    int last_faces_to_delete_num = 0;
    bool if_delete = false;
    std::mutex mtx_ptcl_grid;



    MeshFragment()
    {
        // 初始化
        ptcl_all = boost::make_shared<PointCloudXYZI>();
        ptcl_edge = boost::make_shared<PointCloudXYZI>();
        ptcl_grid = boost::make_shared<PointCloudXYZI>();
//        rootNode = new QuadtreeTreeNode();
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
