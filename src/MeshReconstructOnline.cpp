//
// Created by neo on 2024/1/8.
//
// 订阅odometry.cpp发布的word坐标系下的一帧一帧点云，即连续获取的点云，增量构建mesh

/***  include  ***/
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
#include <sys/types.h>
#include <sys/stat.h>
// opencv
#include <opencv2/opencv.hpp>
#include <cv_bridge/cv_bridge.h>
// pcl
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
// ros
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <geometry_msgs/Vector3.h>
#include <image_transport/image_transport.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>
#include <tf2_msgs/TFMessage.h>
#include <livox_ros_driver/CustomMsg.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
// tbb
#include <tbb/tbb.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
// self
//#include <voxel_map/States.h>
#include <common_lib.h>
#include <so3_math.h>
#include <ikd-Tree/ikd_Tree.h>
#include "MeshFragment.h"
#include "NonPlaneMesh.h"
#include "GL_gui.h"

/***  预定义  ***/
//typedef pcl::PointXYZINormal PointType;
//typedef pcl::PointCloud<PointType> PointCloudXYZI;

#define PI_M (3.14159265358)

/***  参数定义  ***/
// 系统参数
mutex mtx_buffer_ptcl, mtx_buffer_odo;
condition_variable m_sig_buffer;
bool flg_exit = false;

// cbk参数
string ptcl_topic, odo_topic;
deque<PointCloudXYZI::Ptr> ptcl_buffer;
deque<double> time_buffer;
deque<nav_msgs::Odometry::Ptr> odo_buffer;
double last_timestamp_ptcl = -1.0;

// 点云map参数
string ptcl_save_path, dataset_name;
double frame_beg_time = 0, frame_end_time = 0;
bool if_save_raw_points;
double fov_horizon, fov_vertical_up, fov_vertical_down, angle_resolution;
KD_TREE ikdtree_map;

// 位姿参数
Eigen::Quaterniond rot_q;
Eigen::Matrix3d rot_mat;
Eigen::Vector3d pos_vec;

// 平面检测参数
int r_kNN = 3;
int plane_points_threshold;
double depth_th_normal;
double eigenvalue_threshold, angle_threshold=5, cosAngle_threshold, depth_threshold;
double angle_th_spd, dis_th_spd, cosAngle_th_spd;
double ptsNum_threshold;
mutex mtx_normal, mtx_hash, mtx_rest;

// 平面网格构建参数
double dis_resolution;
int grid_ptsNum_threshold, ptsNum_decimate_threshold;
bool if_quadtree_decimate, if_debug;

// 非平面处理参数
double voxel_resolution;
double r_np, min_pt_dis;
int N_np;

// hole_process参数
int hole_every_n_frame;

// log参数
FILE *fpTime;
FILE *fpTime_mesh;
FILE *fpTime_hole;

// gui参数
GL_gui m_gui;
vec_3f        g_axis_min_max[ 2 ];
bool if_gui;

std::unordered_map<int, std::shared_ptr<MeshFragment>> plane_idx2ptr;
std::vector<std::shared_ptr<MeshFragment>> plane_map;         // 各不相同的plane
int all_plane_num = 1;
std::shared_ptr<NonPlaneMesh> non_plane_map = std::make_shared<NonPlaneMesh>();
std::vector<int> new_plane_list;
PointCloudXYZI::Ptr pts_list_rest = boost::make_shared<PointCloudXYZI>();

/***  函数定义  ***/
void SigHandle(int sig)
{
    flg_exit = true;
    ROS_WARN("catch sig %d", sig);
    m_sig_buffer.notify_all();
}

void world_pcl_cbk(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
    mtx_buffer_ptcl.lock();
    // std::cout<<"got feature"<<endl;
    if (msg->header.stamp.toSec() < last_timestamp_ptcl) {
        ROS_ERROR("lidar pointcloud loop back, clear buffer");
        ptcl_buffer.clear();
    }
    // ROS_INFO("get point cloud at time: %.6f", msg->header.stamp.toSec());
    PointCloudXYZI::Ptr ptr(new PointCloudXYZI());
    pcl::fromROSMsg( *msg, *ptr );
    ptcl_buffer.push_back(ptr);
    time_buffer.push_back(msg->header.stamp.toSec());
    last_timestamp_ptcl = msg->header.stamp.toSec();
    mtx_buffer_ptcl.unlock();

    // condition_variable是C++中的一个线程同步原语，用于在多个线程之间进行条件变量的等待和通知。
    // notify_all()是condition_variable类的一个成员函数，用于通知所有等待在该条件变量上的线程。
    m_sig_buffer.notify_all();
}

void odometry_cbk(const nav_msgs::Odometry::ConstPtr& msg)
{
    nav_msgs::Odometry::Ptr new_msg(new nav_msgs::Odometry(*msg));
    mtx_buffer_odo.lock();

    odo_buffer.push_back(new_msg);

    mtx_buffer_odo.unlock();
    m_sig_buffer.notify_all();
}

// 判断两个平面是否实际上是同一个平面，
// 先判断法向是否相似，然后判断中心点连线的向量沿法向投影的大小（平行的平面），然后（栅栏式平面）
bool if_same_plane( const std::shared_ptr<MeshFragment>& plane1, const std::shared_ptr<MeshFragment>& plane2,
                    double thre_cos0, double thre_dis)
{
    AABB index1 = plane1->rec_index;
    AABB index2 = plane2->rec_index;;
    Eigen::Vector3d norm_vec1 = plane1->normal_vector;
    Eigen::Vector3d norm_vec2 = plane2->normal_vector;

    // AABB3D相交检测
    for (int i = 0; i < 3; i++)
    {
        if (index1.min_pt(i) > index2.max_pt(i) || index1.max_pt(i) < index2.min_pt(i))
        {
            return false;
        }
    }

    double dotProduct = norm_vec1.dot(norm_vec2);
    // 法向相似性，-1->1，越1越相似
    if ( dotProduct < thre_cos0 )
    {
        return false;
    }

    Eigen::Vector3d center_pt1 = index1.center;
    Eigen::Vector3d center_pt2 = index2.center;
    Eigen::Vector3d norm_avg = (norm_vec1 + norm_vec2).normalized();
    Eigen::Vector3d p2q = center_pt1 - center_pt2;
    double norm_dis = norm_avg.dot(p2q);
    // 中心点连线与法向点积，-1->1，越0越相似
    if ( std::abs(norm_dis) > thre_dis )
    {
        return false;
    }

    return true;
}

// plane1是否与plane2是同一平面
bool if_same_plane_2( const std::shared_ptr<MeshFragment>& plane1, const std::shared_ptr<MeshFragment>& plane2 )
{
    Eigen::Vector3d norm_vec1 = plane1->normal_vector;
    Eigen::Vector3d norm_vec2 = plane2->normal_vector;
    AABB index1 = plane1->rec_index;
    AABB index2 = plane2->rec_index;

    double dotProduct = norm_vec1.dot(norm_vec2);
    // 法向相似性，-1->1，越1越相似
    // TODO 阈值0-5度 是否合适？
    double thre_cos0 = 0.984808;   // cos5=0.996195 cos10=0.984808
    if ( dotProduct < thre_cos0 )
    {
        return false;
    }

    if (plane2->if_first_decimate)
    {
        // 未进行过降采样
        Eigen::Vector3d center_pt1 = index1.center;
        Eigen::Vector3d center_pt2 = index2.center;
        Eigen::Vector3d norm_avg = (norm_vec1 + norm_vec2).normalized();
        Eigen::Vector3d p2q = center_pt1 - center_pt2;
        double norm_dis = norm_avg.dot(p2q);
        // 中心点连线与法向点积，-1->1，越0越相似
        // TODO 阈值0.05m 是否合适？
        double thre_dis = 0.03;
        if ( std::abs(norm_dis) > thre_dis )
        {
            return false;
        }
    }
    else
    {
        Eigen::Vector3d center_pt1 = index1.center;
        Eigen::Vector3d pt1_pro_plane2 = plane2->quadtree_axes.transpose() * (center_pt1 - plane2->quadtree_center);
        double y_start = plane2->y0 - dis_resolution * plane2->leftW;
        double z_start = plane2->z0 + dis_resolution * plane2->topH;
        int u = int((z_start - pt1_pro_plane2(2)) / dis_resolution);
        int v = int((pt1_pro_plane2(1) - y_start) / dis_resolution);
        int H = plane2->bottomH + plane2->topH;
        int W = plane2->leftW + plane2->rightW;
        if (u > H - 1) { u = H - 1; }
        if (u < 0) { u = 0; }
        if (v > W - 1) { v = W - 1; }
        if (v < 0) { u = 0; }
        double x_plane1 = pt1_pro_plane2(0);
        double x_plane2 = -10.0;
        while (x_plane2 == -10.0)
        {
            int range = std::max(H, W);
            for (int i = 0; i < range; i++)
            {
                for (int j = -i; j <= i; j++)
                {
                    int temp_u = u - i;
                    int temp_v = v + j;
                    if (temp_u < 0 || temp_u > H-1 || temp_v < 0 || temp_v > W-1)
                    {
                        continue;
                    }
                    if (plane2->uniform_grid[temp_u][temp_v].node_type != 0)
                    {
                        x_plane2 = plane2->uniform_grid[temp_u][temp_v].x_avg;
                        break;
                    }
                }
                if (x_plane2 != -10.0) { break; }
                for (int j = -i; j <= i; j++)
                {
                    int temp_u = u + i;
                    int temp_v = v + j;
                    if (temp_u < 0 || temp_u > H-1 || temp_v < 0 || temp_v > W-1)
                    {
                        continue;
                    }
                    if (plane2->uniform_grid[temp_u][temp_v].node_type != 0)
                    {
                        x_plane2 = plane2->uniform_grid[temp_u][temp_v].x_avg;
                        break;
                    }
                }
                if (x_plane2 != -10.0) { break; }
                for (int j = -i; j <= i; j++)
                {
                    int temp_u = u + j;
                    int temp_v = v - i;
                    if (temp_u < 0 || temp_u > H-1 || temp_v < 0 || temp_v > W-1)
                    {
                        continue;
                    }
                    if (plane2->uniform_grid[temp_u][temp_v].node_type != 0)
                    {
                        x_plane2 = plane2->uniform_grid[temp_u][temp_v].x_avg;
                        break;
                    }
                }
                if (x_plane2 != -10.0) { break; }
                for (int j = -i; j <= i; j++)
                {
                    int temp_u = u + j;
                    int temp_v = v + i;
                    if (temp_u < 0 || temp_u > H-1 || temp_v < 0 || temp_v > W-1)
                    {
                        continue;
                    }
                    if (plane2->uniform_grid[temp_u][temp_v].node_type != 0)
                    {
                        x_plane2 = plane2->uniform_grid[temp_u][temp_v].x_avg;
                        break;
                    }
                }
                if (x_plane2 != -10.0) { break; }

                if ( i == range-1 )
                {
                    x_plane2 = 0;
                }
            }
        }
        double thre_dis = 0.03;
        if ( std::abs(std::abs(x_plane1)-std::abs(x_plane2)) > thre_dis )
        {
            return false;
        }
    }

    // AABB3D相交检测
    for (int i = 0; i < 3; i++)
    {
        if (index1.min_pt(i) > index2.max_pt(i) || index1.max_pt(i) < index2.min_pt(i))
        {
            return false;
            break;
        }
    }
    return true;
}

// 判断某个世界坐标系下的点是否投影在某个平面的内部
bool if_project_inner_plane(const Eigen::Vector3d& pt, const std::shared_ptr<MeshFragment>& plane)
{
    AABB index = plane->rec_index;
    // AABB包围框相交检测
    for (int i = 0; i < 3; i++)
    {
        if (pt(i) > index.max_pt(i) || pt(i) < index.min_pt(i))
        {
            // 点不在平面的矩形框内
            return false;
        }
    }

    // 投影检测
    if (plane->if_first_decimate)
    {
        // 平面若还没降采样，相近非平面点不处理
        return false;
    }
    else
    {
        Eigen::Vector3d pt_pro_plane = plane->quadtree_axes.transpose() * (pt - plane->quadtree_center);
        double y_start = plane->y0 - dis_resolution * plane->leftW;
        double z_start = plane->z0 + dis_resolution * plane->topH;
        int u = int((z_start - pt_pro_plane(2)) / dis_resolution);
        int v = int((pt_pro_plane(1) - y_start) / dis_resolution);
        int H = plane->bottomH + plane->topH;
        int W = plane->leftW + plane->rightW;
        if (u < 0 || u > H - 1 || v < 0 || v > W - 1)
        {
            // 投影在平面外
            return false;
        }
        if (plane->uniform_grid[u][v].node_type != 1)
        {
            // 投影网格还不是内部点
            return false;
        }
        if (std::abs(pt_pro_plane(0) - plane->uniform_grid[u][v].x_avg) > 0.2)
        {
            // 虽投影在平面内，但与平面法向距离太大
            return false;
        }
    }
    return true;
}

// 保存点云到.bin文件
void saveCloudToBin(const pcl::PointCloud<pcl::PointXYZI>::Ptr &cloud, const std::string &filename) {
    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile.is_open()) {
        std::cout << "Couldn't open file to write\n";
        return;
    }

    for (const auto &point : cloud->points) {
        outFile.write(reinterpret_cast<const char *>(&point.x), sizeof(float));
        outFile.write(reinterpret_cast<const char *>(&point.y), sizeof(float));
        outFile.write(reinterpret_cast<const char *>(&point.z), sizeof(float));
        outFile.write(reinterpret_cast<const char *>(&point.intensity), sizeof(float));
    }

    outFile.close();
}

/***  主函数  ***/
int main(int argc, char **argv) {
    ros::init(argc, argv, "MeshReconstructOnline");
    ros::NodeHandle m_node;

    // common params
    m_node.param<string>("mesh/ptcl_topic", ptcl_topic, "/cloud_registered");
    m_node.param<string>("mesh/odo_topic", odo_topic, "/Odometry");

    // 点云map参数
    m_node.param<string>("mesh/ptcl_save_path", ptcl_save_path, "");
    m_node.param<string>("mesh/dataset_name", dataset_name, "unknown");
    m_node.param<bool>("mesh/if_save_raw_points", if_save_raw_points, false);
    m_node.param<double>("mesh/angle_resolution", angle_resolution, 1);
    m_node.param<double>("mesh/fov_horizon", fov_horizon, 80);
    m_node.param<double>("mesh/fov_vertical_up", fov_vertical_up, 80);
    m_node.param<double>("mesh/fov_vertical_down", fov_vertical_down, 80);

    // 平面点提取参数
    m_node.param<int>("mesh/r_kNN", r_kNN, 3);                            // 计算法向量时邻域大小，2*r_kNN+1
    m_node.param<double>("mesh/depth_th_normal", depth_th_normal, 0.5);
    m_node.param<double >("mesh/eigenvalue_threshold", eigenvalue_threshold, 0.01);  // 最小特征值阈值 -> 判断平面点
    m_node.param<int>("mesh/plane_points_threshold", plane_points_threshold, 20);    // 某个平面包含的点的数量超过阈值才计数--未用
    m_node.param<double>("mesh/angle_threshold", angle_threshold, 5);                // 法向相似性阈值
    m_node.param<double>("mesh/depth_threshold", depth_threshold, 0.05);             // 深度相似性阈值
    m_node.param<double>("mesh/angle_th_spd", angle_th_spd, 5);
    m_node.param<double>("mesh/dis_th_spd", dis_th_spd, 0.05);
    m_node.param<double>("mesh/ptsNum_threshold", ptsNum_threshold, 50);             // 平面点数量阈值
    m_node.param<double>("mesh/dis_resolution", dis_resolution, 0.1);                // 平面划分网格的边长
    m_node.param<int>("mesh/grid_ptsNum_threshold", grid_ptsNum_threshold, 4);
    m_node.param<int>("mesh/ptsNum_decimate_threshold", ptsNum_decimate_threshold, 1000);
    m_node.param<bool>("mesh/if_quadtree_decimate", if_quadtree_decimate, false);    // true时，整个平面维持一个x值
    m_node.param<bool>("mesh/if_gui", if_gui, false);
    m_node.param<bool>("mesh/if_debug", if_debug, false);

    // 非平面点处理参数
    m_node.param<double>("mesh/voxel_resolution", voxel_resolution, 0.03);           // 体素大小
    m_node.param<double>("mesh/r_nonplane", r_np, 0.1);                              // 查找最近点的距离阈值
    m_node.param<int>("mesh/N_nonplane", N_np, 5);                                   // 拟合两个三角形的点数阈值
    m_node.param<double>("mesh/points_minimum_scale", min_pt_dis, 0.01);             // 任意两个三角形顶点的最小距离阈值

    // hole_process参数
    m_node.param<int>("mesh/hole_every_n_frame", hole_every_n_frame, 10);

    // 创建以系统当前时间命名的文件夹
    mkdir(ptcl_save_path.c_str(), 0777);

    auto now = std::chrono::system_clock::now();
    std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm now_tm = *std::localtime(&now_time_t);
    std::stringstream ss;
    ss << std::put_time(&now_tm, "%Y-%m-%d_%H-%M-%S");
    std::string file_name_str = ss.str();
    ptcl_save_path += "/" + dataset_name + "_" + file_name_str + "/";
    mkdir(ptcl_save_path.c_str(), 0777);


    // 读取所有参数并保存
    FILE *fp_ros_param;
    std::string all_ros_param_save_path = ptcl_save_path + "all_ros_param.txt";
    fp_ros_param = fopen(all_ros_param_save_path.c_str(), "w");
    // 获取当前系统时间点
    std::stringstream ss_;
    ss_ << std::put_time(&now_tm, "%Y-%m-%d %H:%M:%S");
    std::string time_str = ss_.str();
    fprintf(fp_ros_param, "Current Time: %s\n", time_str.c_str());

    ros::V_string ros_keys;
    m_node.getParamNames(ros_keys);
    for (const auto& param_name : ros_keys) {
        // 尝试读取不同类型的参数
        double double_value;
        bool bool_value;
        std::string param_value;

        // 读取参数，根据实际情况选择读取方式
        if (ros::param::get(param_name, double_value)) {
            fprintf(fp_ros_param, "%s: %f \n", param_name.c_str(), double_value);
        } else if (ros::param::get(param_name, bool_value)) {
            fprintf(fp_ros_param, "%s: %s \n", param_name.c_str(), bool_value ? "true" : "false");
        } else {
            // 如果不是上述类型，尝试读取为字符串
            if (ros::param::get(param_name, param_value)) {
                fprintf(fp_ros_param, "%s: %s \n", param_name.c_str(), param_value.c_str());
            } else {
                ROS_WARN("Parameter type for '%s' not supported or not found.", param_name.c_str());
            }
        }
    }
    fclose(fp_ros_param);
    std::cout << "All parameters have been saved to " << all_ros_param_save_path << std::endl;


    // 订阅器
    ros::Subscriber sub_ptcl = m_node.subscribe(ptcl_topic, 200000, world_pcl_cbk);
    ros::Subscriber sub_odo = m_node.subscribe(odo_topic, 100000, odometry_cbk);

    // 参数
    PointCloudXYZI::Ptr pcl_wait_save(new PointCloudXYZI());
    int frame_idx = 0;
    cosAngle_threshold = std::cos(angle_threshold / 180 * 3.1415926);
    cosAngle_th_spd = std::cos(angle_th_spd / 180 * 3.1415926);
    bool is_first_frame = true;
    int H = int((fov_vertical_up + fov_vertical_down) / angle_resolution) + 2 * r_kNN + 12;
    int H_up = int(fov_vertical_up / angle_resolution) + r_kNN + 6;
    int W = int(fov_horizon  / angle_resolution) + 2 * r_kNN + 12;
    double rad2angle = 180.0 / 3.1415926;

    non_plane_map->voxel_resolution = voxel_resolution;
    non_plane_map->r_dis = r_np;
    non_plane_map->N_pts = N_np;
    non_plane_map->minimum_pt_dis = min_pt_dis;

    string log_save_path3 = ptcl_save_path + "log_time.txt";
    fpTime = fopen(log_save_path3.c_str(), "w");

    // 打开GUI界面，放在if判断中会报错
    std::thread thr_gui = std::thread( &GL_gui::display, &m_gui );
    PointCloudXYZI::Ptr pcl_normal(new PointCloudXYZI());

//------------------------------------------------------------------------------------------------------
    signal(SIGINT, SigHandle);      // 收到退出(ctrl-C)信号时执行 SigHandle()函数
    ros::Rate rate(5000);           // 通过rate.sleep()函数控制节点按照指定的帧率运行
    bool status = ros::ok();        // 用于检查ROS节点是否仍在运行。它返回一个布尔值，如果节点仍在运行，则返回true，否则返回false。

    while (status) {
        if (flg_exit)
            break;
        ros::spinOnce();                // 读取ROS待处理队列消息

        if ( ptcl_buffer.empty() || time_buffer.empty() || odo_buffer.empty() ) {
            continue;
        }

        auto total_start = std::chrono::high_resolution_clock::now();

        cout << "***Mesh Reconstruction*** points buffer size is " << ptcl_buffer.size() << endl;
        PointCloudXYZI::Ptr ptcl_frame = ptcl_buffer.front();
        frame_beg_time = time_buffer.front();
        frame_end_time = frame_beg_time + ptcl_frame->points.back().curvature / double(1000);
        time_buffer.pop_front();
        ptcl_buffer.pop_front();

        nav_msgs::Odometry::Ptr odo_frame = odo_buffer.front();
        odo_buffer.pop_front();
        rot_q.x() = odo_frame->pose.pose.orientation.x;
        rot_q.y() = odo_frame->pose.pose.orientation.y;
        rot_q.z() = odo_frame->pose.pose.orientation.z;
        rot_q.w() = odo_frame->pose.pose.orientation.w;
        rot_mat = rot_q.toRotationMatrix();
        pos_vec(0) = odo_frame->pose.pose.position.x;
        pos_vec(1) = odo_frame->pose.pose.position.y;
        pos_vec(2) = odo_frame->pose.pose.position.z;

        // gui处理
        Eigen::Matrix< double, 7, 1 > pose_vec;
        pose_vec.head< 4 >() = rot_q.coeffs();
        pose_vec.block( 4, 0, 3, 1 ) = pos_vec;
        // 每帧的点云xyzi + pose
        for ( int i = 0; i < ptcl_frame->points.size(); i++ )
        {
            m_gui.g_eigen_vec_vec[ frame_idx ].first.emplace_back( ptcl_frame->points[ i ].x, ptcl_frame->points[ i ].y, ptcl_frame->points[ i ].z,
                                                             ptcl_frame->points[ i ].intensity );
        }
        m_gui.g_eigen_vec_vec[ frame_idx ].second = pose_vec;
        m_gui.g_current_frame = frame_idx;


        double time_diff = std::abs( frame_beg_time - odo_frame->header.stamp.toSec());
        if (time_diff >= 0.1)
        {
            std::cout << "\033[31m Warning: time_diff is  \033[0m" << time_diff << std::endl;
        }

        // 保存点云数据, pcd格式
        if ( if_save_raw_points ) {
            *pcl_wait_save += *ptcl_frame;
        }

        fprintf(fpTime, "frame: %d ", frame_idx);
        std::cout << "***Mesh Reconstruction*** frame_idx: " << frame_idx << "  pts_size: " << ptcl_frame->points.size() << std::endl;

        /****** 点云投影 ******/
        int ptcl_size = ptcl_frame->points.size();
        int u, v;
        Tensor<double, 3> projectImage(H, W, 11);  // x y z xx xy xz yy yz zz N intensity
        projectImage.setConstant(0.0f);

        // 世界坐标系 -> 体坐标系
        Eigen::MatrixXd frame_world_pt_eigen(ptcl_size, 3);
        Eigen::MatrixXd frame_body_pt_eigen(ptcl_size, 3);
        for (int i = 0; i < ptcl_size; ++i)
        {
            frame_world_pt_eigen(i, 0) = ptcl_frame->points[i].x;
            frame_world_pt_eigen(i, 1) = ptcl_frame->points[i].y;
            frame_world_pt_eigen(i, 2) = ptcl_frame->points[i].z;
        }
        frame_body_pt_eigen = (frame_world_pt_eigen.rowwise() - pos_vec.transpose()) * rot_mat;

        for ( int i = 0; i < ptcl_size; i++ )
        {
            // world -> body
            double temp_x, temp_y, temp_z;
            temp_x = frame_body_pt_eigen(i, 0);
            temp_y = frame_body_pt_eigen(i, 1);
            temp_z = frame_body_pt_eigen(i, 2);

            // 根据lidar的FOV，先柱面后展开
            double temp_x2 = temp_x * temp_x;
            double temp_y2 = temp_y * temp_y;
            double temp_z2 = temp_z * temp_z;
            double point_depth = std::sqrt(temp_x2 + temp_y2 + temp_z2);
            double xy_sqrt = std::sqrt(temp_x2 + temp_y2);
            double point_theti = std::acos(temp_x / xy_sqrt) * rad2angle;
            double point_alphi = std::atan(temp_z / xy_sqrt) * rad2angle;
            double angle_sign;
            if (temp_y >= 0) { angle_sign = -1.0; }
            else { angle_sign = 1.0; }
            u = int(angle_sign * point_theti / angle_resolution) + W / 2;
            v = int(-1.0 * point_alphi / angle_resolution) + H_up;

            if ( u<0 || u>=W || v<0 || v>=H)
            {
                cout << "\033[31m project error! \033[0m" << H << " " << W << " " << u << " " << v << endl;   // 输出红色字体
                continue;
            }

            if ( projectImage(v, u, 3) == 0 )
            {
                projectImage(v, u, 0) = temp_x;
                projectImage(v, u, 1) = temp_y;
                projectImage(v, u, 2) = temp_z;
                projectImage(v, u, 3) = point_depth;
                projectImage(v, u, 10) = ptcl_frame->points[ i ].intensity;
                projectImage(v, u, 9) = 1;
            }
            else
            {
                if ( projectImage(v, u, 3) > point_depth )
                {
                    projectImage(v, u, 0) = temp_x;
                    projectImage(v, u, 1) = temp_y;
                    projectImage(v, u, 2) = temp_z;
                    projectImage(v, u, 3) = point_depth;
                    projectImage(v, u, 10) = ptcl_frame->points[ i ].intensity;
                    projectImage(v, u, 9) = 1;
                }
            }
        }

        auto end_1 = std::chrono::high_resolution_clock::now();
        auto time_1 = std::chrono::duration_cast<std::chrono::duration<double>>(end_1 - total_start).count() * 1000;
        fprintf(fpTime, "Project: %.4f ms ", time_1);

        /****** 一般方法计算法向，并行 ******/
        Tensor<double, 3> normalImage(H, W, 4);
        normalImage.setConstant(-1.0f);
        int parallel_num = (H - r_kNN * 2 - 2) * (W - r_kNN * 2 - 2);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, parallel_num), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t ind = r.begin(); ind != r.end(); ++ind)
            {
                int i = ind / (W - r_kNN * 2 - 2) + r_kNN + 1;
                int j = ind % (W - r_kNN * 2 - 2) + r_kNN + 1;

                // 当前位置没有点投影
                if ( projectImage(i, j, 9) == 0 ) { continue; }

                double cur_x = projectImage(i, j, 3);
                std::vector<Eigen::Vector2i> temp_Pt_uv;
                for (int m = i-r_kNN; m <= i+r_kNN; m++)
                {
                    for (int n = j-r_kNN; n <= j+r_kNN; n++)
                    {
                        if (projectImage(m, n, 9) == 0) { continue; }
                        // 深度差太大说明属于同一个平面的概率比较小
                        if (std::abs(projectImage(m, n, 3) - cur_x) > depth_th_normal) { continue; }

                        Eigen::Vector2i temp_pt_uv = {m, n};
                        temp_Pt_uv.push_back(temp_pt_uv);
                    }
                }

                if (temp_Pt_uv.size() < 10/*r_kNN * 2*/) { continue; }  // 最近邻点的数量小于阈值则不计算法向

                Eigen::MatrixXd eigen_cloud(temp_Pt_uv.size(), 3);
                int temp_i = 0;
                for (auto it : temp_Pt_uv)
                {
                    eigen_cloud(temp_i, 0) = projectImage(it(0), it(1), 0);
                    eigen_cloud(temp_i, 1) = projectImage(it(0), it(1), 1);
                    eigen_cloud(temp_i, 2) = projectImage(it(0), it(1), 2);
                    temp_i++;
                }
                Eigen::Vector3d center_point = eigen_cloud.colwise().mean();
                Eigen::Matrix3d point_cov = (eigen_cloud.rowwise() - center_point.transpose()).transpose() *
                                            (eigen_cloud.rowwise() - center_point.transpose()) / eigen_cloud.rows();

                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver( point_cov );
                Eigen::Vector3d eigenvalues = eigensolver.eigenvalues();
                Eigen::Matrix3d eigenvectors = eigensolver.eigenvectors();

                // 统一平面法向指向外侧
                Eigen::Vector3d lidar_point_beam;
                lidar_point_beam(0) = projectImage(i, j, 0);
                lidar_point_beam(1) = projectImage(i, j, 1);
                lidar_point_beam(2) = projectImage(i, j, 2);
                Eigen::Vector3d norm_vec = eigenvectors.col(0);
                mtx_normal.lock();
                if ( lidar_point_beam.dot(norm_vec) > 0)
                {
                    normalImage(i, j, 1) = -eigenvectors(0, 0);     // 列向量
                    normalImage(i, j, 2) = -eigenvectors(1, 0);
                    normalImage(i, j, 3) = -eigenvectors(2, 0);
                }
                else
                {
                    normalImage(i, j, 1) = eigenvectors(0, 0);     // 列向量
                    normalImage(i, j, 2) = eigenvectors(1, 0);
                    normalImage(i, j, 3) = eigenvectors(2, 0);
                }
                normalImage(i, j, 0) = eigenvalues(0);         // 从小到大排列
                mtx_normal.unlock();
            }
        });


        auto end_3 = std::chrono::high_resolution_clock::now();
        auto time_3 = std::chrono::duration_cast<std::chrono::duration<double>>(end_3 - end_1).count() * 1000;
        fprintf(fpTime, "Normal: %.4f ms ", time_3);

        /****** 区域生长/广度优先->平面分割/分类 ******/
        std::vector<std::vector<Eigen::Vector2i>> plane_points_uv;
        Tensor<double, 2> seg_result(H, W);
        seg_result = projectImage.chip(9, 2);
        std::queue<Eigen::Vector2i> points_q;
        double plane_num = 2;

        for ( int i = r_kNN+1; i < H-r_kNN-1; i++ )
        {
            for (int j = r_kNN+1; j < W-r_kNN-1; j++)
            {
                if ( normalImage(i, j, 0) < 0 || normalImage(i, j, 0) > eigenvalue_threshold
                    || seg_result(i, j) > 1 )  // || abs(normalImage(k, l, 1)) < 0.1    // TODO 考虑不同深度平面接缝？
                {
                    continue;
                }

                std::vector<Eigen::Vector2i> temp_planePt_uv;
                Eigen::Vector3d norm_ij, seg_center;
                norm_ij(0) = normalImage(i, j, 1);
                norm_ij(1) = normalImage(i, j, 2);
                norm_ij(2) = normalImage(i, j, 3);
                seg_center(0) = projectImage(i, j, 0);
                seg_center(1) = projectImage(i, j, 1);
                seg_center(2) = projectImage(i, j, 2);

                for (int m = i-r_kNN; m <= i+r_kNN; m++)
                {
                    for (int n = j-r_kNN; n <= j+r_kNN; n++)
                    {
                        if (seg_result(m, n) == 0 || seg_result(m, n) > 1) { continue; }  // TODO ？

                        // 距离法向投影相似性
                        Eigen::Vector3d p2q;
                        p2q(0) = seg_center(0) - projectImage(m, n, 0);
                        p2q(1) = seg_center(1) - projectImage(m, n, 1);
                        p2q(2) = seg_center(2) - projectImage(m, n, 2);
                        if (abs(p2q.dot(norm_ij)) > depth_threshold) { continue; }
                        // 法向相似性
                        Eigen::Vector3d norm_mn;
                        norm_mn(0) = normalImage(m, n, 1);
                        norm_mn(1) = normalImage(m, n, 2);
                        norm_mn(2) = normalImage(m, n, 3);
                        if (norm_ij.dot(norm_mn) < cosAngle_threshold) { continue; }

                        seg_result(m, n) = plane_num;

                        Eigen::Vector2i temp_pt_uv = {m, n};
                        temp_planePt_uv.push_back(temp_pt_uv);

                        if (m != i || n != j) {
                            Eigen::Vector2i point_uv = {m, n};
                            points_q.push(point_uv);     // 候选平面点
                        }
                    }
                }

                while (!points_q.empty())
                {
                    Eigen::Vector2i currentPoint = points_q.front();
                    points_q.pop();
                    int k = currentPoint(0);
                    int l = currentPoint(1);
                    for (int m = k-r_kNN; m <= k+r_kNN; m++)
                    {
                        for (int n = l-r_kNN; n <= l+r_kNN; n++)
                        {
                            if (seg_result(m, n) == 0 || seg_result(m, n) > 1) { continue; }

                            Eigen::Vector3d p2q;
                            p2q(0) = seg_center(0) - projectImage(m, n, 0);
                            p2q(1) = seg_center(1) - projectImage(m, n, 1);
                            p2q(2) = seg_center(2) - projectImage(m, n, 2);
                            if (abs(p2q.dot(norm_ij)) > depth_threshold) { continue; }

                            Eigen::Vector3d norm_mn;
                            norm_mn(0) = normalImage(m, n, 1);
                            norm_mn(1) = normalImage(m, n, 2);
                            norm_mn(2) = normalImage(m, n, 3);
                            if (norm_ij.dot(norm_mn) < cosAngle_threshold) { continue; }

                            seg_result(m, n) = plane_num;

                            Eigen::Vector2i temp_pt_uv = {m, n};
                            temp_planePt_uv.push_back(temp_pt_uv);

                            if (m != k || n != l) {
                                Eigen::Vector2i point_uv = {m, n};
                                points_q.push(point_uv);
                            }
                        }
                    }
                }
                plane_points_uv.push_back(temp_planePt_uv);
                plane_num++;
            }
        }

        auto end_4 = std::chrono::high_resolution_clock::now();
        auto time_4 = std::chrono::duration_cast<std::chrono::duration<double>>(end_4 - end_3).count() * 1000;
        fprintf(fpTime, "Segment: %.4f ms ", time_4);

        /****** 构建平面mesh ******/
        /* 第一帧直接新建，添加平面点、计算平面参数，添加到平面vector中
         * 之后帧，先检测相似平面，同一帧不判断是否同一平面，只与之前帧的所有平面比较
         * 若没有则新建，若有则将所有属于同一平面的mesh片段融合为一个
        */

        new_plane_list.clear();
        for ( int i = 0; i < plane_points_uv.size(); i++ )
        {
            // 点数太少不处理
            if ( plane_points_uv[i].size() < ptsNum_threshold )
            {
                // 改为非平面点
                for (auto it : plane_points_uv[i])
                {
                    seg_result(it(0), it(1)) = 1;
                }
                continue;
            }

            std::shared_ptr<MeshFragment> temp_mf = std::make_shared<MeshFragment>();
            temp_mf->ptsNum_thres = ptsNum_decimate_threshold;
            temp_mf->grid_ptsNum_threshold = grid_ptsNum_threshold;
            temp_mf->give_plane_point(projectImage, plane_points_uv[i], rot_mat, pos_vec);
            temp_mf->compute_plane_parameter(pos_vec);

            // 判断平面质量   TODO 阈值
            if (temp_mf->plane_quality > eigenvalue_threshold * 10)
            {
                for (auto it : plane_points_uv[i])
                {
                    seg_result(it(0), it(1)) = 1;
                }
                continue;
            }

            if (frame_idx == 0)
            {
                // 降采样
                temp_mf->grid_decimate(dis_resolution, fpTime_mesh);
                // 传递信息给np
                if (!temp_mf->if_first_decimate)
                {
                    non_plane_map->mesh_update_from_plane_map(temp_mf);
                }

                temp_mf->plane_idx = all_plane_num;
                plane_idx2ptr[all_plane_num] = temp_mf;         // 为了后面查找该平面指针
                new_plane_list.push_back(all_plane_num);        // 用于np_p_connect
                all_plane_num++;
                plane_map.push_back(temp_mf);
            }
            else
            {
                int indexJ = 0;
                std::vector<int> same_plane_list;
                for (auto plane : plane_map)
                {
                    if ( if_same_plane(temp_mf, plane, cosAngle_th_spd, dis_th_spd) )
//                    if ( if_same_plane_2(temp_mf, plane) )
                    {
                        same_plane_list.push_back(indexJ);
                    }
                    indexJ++;
                }

                if (same_plane_list.size() > 0)
                {
                    // 更新平面点、平面参数
                    plane_map[same_plane_list[0]]->update_plane(temp_mf);
                    new_plane_list.push_back(plane_map[same_plane_list[0]]->plane_idx);
                    for ( int j = 1; j < same_plane_list.size(); j++ )
                    {
                        plane_map[same_plane_list[0]]->update_plane(plane_map[same_plane_list[j]]);
                        new_plane_list.push_back(plane_map[same_plane_list[j]]->plane_idx);
                    }

                    // 降采样新加入的点
                    plane_map[same_plane_list[0]]->grid_decimate(dis_resolution, fpTime_mesh);
                    // 传递信息给np
                    if (!plane_map[same_plane_list[0]]->if_first_decimate)
                    {
                        // TODO 删除的非平面能否分给平面？
                        non_plane_map->mesh_update_from_plane_map(plane_map[same_plane_list[0]]);
                    }
                    for ( int j = 1; j < same_plane_list.size(); j++ )
                    {
                        non_plane_map->mesh_delete_plane_map(plane_map[same_plane_list[j]]);
                    }

                    // 删除plane_map中除第一个相同平面mesh片段外的所有其他相同平面mesh片段
                    if (same_plane_list.size() > 1)
                    {
                        // 删除已经合并的平面
                        for ( int j = 1; j < same_plane_list.size(); j++ )
                        {
                            plane_map[same_plane_list[j]]->if_delete = true;
                        }

                        for ( int j = 1; j < same_plane_list.size(); j++ )
                        {
                            std::vector<std::shared_ptr<MeshFragment>>::size_type st = same_plane_list[j]-(j-1);
                            plane_map.erase(plane_map.begin() + st);
                        }
                    }
                }
                else
                {
                    // 降采样
                    temp_mf->grid_decimate(dis_resolution, fpTime_mesh);
                    // 传递信息给np
                    if (!temp_mf->if_first_decimate)
                    {
                        non_plane_map->mesh_update_from_plane_map(temp_mf);
                    }

                    temp_mf->plane_idx = all_plane_num;
                    plane_idx2ptr[all_plane_num] = temp_mf;
                    new_plane_list.push_back(all_plane_num);
                    all_plane_num++;
                    plane_map.push_back(temp_mf);
                }
            }
        }

        auto end_5 = std::chrono::high_resolution_clock::now();
        auto time_5 = std::chrono::duration_cast<std::chrono::duration<double>>(end_5 - end_4).count() * 1000;
        fprintf(fpTime, "planePoints mesh: %.4f ms ", time_5);

        std::cout << "***Mesh Reconstruction*** plane num is " << plane_map.size() << std::endl;

        /****** 非平面点构建mesh ******/
        // 建立体素的hash索引
        std::vector<std::array<long, 3>> voxel_to_update;
        int parallel_num_1 = H * W;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, parallel_num_1), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t ind = r.begin(); ind != r.end(); ++ind)
            {
                int i = ind / W;
                int j = ind % W;

                // 平面外的点
                if (seg_result(i, j) != 1 /*|| normalImage(i, j, 0) == -1*/) { continue; }

                Eigen::Vector3d p_body;
                p_body(0) = projectImage(i, j, 0);
                p_body(1) = projectImage(i, j, 1);
                p_body(2) = projectImage(i, j, 2);
                Eigen::Vector3d p_global(rot_mat * (p_body) + pos_vec);

                // TODO 跳过的非平面点可否分给平面？，遍历平面耗时长
                bool if_pro_in_plane = false;
                for (auto plane : plane_map)
                {
                    if ( non_plane_map->if_project_inner_plane(p_global, plane) )
                    {
                        if_pro_in_plane = true;
                        break;
                    }
                }

                pcl::PointXYZINormal temp_pt;
                temp_pt.x = p_global(0);
                temp_pt.y = p_global(1);
                temp_pt.z = p_global(2);
                temp_pt.normal_x = normalImage(i, j, 1);
                temp_pt.normal_y = normalImage(i, j, 2);
                temp_pt.normal_z = normalImage(i, j, 3);

                if (if_pro_in_plane) {
                    mtx_rest.lock();
                    pts_list_rest->points.push_back(temp_pt);
                    mtx_rest.unlock();
                    continue; }

                mtx_hash.lock();
                non_plane_map->pts_list->push_back(temp_pt);
                non_plane_map->pts_state.push_back(0);
                non_plane_map->pts_process_num.push_back(0);
                non_plane_map->hole_process_num.push_back(0);
                non_plane_map->hole_v_process.push_back(0);
                non_plane_map->pt2face_idx.resize(non_plane_map->pts_list->size());

                int current_idx = non_plane_map->pts_list->size() - 1;
                long voxel_x = std::round(p_global(0) / voxel_resolution);
                long voxel_y = std::round(p_global(1) / voxel_resolution);
                long voxel_z = std::round(p_global(2) / voxel_resolution);
                if (non_plane_map->voxel_pts.if_exist(voxel_x, voxel_y, voxel_z))
                {
                    non_plane_map->voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->all_pts_idx.push_back(current_idx);
                }
                else
                {
                    std::shared_ptr<m_voxel> temp_voxel = std::make_shared<m_voxel>();
                    non_plane_map->voxel_pts.insert(voxel_x, voxel_y, voxel_z, temp_voxel);
                    non_plane_map->voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->all_pts_idx.push_back(current_idx);
                }

                if (!non_plane_map->voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->if_need_update)
                {
                    std::array<long, 3> temp_voxel_idx;
                    temp_voxel_idx[0] = voxel_x;
                    temp_voxel_idx[1] = voxel_y;
                    temp_voxel_idx[2] = voxel_z;
                    voxel_to_update.push_back(temp_voxel_idx);
                    non_plane_map->voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->if_need_update = true;
                }
                mtx_hash.unlock();
            }
        });

        auto end_6 = std::chrono::high_resolution_clock::now();
        auto time_6 = std::chrono::duration_cast<std::chrono::duration<double>>(end_6 - end_5).count() * 1000;
        fprintf(fpTime, "np_hash_1: %.4f ms ", time_6);

        // 体素按空间位置分块，间隔区域大小保证各分块之间不会互相影响
        std::vector<std::vector<std::array<long, 3>>> voxel_to_update_parallel;
        std::vector<std::vector<std::array<long, 3>>> edge_voxel_to_update_1;      // 按照x轴-》y轴-》z轴的顺序
        std::vector<std::vector<std::array<long, 3>>> edge_voxel_to_update_2;
        std::vector<std::vector<std::array<long, 3>>> edge_voxel_to_update_3;
        Hash_map_3d_base< long, int > center_voxel_xyz_to_index;
        Hash_map_3d_base< long, int > edge_voxel_x_to_index;
        Hash_map_3d_base< long, int > edge_voxel_xy_to_index;
        Hash_map_3d_base< long, int > edge_voxel_xyz_to_index;
        int parallel_voxel_num = 3;

        for (auto it : voxel_to_update)
        {
            long p_voxel[3];
            long mod_voxel[3];
            for (int i = 0; i < 3; i++)
            {
                if (it[i] >= 0)
                {
                    p_voxel[i] = it[i] / parallel_voxel_num;
                    mod_voxel[i] = it[i] % parallel_voxel_num;
                }
                else
                {
                    p_voxel[i] = (it[i] + 1) / parallel_voxel_num - 1;
                    mod_voxel[i] = it[i] % parallel_voxel_num + parallel_voxel_num;
                }
            }

            // 间隔区域
            if (mod_voxel[0] == parallel_voxel_num - 1 || mod_voxel[0] == parallel_voxel_num - 2)
            {
                if (edge_voxel_x_to_index.if_exist(p_voxel[0], 0, 0))
                {
                    int temp_ind = edge_voxel_x_to_index.m_map_3d_hash_map[p_voxel[0]][0][0];
                    std::array<long, 3> temp_voxel_idx;
                    temp_voxel_idx[0] = it[0];
                    temp_voxel_idx[1] = it[1];
                    temp_voxel_idx[2] = it[2];
                    edge_voxel_to_update_1[temp_ind].push_back(temp_voxel_idx);
                }
                else
                {
                    int cur_ind = edge_voxel_to_update_1.size();
                    edge_voxel_x_to_index.insert(p_voxel[0], 0, 0, cur_ind);
                    std::array<long, 3> temp_voxel_idx;
                    temp_voxel_idx[0] = it[0];
                    temp_voxel_idx[1] = it[1];
                    temp_voxel_idx[2] = it[2];
                    std::vector<std::array<long, 3>> temp_voxel_region;
                    temp_voxel_region.push_back(temp_voxel_idx);
                    edge_voxel_to_update_1.push_back(temp_voxel_region);
                }
            }
            else if (mod_voxel[1] == parallel_voxel_num - 1 || mod_voxel[1] == parallel_voxel_num - 2)
            {
                if (edge_voxel_xy_to_index.if_exist(p_voxel[0], p_voxel[1], 0))
                {
                    int temp_ind = edge_voxel_xy_to_index.m_map_3d_hash_map[p_voxel[0]][p_voxel[1]][0];
                    std::array<long, 3> temp_voxel_idx;
                    temp_voxel_idx[0] = it[0];
                    temp_voxel_idx[1] = it[1];
                    temp_voxel_idx[2] = it[2];
                    edge_voxel_to_update_2[temp_ind].push_back(temp_voxel_idx);
                }
                else
                {
                    int cur_ind = edge_voxel_to_update_2.size();
                    edge_voxel_xy_to_index.insert(p_voxel[0], p_voxel[1], 0, cur_ind);
                    std::array<long, 3> temp_voxel_idx;
                    temp_voxel_idx[0] = it[0];
                    temp_voxel_idx[1] = it[1];
                    temp_voxel_idx[2] = it[2];
                    std::vector<std::array<long, 3>> temp_voxel_region;
                    temp_voxel_region.push_back(temp_voxel_idx);
                    edge_voxel_to_update_2.push_back(temp_voxel_region);
                }
            }
            else if (mod_voxel[2] == parallel_voxel_num - 1 || mod_voxel[2] == parallel_voxel_num - 2)
            {
                if (edge_voxel_xyz_to_index.if_exist(p_voxel[0], p_voxel[1], p_voxel[2]))
                {
                    int temp_ind = edge_voxel_xyz_to_index.m_map_3d_hash_map[p_voxel[0]][p_voxel[1]][p_voxel[2]];
                    std::array<long, 3> temp_voxel_idx;
                    temp_voxel_idx[0] = it[0];
                    temp_voxel_idx[1] = it[1];
                    temp_voxel_idx[2] = it[2];
                    edge_voxel_to_update_3[temp_ind].push_back(temp_voxel_idx);
                }
                else
                {
                    int cur_ind = edge_voxel_to_update_3.size();
                    edge_voxel_xyz_to_index.insert(p_voxel[0], p_voxel[1], p_voxel[2], cur_ind);
                    std::array<long, 3> temp_voxel_idx;
                    temp_voxel_idx[0] = it[0];
                    temp_voxel_idx[1] = it[1];
                    temp_voxel_idx[2] = it[2];
                    std::vector<std::array<long, 3>> temp_voxel_region;
                    temp_voxel_region.push_back(temp_voxel_idx);
                    edge_voxel_to_update_3.push_back(temp_voxel_region);
                }
            }
            // 内部区域
            else
            {
                if (center_voxel_xyz_to_index.if_exist(p_voxel[0], p_voxel[1], p_voxel[2]))
                {
                    int temp_ind = center_voxel_xyz_to_index.m_map_3d_hash_map[p_voxel[0]][p_voxel[1]][p_voxel[2]];
                    std::array<long, 3> temp_voxel_idx;
                    temp_voxel_idx[0] = it[0];
                    temp_voxel_idx[1] = it[1];
                    temp_voxel_idx[2] = it[2];
                    voxel_to_update_parallel[temp_ind].push_back(temp_voxel_idx);
                }
                else
                {
                    int cur_ind = voxel_to_update_parallel.size();
                    center_voxel_xyz_to_index.insert(p_voxel[0], p_voxel[1], p_voxel[2], cur_ind);
                    std::array<long, 3> temp_voxel_idx;
                    temp_voxel_idx[0] = it[0];
                    temp_voxel_idx[1] = it[1];
                    temp_voxel_idx[2] = it[2];
                    std::vector<std::array<long, 3>> temp_voxel_region;
                    temp_voxel_region.push_back(temp_voxel_idx);
                    voxel_to_update_parallel.push_back(temp_voxel_region);
                }
            }
        }

        auto end_66 = std::chrono::high_resolution_clock::now();
        auto time_66 = std::chrono::duration_cast<std::chrono::duration<double>>(end_66 - end_6).count() * 1000;
        fprintf(fpTime, "np_hash_2: %.4f ms ", time_66);

        // 逐体素、逐点处理，孤立体素新建mesh，体素内的点先处理法向投影在三角形内部的，再处理投影在外部
        // 对于投影在外部的点，根据边缘点查找最近然后扩展该点
//        non_plane_map->mesh_update(voxel_to_update, pos_vec);

        // 并行处理体素块
        tbb::parallel_for_each( voxel_to_update_parallel.begin(), voxel_to_update_parallel.end(),
                                [ & ]( const std::vector<std::array<long, 3>> &voxel_region ) {

            non_plane_map->mesh_update(voxel_region, pos_vec);
        } );

        auto end_10 = std::chrono::high_resolution_clock::now();
        auto time_10 = std::chrono::duration_cast<std::chrono::duration<double>>(end_10 - end_66).count() * 1000;
        fprintf(fpTime, "np_mesh_p: %.4f ms ", time_10);

        // 最后处理间隔区域
        tbb::parallel_for_each( edge_voxel_to_update_1.begin(), edge_voxel_to_update_1.end(),
                                [ & ]( const std::vector<std::array<long, 3>> &voxel_region ) {

            non_plane_map->mesh_update(voxel_region, pos_vec);
        } );

        tbb::parallel_for_each( edge_voxel_to_update_2.begin(), edge_voxel_to_update_2.end(),
                                [ & ]( const std::vector<std::array<long, 3>> &voxel_region ) {

            non_plane_map->mesh_update(voxel_region, pos_vec);
        } );

        tbb::parallel_for_each( edge_voxel_to_update_3.begin(), edge_voxel_to_update_3.end(),
                                [ & ]( const std::vector<std::array<long, 3>> &voxel_region ) {

            non_plane_map->mesh_update(voxel_region, pos_vec);
        } );

        auto end_7 = std::chrono::high_resolution_clock::now();
        auto time_7 = std::chrono::duration_cast<std::chrono::duration<double>>(end_7 - end_10).count() * 1000;
        fprintf(fpTime, "np_mesh_l: %.4f ms ", time_7);

        // 处理空洞
        if (frame_idx % hole_every_n_frame == 0)
        {
            non_plane_map->mesh_hole_process(fpTime_hole);
        }

        auto end_8 = std::chrono::high_resolution_clock::now();
        auto time_8 = std::chrono::duration_cast<std::chrono::duration<double>>(end_8 - end_7).count() * 1000;
        fprintf(fpTime, "np_hole: %.4f ms ", time_8);

        // 连接平面非平面
        if (frame_idx % hole_every_n_frame == 0)
        {
            non_plane_map->plane_and_noplane_connect();
        }

        auto end_9 = std::chrono::high_resolution_clock::now();
        auto time_9 = std::chrono::duration_cast<std::chrono::duration<double>>(end_9 - end_8).count() * 1000;
        fprintf(fpTime, "npp_connect: %.4f ms ", time_9);

        // g_triangles_manager更新
        if (if_gui)
        {
            m_gui.triangle_manager_update(new_plane_list);
            m_gui.service_refresh_and_synchronize_triangle();
        }

        auto end_12 = std::chrono::high_resolution_clock::now();
        auto time_12 = std::chrono::duration_cast<std::chrono::duration<double>>(end_12 - end_9).count() * 1000;
        fprintf(fpTime, "gui: %.4f ms ", time_12);

        frame_idx++;

        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::duration<double>>(total_end - total_start).count() * 1000;
        cout << "***Mesh Reconstruction*** Total cost time is " << total_time << " ms!" << endl;
        cout << "---------------------------------------------------------------" << endl;
        fprintf(fpTime, "Total: %.4f ms \n", total_time);

        status = ros::ok();
        rate.sleep();
    }

    non_plane_map->mesh_hole_process(fpTime_hole);
    non_plane_map->plane_and_noplane_connect();

    /**************** 保存带法向的点云 ****************/
//    string normal_path = ptcl_save_path + "00_normal_mylio.ply";
//    pcl::PolygonMesh mesh;
//    mesh.header = pcl::PCLHeader();
//    mesh.cloud = pcl::PCLPointCloud2();
//    mesh.polygons = std::vector<pcl::Vertices>();
//    pcl::toPCLPointCloud2(*pcl_normal, mesh.cloud);
//    pcl::io::savePLYFileBinary(normal_path, mesh);

    /**************** save subscribe map ****************/
    /* 1. make sure you have enough memories
    /* 2. pcd save will largely influence the real-time performences **/
    if ( if_save_raw_points )
    {
        string file_name = string("points.pcd");
        string all_points_dir(ptcl_save_path + file_name);
        pcl::PCDWriter pcd_writer;
        std::cout << "Points saved to " << all_points_dir << endl << endl;
        pcd_writer.writeBinary(all_points_dir, *pcl_wait_save);
    }

    /****** 保存合并后的平面mesh ******/
    string save_dir = ptcl_save_path + "p.ply";
    string mp_save_dir = ptcl_save_path + "m_p.ply";
    std::shared_ptr<MeshFragment> merged_plane = std::make_shared<MeshFragment>();;
    for (int i = 0; i < plane_map.size(); i++)
    {
        auto mf = plane_map[i];
        if ( mf->ptcl_grid->points.size() > 0)
        {
            merged_plane->vertex_and_face_list_merge(mf);
        }
    }

    if (if_debug)
    {
        // 去除不连接三角形的点
        merged_plane->save_to_ply_without_redundancy(save_dir);
    }

    /****** 四叉树降采样 quadtree decimate ******/
    std::shared_ptr<MeshFragment> merged_plane_1 = std::make_shared<MeshFragment>();
    if (if_quadtree_decimate)
    {
        // 并行
//        tbb::parallel_for_each( plane_map.begin(), plane_map.end(),
//                                [ & ]( const std::shared_ptr<MeshFragment> &plane ) {
//                                    plane->quadtree_decimate();
//                                } );
        // 串行
        for (auto it : plane_map)
        {
            it->quadtree_decimate();
        }

        // 保存合并起来的平面mesh
        for (int i = 0; i < plane_map.size(); i++)
        {
            auto mf = plane_map[i];
            if ( mf->ptcl_grid->points.size() > 0)
            {
                merged_plane_1->vertex_and_face_list_merge(mf);
            }
        }
        if (if_debug)
        {
            string save_dir_pQd = ptcl_save_path + "p_Qdecimate.ply";
            // 去除不连接三角形的点
            merged_plane_1->save_to_ply_without_redundancy(save_dir_pQd);
        }
    }

    /****** 保存非平面mesh、自定义mesh、点云 ******/
    string mesh_save_dir = ptcl_save_path + "np.ply";
    if (if_debug)
    {
        // 去除不连接三角形的点
        non_plane_map->save_to_ply_without_redundancy(mesh_save_dir);
    }

    if (if_quadtree_decimate)
    {
        for (auto it : plane_map)
        {
            non_plane_map->plane_vertex_adjust(it);
        }

        if (if_debug)
        {
            // 去除不连接三角形的点
            string mesh_save_dir_pQd = ptcl_save_path + "np_Qdecimate.ply";
            non_plane_map->save_to_ply_without_redundancy(mesh_save_dir_pQd);
        }
    }

    /****** 将平面mesh、非平面mesh合并保存 ******/
    if (if_quadtree_decimate)
    {
        int last_pts_num = non_plane_map->pts_list->points.size();
        int last_tri_num = non_plane_map->faces_list_all.size();
        *non_plane_map->pts_list += *(merged_plane_1->ptcl_grid);

        for (int i = 0; i < merged_plane_1->faces_list.size(); i++)
        {
            std::array<int, 3> new_tri_idx;
            new_tri_idx[0] = merged_plane_1->faces_list[i][0] + last_pts_num;
            new_tri_idx[1] = merged_plane_1->faces_list[i][1] + last_pts_num;
            new_tri_idx[2] = merged_plane_1->faces_list[i][2] + last_pts_num;
            non_plane_map->faces_list_all.push_back(new_tri_idx);
        }

        for (int i = 0; i < merged_plane_1->faces_to_delete.size(); i++)
        {
            int delete_tri_idx = merged_plane_1->faces_to_delete[i];
            non_plane_map->faces_to_delete.push_back(delete_tri_idx + last_tri_num);
        }
    }
    else
    {
        int last_pts_num = non_plane_map->pts_list->points.size();
        int last_tri_num = non_plane_map->faces_list_all.size();
        *non_plane_map->pts_list += *(merged_plane->ptcl_grid);

        for (int i = 0; i < merged_plane->faces_list.size(); i++)
        {
            std::array<int, 3> new_tri_idx;
            new_tri_idx[0] = merged_plane->faces_list[i][0] + last_pts_num;
            new_tri_idx[1] = merged_plane->faces_list[i][1] + last_pts_num;
            new_tri_idx[2] = merged_plane->faces_list[i][2] + last_pts_num;
            non_plane_map->faces_list_all.push_back(new_tri_idx);
        }

        for (int i = 0; i < merged_plane->faces_to_delete.size(); i++)
        {
            int delete_tri_idx = merged_plane->faces_to_delete[i];
            non_plane_map->faces_to_delete.push_back(delete_tri_idx + last_tri_num);
        }
    }

    // 去除不连接三角形的点
    string mesh_save_dir_all = ptcl_save_path + "all.ply";
    non_plane_map->save_to_ply_without_redundancy(mesh_save_dir_all);


    std::cout << "over!" << std::endl;

    fclose(fpTime);

    return 0;
}

































