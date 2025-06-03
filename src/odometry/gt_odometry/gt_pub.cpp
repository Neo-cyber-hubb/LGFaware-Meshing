//
// Created by neo on 2024/9/20.
/*
 * 读取点云（每帧以单独文件保存），真值pose（4*4矩阵的前三行）
 * 然后publish
 * */

#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <pcl/register_point_struct.h>
#include <pcl/filters/voxel_grid.h>
//#include <pcl/io/txt_io.h>
#include <fstream>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Dense>
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
#include <thread>



// 读取txt格式的pose
void Mtrix2Pose(const Eigen::Matrix<double, 12, 1>& M, Eigen::Matrix<double, 3, 3>& R, Eigen::Vector3d& t)
{
    R(0, 0) = M(0);
    R(0, 1) = M(1);
    R(0, 2) = M(2);
    R(1, 0) = M(4);
    R(1, 1) = M(5);
    R(1, 2) = M(6);
    R(2, 0) = M(8);
    R(2, 1) = M(9);
    R(2, 2) = M(10);

    t(0) = M(3);
    t(1) = M(7);
    t(2) = M(11);
}

// 读取ply格式的点云
bool readPointCloudFromFile(const std::string& filename, pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud) {
    if (pcl::io::loadPLYFile<pcl::PointXYZ>(filename, *cloud) == -1) {
        ROS_ERROR("Failed to read PLY file: %s", filename.c_str());
        return false;
    }
    return true;
}

// 读取bin格式的点云
bool read_bin(const std::string & file_path, pcl::PointCloud<pcl::PointXYZ> & laser_cloud)
{
    std::cout<<std::setprecision(7)<<setiosflags(std::ios::fixed);
    std::ifstream bin_point_cloud_file(file_path, std::ifstream::in | std::ifstream::binary);
    if(!bin_point_cloud_file.good()){
        return false;
    }
    bin_point_cloud_file.seekg(0, std::ios::end);
    const size_t num_elements = bin_point_cloud_file.tellg() / sizeof(float);
    bin_point_cloud_file.seekg(0, std::ios::beg);
    std::vector<float> lidar_data(num_elements);
    bin_point_cloud_file.read(reinterpret_cast<char*>(&lidar_data[0]), num_elements * sizeof(float));
//    std::cout << "totally " << int(lidar_data.size() / 4.0) << " points in this lidar frame " << line_num << "\n";

    std::vector<Eigen::Vector3d> lidar_points;
    std::vector<float> lidar_intensities;
    for (std::size_t i = 0; i < lidar_data.size(); i += 4)
    {
        lidar_points.emplace_back(lidar_data[i], lidar_data[i+1], lidar_data[i+2]);
        lidar_intensities.push_back(lidar_data[i+3]);

        pcl::PointXYZ point;
        point.x = lidar_data[i];
        point.y = lidar_data[i + 1];
        point.z = lidar_data[i + 2];
        //if(point.z > -2.5){//there are some underground outliers in kitti dataset, remove them
        laser_cloud.push_back(point);
        //}
    }
    return true;
}

bool read_pcd(const std::string & file_path, pcl::PointCloud<pcl::PointXYZ> & laser_cloud)
{
    if (pcl::io::loadPCDFile<pcl::PointXYZ>(file_path, laser_cloud) == -1)
    {
        return false;
    }
    return true;
}

// 点云坐标系转换
void pointBodyToWorld(pcl::PointXYZ const *const pi, pcl::PointXYZ *const po,
                      const Eigen::Matrix3d& R, const Eigen::Vector3d& t)
{
    Eigen::Vector3d p_body(pi->x, pi->y, pi->z);
    Eigen::Vector3d p_global(R * p_body + t);
    po->x = p_global(0);
    po->y = p_global(1);
    po->z = p_global(2);
}

int main(int argc, char** argv) {
    ros::init(argc, argv, "pose_and_pointcloud_publisher");
    ros::NodeHandle nh;

    // 定义参数
    std::string data_dir, pose_dir, calib_dir, ply_name_dir;
    std::string pose_topic, ptcl_topic;
    int dataset_num;
    int frame_max_points_num;
    int pub_n_frame = -1;
    nav_msgs::Path path;
    nav_msgs::Odometry odomAftMapped;
    geometry_msgs::Quaternion geoQuat;
//    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_all(new pcl::PointCloud<pcl::PointXYZ>);

    // 读取参数
    nh.param<std::string>("odometry/data_dir", data_dir, "/xxxx");
    nh.param<std::string>("odometry/lid_topic", ptcl_topic, "/cloud_registered");
    nh.param<std::string>("odometry/pose_topic", pose_topic, "/aft_mapped_to_init");
    nh.param<int>("odometry/dataset_num", dataset_num, 2);
    nh.param<int>("odometry/frame_max_points_num", frame_max_points_num, 15000);
    nh.param<int>("odometry/pub_n_frame", pub_n_frame, -1);

    ros::Publisher pubGlobalPoints = nh.advertise<sensor_msgs::PointCloud2>(ptcl_topic, 100000);
    ros::Publisher pubPose = nh.advertise<nav_msgs::Odometry>(pose_topic, 100000);

    // 读取全部 pose
    pose_dir = data_dir + "/poses.txt";
    calib_dir = data_dir + "/calib.txt";
    std::ifstream file_pose(pose_dir);
    if (!file_pose.is_open()) {
        ROS_ERROR("Failed to open pose file: %s", pose_dir.c_str());
    }
    std::string line;
    std::vector<Eigen::Matrix<double, 12, 1>> pose_list;
    while (std::getline(file_pose, line)) {
        std::stringstream ss(line);
        std::vector<double> values;
        double value;

        while (ss >> value) {
            values.push_back(value);
        }

        Eigen::Matrix<double, 12, 1> pose_frame;
        for (int i = 0; i < values.size(); ++i) {
            pose_frame(i) = values[i];
        }
        pose_list.push_back(pose_frame);
    }
    file_pose.close();

    // 读取外参
    Eigen::Matrix4d calib_T;
    if (dataset_num == 1 || dataset_num == 2 || dataset_num == 3)
    {
        std::ifstream file_calib(calib_dir);
        std::string line_tr;
        while (std::getline(file_calib, line_tr)) {
            if (line_tr.find("Tr:") != std::string::npos) {
                std::istringstream iss(line_tr.substr(line_tr.find("Tr:") + 3));
                std::vector<double> values;
                double value;
                while (iss >> value) {
                    values.push_back(value);
                }

                calib_T(0, 0) = values[0];
                calib_T(0, 1) = values[1];
                calib_T(0, 2) = values[2];
                calib_T(0, 3) = values[3];
                calib_T(1, 0) = values[4];
                calib_T(1, 1) = values[5];
                calib_T(1, 2) = values[6];
                calib_T(1, 3) = values[7];
                calib_T(2, 0) = values[8];
                calib_T(2, 1) = values[9];
                calib_T(2, 2) = values[10];
                calib_T(2, 3) = values[11];
                calib_T(3, 0) = 0;
                calib_T(3, 1) = 0;
                calib_T(3, 2) = 0;
                calib_T(3, 3) = 1;

                break;
            }
        }
        std::cout << "***GT odometry*** extrinsic matrix:\n" << calib_T << std::endl;
    }

    // 等待3s 以保证mesh_module初始化完毕
    std::this_thread::sleep_for( std::chrono::milliseconds( 3000 ) );

    // 创建一个循环频率为10Hz的Rate对象
    ros::Rate loop_rate(10);
    int frame_idx = 0;

    while (ros::ok())
    {
        // 只 pub 前 pub_n_frame 帧
        if (pub_n_frame != -1 && frame_idx >= pub_n_frame)
        {
            break;
        }

        // 读取一帧lidar数据
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        std::stringstream point_cloud_path;

        if(dataset_num == 1) //kitti
        {
            point_cloud_path << data_dir << "/velodyne/" << std::setfill('0') << std::setw(6) << frame_idx << ".bin";
            if (!read_bin(point_cloud_path.str(), *cloud)) {
                std::cout << "***GT odometry*** No more data or path error!" << std::endl;
                break;
            }
        }
        else if(dataset_num == 2) //mai_city
        {
            point_cloud_path << data_dir << "/velodyne/" << std::setfill('0') << std::setw(5) << frame_idx << ".bin";
            if (!read_bin(point_cloud_path.str(), *cloud)) {
                std::cout << "***GT odometry*** No more data or path error!" << std::endl;
                break;
            }
        }
        else if(dataset_num == 3) //ncd
        {
            int ncd_frame_idx = frame_idx + 500;
            point_cloud_path << data_dir << "/pcd/" << std::setfill('0') << std::setw(5) << ncd_frame_idx << ".pcd";
            if (!read_pcd(point_cloud_path.str(), *cloud)) {
                std::cout << "***GT odometry*** No more data or path error!" << std::endl;
                break;
            }
        }
        else if(dataset_num == 5) //StreetRecon
        {
            point_cloud_path << data_dir << "/lidars/" << std::setfill('0') << std::setw(8) << frame_idx << ".pcd";
            if (!read_pcd(point_cloud_path.str(), *cloud)) {
                std::cout << "***GT odometry*** No more data or path error!" << std::endl;
                break;
            }
        }
        else
        {
            std::cout << "***GT odometry*** dataset_num false!" << std::endl;
            break;
        }
//        std::cout << "***GT odometry*** Load " << cloud->points.size() << " points from" << point_cloud_path.str() << std::endl;

        // 获取对应的pose
        Eigen::Vector3d t;
        Eigen::Matrix<double, 3, 3> rotation_matrix;
        Mtrix2Pose(pose_list[frame_idx], rotation_matrix, t);

        if (dataset_num == 1 || dataset_num == 2 || dataset_num == 3)
        {
            Eigen::Matrix4d pose_T;
            pose_T.block<3, 3>(0, 0) = rotation_matrix;
            pose_T.block<3, 1>(0, 3) = t;
            pose_T(3, 0) = 0;
            pose_T(3, 1) = 0;
            pose_T(3, 2) = 0;
            pose_T(3, 3) = 1;

            Eigen::Matrix4d pose_T_ = pose_T * calib_T;
            rotation_matrix = pose_T_.block<3, 3>(0, 0);
            t = pose_T_.block<3, 1>(0, 3);
        }

        // 点云转到世界坐标系下
        int size = cloud->points.size();
//        int point_step = int(size / frame_max_points_num);

//        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_world( new pcl::PointCloud<pcl::PointXYZ>(size, 1));
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_world( new pcl::PointCloud<pcl::PointXYZ>);
        for (int i = 0; i < size; i++)
        {
//            pointBodyToWorld(&cloud->points[i], &cloud_world->points[i], rotation_matrix, t);
            if (dataset_num == 3)
            {
                if (std::abs(cloud->points[i].x) < 0.5 ||
                    std::abs(cloud->points[i].y) < 0.5 ||
                    std::abs(cloud->points[i].z) < 0.5 )
                {
                    continue;
                }
            }
            pcl::PointXYZ temp_pt;
            pointBodyToWorld(&cloud->points[i], &temp_pt, rotation_matrix, t);
            cloud_world->points.push_back(temp_pt);
        }

        // 点云体素降采样
        if (cloud->points.size() > frame_max_points_num && 0)
        {
            pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);
            // 创建VoxelGrid滤波器
            pcl::VoxelGrid<pcl::PointXYZ> vg;
            vg.setInputCloud(cloud_world);
            vg.setLeafSize(0.01f, 0.01f, 0.01f); // 设置体素大小为1cm
            vg.filter(*cloud_filtered);
            cloud_world = cloud_filtered;
            // 输出降采样前后的点云信息
            std::cout << "***GT odometry*** PointCloud's number after voxel downsample: " << cloud_world->points.size() << std::endl;
        }

//        *cloud_all += *cloud_world;

        // 发布点云
        auto current_frame_header_stamp = ros::Time::now();
        sensor_msgs::PointCloud2 laserCloudmsg;
        pcl::toROSMsg(*cloud_world, laserCloudmsg);
        laserCloudmsg.header.stamp = current_frame_header_stamp;
        laserCloudmsg.header.frame_id = "world";
        pubGlobalPoints.publish(laserCloudmsg);

        // 发布pose
//        Eigen::Vector3d euler_cur = RotMtoEuler( rotation_matrix );
//        geoQuat = tf::createQuaternionMsgFromRollPitchYaw( euler_cur( 0 ), euler_cur( 1 ), euler_cur( 2 ) );
        Eigen::Quaterniond q = Eigen::Quaterniond(rotation_matrix);
        q.normalize();
        odomAftMapped.header.frame_id = "world";
        odomAftMapped.child_frame_id = "body";
        odomAftMapped.header.stamp = current_frame_header_stamp;
        odomAftMapped.pose.pose.orientation.x = q.x();
        odomAftMapped.pose.pose.orientation.y = q.y();
        odomAftMapped.pose.pose.orientation.z = q.z();
        odomAftMapped.pose.pose.orientation.w = q.w();
        odomAftMapped.pose.pose.position.x = t( 0 );
        odomAftMapped.pose.pose.position.y = t( 1 );
        odomAftMapped.pose.pose.position.z = t( 2 );

        pubPose.publish( odomAftMapped );

        // 发布父坐标系到子坐标系的变换关系
        static tf::TransformBroadcaster br;
        tf::Transform                   transform;
        tf::Quaternion                  q_;
        transform.setOrigin(
                tf::Vector3( odomAftMapped.pose.pose.position.x, odomAftMapped.pose.pose.position.y, odomAftMapped.pose.pose.position.z ) );
        q_.setW( odomAftMapped.pose.pose.orientation.w );
        q_.setX( odomAftMapped.pose.pose.orientation.x );
        q_.setY( odomAftMapped.pose.pose.orientation.y );
        q_.setZ( odomAftMapped.pose.pose.orientation.z );
        transform.setRotation( q_ );
        br.sendTransform( tf::StampedTransform( transform, current_frame_header_stamp, "world", "body" ) );


        // 控制发布频率10Hz
        frame_idx++;
//        if (frame_idx > 200)
//        {
//            break;
//        }

        ros::spinOnce();
        loop_rate.sleep();
    }

//    std::string all_points_dir = "/home/neo/data_ncd/points.pcd";
//    pcl::PCDWriter pcd_writer;
//    std::cout << "Points saved to " << all_points_dir << std::endl;
//    pcd_writer.writeBinary(all_points_dir, *cloud_all);

    return 0;
}