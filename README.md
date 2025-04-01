## 1. Introduction
Paper

LGFaware-Meshing: Online Mesh Reconstruction from LiDAR Point Cloud with Awareness of Local Geometric Features


## 2. Install
### 2.1 **Ubuntu** and **ROS**
Ubuntu = 20.04

ROS    = Noetic  鱼香ROS一键安装
```
wget http://fishros.com/install -O fishros && . fishros
```
Its additional package
```
sudo apt-get install ros-noetic-cv-bridge ros-noetic-tf ros-noetic-message-filters ros-noetic-image-transport*
```

### 2.2. **Dependency**
PCL    >= 1.8,   Follow [PCL Installation](http://www.pointclouds.org/downloads/linux.html).

Eigen  >= 3.3.4, Follow [Eigen Installation](http://eigen.tuxfamily.org/index.php?title=Main_Page).

Optional: Use ```cmake -D CMAKE_INSTALL_PREFIX=/usr/local ..``` to define customized installation path
```
cd eigen-3.3.9
mkdir build && cd build
cmake ..
sudo make install
```

CGAL and OpenGL
```
sudo apt-get install -y libcgal-dev pcl-tools
sudo apt-get install -y libgl-dev libglm-dev libglfw3-dev libglew-dev libglw1-mesa-dev 
sudo apt-get install -y libcgal-dev libxkbcommon-x11-dev
```


### 2.3. **livox_ros_driver**
Follow [livox_ros_driver Installation](https://github.com/Livox-SDK/livox_ros_driver).

### 2.4. **Build OnlineMeshing on ROS**
Note: disable anaconda3, for example, rename.
```
mkdir ws_onlinemesh && cd ws_onlinemesh
mkdir src && cd src
git clone https://github.com/Neo-cyber-hubb/onlinemesh_private.git
cd ../
source $Livox_ros_driver_dir$/devel/setup.bash
catkin_make 
source devel/setup.bash
```

## 3. Run
Firstly, run ```source devel/setup.bash```

### 3.1. **Mai City dataset**
[Mai City dataset download]()

Run with odometry (voxelmap), set saved path of results ``ptcl_save_path`` in ```&project_path&/config/velodyne_maicity.yaml```. Then open a terminal under path ```ws_onlinemesh``` and run:
```
roslaunch online_mesh mapping_maicity.launch
```
Open another terminal under the path of Mai City dataset and run:
```
rosbag play 01.bag
```
Run with ground truth pose, set dataset path ```data_dir``` and saved path of results ``ptcl_save_path`` in ```&project_path&/config/maicity_gt_mesh.yaml```. Then open a terminal under path ```ws_onlinemesh``` and run:
```
roslaunch online_mesh mapping_maicity_gtpose.launch
```

### 3.2. **KITTI dataset**
[KITTI dataset download]()

Run with odometry (voxelmap), set saved path of results ``ptcl_save_path`` in ```&project_path&/config/velodyne_kitti.yaml```. Then open a terminal under path ```ws_onlinemesh``` and run:
```
roslaunch online_mesh mapping_kitti.launch
```
Convert ```.bin``` to bag file. Open another terminal under the path of KITTI dataset and run:
```
rosbag play 07.bag
```
Run with ground truth pose, set dataset path ```data_dir``` and saved path of results ``ptcl_save_path`` in ```&project_path&/config/kitti_gt_mesh.yaml```. Then open a terminal under path ```ws_onlinemesh``` and run:
```
roslaunch online_mesh mapping_kitti_gtpose.launch
```

### 3.2. **R3live dataset**
[R3live dataset download]()

Run with odometry (voxelmap or fastlio), set saved path of results ``ptcl_save_path`` in ```&project_path&/config/avia_r3live.yaml```. Then open a terminal under path ```ws_onlinemesh``` and run:
```
roslaunch online_mesh mapping_r3live.launch
```
Open another terminal under the path of Mai City dataset and run:
```
rosbag play hku_campus_seq_00.bag
```

## 4. Evaluation
Evaluate with [mesh-evaluation](https://github.com/Neo-cyber-hubb/mesh-evaluation)

## 5. Note
### 5.1. **Parameters**
| Parameters | Explanation | Note |
|----------|------|------|
| hole_every_n_frame | Carry out hole filling and connection every ```n``` frame|  |
| dis_resolution | Size of plane grid |  |
| grid_ptsNum_threshold | Threshold to define occupancy of plane grid|  |
| voxel_resolution | Size of non-plane voxel |  |
| points_minimum_scale | The minimum distance between any two non-plane points  |  |
| r_nonplane | Nearest neighbor search distance threshold in non-plane process |  |
| N_nonplane | Threshold to fit triangle of isolated points|  |
| if_save_raw_points | Save subscribed points or not |  |
| if_quadtree_decimate | Carry out quadtree decimate or not |  |
| if_gui | Display real-time reconstruction or not | ```true```: can't real-time |
| if_debug | Save np.ply and p.ply or not|  |


## 6. Acknowledgments
Thanks for [ImMesh](https://github.com/hku-mars/ImMesh), 
[FastLIO](https://github.com/hku-mars/FAST_LIO),
[VoxelMap](https://github.com/hku-mars/VoxelMap).