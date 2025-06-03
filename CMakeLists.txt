cmake_minimum_required(VERSION 2.8.3)
project(online_mesh)

SET(CMAKE_BUILD_TYPE "Release")

ADD_COMPILE_OPTIONS(-std=c++14 )

add_definitions(-DROOT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fexceptions" )
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++14 -g -O3 -ltbb -lboost_system -msse2 -msse3 -pthread -w")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -pthread -std=c++0x -std=c++14 -fexceptions")

message("Current CPU archtecture: ${CMAKE_SYSTEM_PROCESSOR}")
if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)" )
  include(ProcessorCount)
  ProcessorCount(N)
  message("Processer number:  ${N}")
  if(N GREATER 5)
    add_definitions(-DMP_EN)
    add_definitions(-DMP_PROC_NUM=4)
    message("core for MP:  4")
  elseif(N GREATER 3)
    math(EXPR PROC_NUM "${N} - 2")
    add_definitions(-DMP_EN)
    add_definitions(-DMP_PROC_NUM="${PROC_NUM}")
    message("core for MP:  ${PROC_NUM}")
  else()
    add_definitions(-DMP_PROC_NUM=1)
  endif()
else()
  add_definitions(-DMP_PROC_NUM=1)
endif()

find_package(OpenMP QUIET)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")

find_package(catkin REQUIRED COMPONENTS
        geometry_msgs
        nav_msgs
        sensor_msgs
        roscpp
        cv_bridge
        rospy
        std_msgs
        pcl_ros
        tf
        livox_ros_driver
        message_generation
        eigen_conversions
        cv_bridge
        )

find_package(Eigen3 REQUIRED)
find_package(PCL 1.8 REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Ceres REQUIRED)
#find_package(OpenCV)
find_package(CGAL REQUIRED)
FIND_PACKAGE(Boost REQUIRED COMPONENTS filesystem iostreams program_options system serialization)
set(Sophus_LIBRARIES libSophus.so)

message(Eigen: ${EIGEN3_INCLUDE_DIR})

include_directories(
        ${catkin_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        ${PCL_INCLUDE_DIRS}
        ${PYTHON_INCLUDE_DIRS}
        ${Sophus_INCLUDE_DIRS})

add_message_files(
        FILES
        Pose6D.msg
        States.msg
)

generate_messages(
        DEPENDENCIES
        geometry_msgs
)

catkin_package(
        CATKIN_DEPENDS geometry_msgs nav_msgs roscpp rospy std_msgs message_runtime
        DEPENDS EIGEN3 PCL OpenCV Sophus
        INCLUDE_DIRS
)

set(FONTPATH ${CMAKE_CURRENT_SOURCE_DIR}/font/)
add_definitions(-DFONT_DIR=\"${FONTPATH}/\")

# 定义了一个预处理器宏 SHADER_DIR，其值为变量 SHADERPATH 的值，后跟一个斜杠。在C++代码中，你可以通过预定义的宏来访问这个路径
# 在CMake脚本中使用这些命令可以方便地在多个目标（如可执行文件或库）之间共享着色器路径，同时也可以方便地修改路径而不需要更改多个文件。
# 这种设置通常用在C++项目中，特别是那些需要加载外部资源（如着色器文件）的项目。通过在编译时定义这些路径，程序可以在运行时知道去哪里查找所需的资源。
set(SHADERPATH ${CMAKE_CURRENT_SOURCE_DIR}/src/shader/)
add_definitions(-DSHADER_DIR=\"${SHADERPATH}/\")

add_executable(my_mesh
        ./src/MeshReconstructOnline.cpp
        ./src/MeshFragment.cpp
        ./src/NonPlaneMesh.cpp
        ./src/GL_gui.cpp

        ./include/ikd-Tree/ikd_Tree.cpp

#        ./src/meshing/mesh_rec_display.cpp
        ./src/tools/openGL_libs/openGL_camera.cpp
        ./src/tools/openGL_libs/openGL_camera_view.cpp
        ./src/tools/openGL_libs/glad.c

        ./src/tools/imgui/imgui.cpp
        ./src/tools/imgui/imgui_draw.cpp
        ./src/tools/imgui/imgui_impl_glfw.cpp
        ./src/tools/imgui/imgui_impl_opengl3.cpp
        ./src/tools/imgui/imgui_tables.cpp
        ./src/tools/imgui/imgui_widgets.cpp
        )
target_link_libraries(my_mesh
        ${catkin_LIBRARIES}
        ${PCL_LIBRARIES}
        ${PYTHON_LIBRARIES}
        CGAL::CGAL
        tbb
        ${OpenCV_LIBRARIES}
        ${Boost_LIBRARIES}
        ${Boost_SYSTEM_LIBRARIES}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SERIALIZATION_LIBRARY} # serialization
        ${CERES_LIBRARIES}
        pcl_common
        pcl_io
        pcl_filters
        pcl_kdtree
        GL GLU glfw X11 Xrandr pthread Xi dl Xinerama Xcursor)
target_include_directories(my_mesh PRIVATE
        ${PYTHON_INCLUDE_DIRS}
        # ===== Meshing =====
        ./include
        ./include/ikd-Tree/
        ./src/meshing/
        ./src/meshing//r3live/
        ./src/meshing//kd_tree_dev/
        ./src/meshing//rgb_map/
        ./src/tools/
        ./src/tools/imgui
        ./src/tools/shader
        ./src/tools/openGL_libs
        ./src)

add_subdirectory(src/odometry/voxelmap)
add_subdirectory(src/odometry/fastlio)
add_subdirectory(src/odometry/gt_odometry)
