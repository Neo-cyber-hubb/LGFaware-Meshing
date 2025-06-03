

add_executable(gt_odom gt_pub.cpp)
target_link_libraries(gt_odom ${catkin_LIBRARIES} ${PCL_LIBRARIES}
        pcl_common
        pcl_io)
target_include_directories(gt_odom PRIVATE ${PYTHON_INCLUDE_DIRS})