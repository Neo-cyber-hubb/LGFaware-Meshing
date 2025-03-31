
add_executable(voxel_mapping_odom
        ./src/voxelMapping.cpp
        ./src/preprocess.cpp)
target_link_libraries(voxel_mapping_odom ${catkin_LIBRARIES} ${PCL_LIBRARIES})
target_include_directories(voxel_mapping_odom
        PRIVATE ${PYTHON_INCLUDE_DIRS}
        ./include)


