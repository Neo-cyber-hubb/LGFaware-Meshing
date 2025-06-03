
find_package(PythonLibs REQUIRED)
find_path(MATPLOTLIB_CPP_INCLUDE_DIRS "matplotlibcpp.h")

add_executable(fastlio_mapping
        ./src/laserMapping.cpp
        ./include/ikd-Tree/ikd_Tree.cpp
        ./src/preprocess.cpp)
target_link_libraries(fastlio_mapping ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${PYTHON_LIBRARIES})
target_include_directories(fastlio_mapping
        PRIVATE ${PYTHON_INCLUDE_DIRS}
        ./include)