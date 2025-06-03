//
// Created by neo on 2024/4/7.
//

#ifndef SRC_DATATYPE_H
#define SRC_DATATYPE_H

#include <iostream>

typedef pcl::PointXYZINormal PointType;
typedef pcl::PointCloud<PointType> PointCloudXYZI;
typedef pcl::PointCloud<pcl::PointXYZ> PointCloudXYZ;


// hash体素
struct m_voxel {
    std::vector<int>                    all_pts_idx;    // 体素内点的索引
    std::vector<std::array<int, 3>>     face_list;      //
    int                                 new_pt_num;
    bool                                if_need_update; // 为了避免重复添加
    m_voxel()
    {
        if_need_update = false;
        new_pt_num = 0;
    }
};

#endif //SRC_DATATYPE_H
