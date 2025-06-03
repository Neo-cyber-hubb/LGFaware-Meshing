//
// Created by neo on 2024/4/7.
//

#include "NonPlaneMesh.h"

#include <tbb/tbb.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>

// 从小到大排序
static bool compare(const std::pair<int, double>& a, const std::pair<int, double>& b)
{
    return a.second < b.second;
}

// 从大到小排序
static bool compare_1(const std::pair<int, double>& a, const std::pair<int, double>& b)
{
    return a.second > b.second;
}

void NonPlaneMesh::search_nearest_vertex_under_distance_r(const int& idx, std::vector<int>& pt_vec, std::vector<double>& dis_vec)
{
    std::unordered_map<int, double> idx_dis_dic;
    Eigen::Vector3d pt_b;
    pt_b << pts_list->points[idx].x, pts_list->points[idx].y, pts_list->points[idx].z;
    long voxel_x_min = std::round((pt_b(0) - r_dis) / voxel_resolution);
    long voxel_x_max = std::round((pt_b(0) + r_dis) / voxel_resolution);
    long voxel_y_min = std::round((pt_b(1) - r_dis) / voxel_resolution);
    long voxel_y_max = std::round((pt_b(1) + r_dis) / voxel_resolution);
    long voxel_z_min = std::round((pt_b(2) - r_dis) / voxel_resolution);
    long voxel_z_max = std::round((pt_b(2) + r_dis) / voxel_resolution);
    for (long i = voxel_x_min; i <= voxel_x_max; i++)
    {
        for (long j = voxel_y_min; j <= voxel_y_max; j++)
        {
            for (long k = voxel_z_min; k <= voxel_z_max; k++)
            {

                if (!voxel_pts.if_exist(i, j, k))
                {
                    continue;
                }

                for (auto it : voxel_pts.m_map_3d_hash_map[i][j][k]->all_pts_idx)
                {
                    if (it == idx) { continue; }
                    // 若不是三角形顶点则跳过
                    if (pts_state[it] != 1 && pts_state[it] != 2)
                    {
                        continue;
                    }
                    Eigen::Vector3d pt_a;
                    pt_a << pts_list->points[it].x, pts_list->points[it].y, pts_list->points[it].z;

                    double dis = (pt_a - pt_b).norm();
                    if (dis <= r_dis /*&& ang >= 0*/)
                    {
                        // 距离小于阈值，法向夹角不超过90度
                        idx_dis_dic[it] = dis;
                    }
                }
            }
        }
    }
    // 按距离从小到大排序
    if (idx_dis_dic.size() > 0)
    {
        std::vector<std::pair<int, double>> indexVec(idx_dis_dic.begin(), idx_dis_dic.end());
        std::sort(indexVec.begin(), indexVec.end(), compare);
        for (const auto& pair : indexVec)
        {
            pt_vec.push_back(pair.first);
            dis_vec.push_back(pair.second);
        }
    }
}

void NonPlaneMesh::search_edge_vertex_under_distance_r(const int& idx, std::vector<int>& vec)
{
    Eigen::Vector3d pt_b;
    pt_b << pts_list->points[idx].x, pts_list->points[idx].y, pts_list->points[idx].z;
    long voxel_x_min = std::round((pt_b(0) - r_dis) / voxel_resolution);
    long voxel_x_max = std::round((pt_b(0) + r_dis) / voxel_resolution);
    long voxel_y_min = std::round((pt_b(1) - r_dis) / voxel_resolution);
    long voxel_y_max = std::round((pt_b(1) + r_dis) / voxel_resolution);
    long voxel_z_min = std::round((pt_b(2) - r_dis) / voxel_resolution);
    long voxel_z_max = std::round((pt_b(2) + r_dis) / voxel_resolution);
    for (long i = voxel_x_min; i <= voxel_x_max; i++)
    {
        for (long j = voxel_y_min; j <= voxel_y_max; j++)
        {
            for (long k = voxel_z_min; k <= voxel_z_max; k++)
            {

                if (!voxel_pts.if_exist(i, j, k))
                {
                    continue;
                }

                for (auto it : voxel_pts.m_map_3d_hash_map[i][j][k]->all_pts_idx)
                {
                    if (it == idx) { continue; }
                    // 若不是三角形顶点则跳过
                    if (pts_state[it] != 2)
                    {
                        continue;
                    }
                    Eigen::Vector3d pt_a;
                    pt_a << pts_list->points[it].x, pts_list->points[it].y, pts_list->points[it].z;

                    double dis = (pt_a - pt_b).norm();
                    if (dis <= r_dis /*&& ang >= 0*/)
                    {
                        // 距离小于阈值，法向夹角不超过90度
                        vec.push_back(it);
                    }
                }
            }
        }
    }
}

void NonPlaneMesh::search_edge_vertex_under_distance_kr(const int& idx, double max_dis, std::vector<int>& vec)
{
    Eigen::Vector3d pt_b;
    pt_b << pts_list->points[idx].x, pts_list->points[idx].y, pts_list->points[idx].z;
    long voxel_x_min = std::round((pt_b(0) - max_dis) / voxel_resolution);
    long voxel_x_max = std::round((pt_b(0) + max_dis) / voxel_resolution);
    long voxel_y_min = std::round((pt_b(1) - max_dis) / voxel_resolution);
    long voxel_y_max = std::round((pt_b(1) + max_dis) / voxel_resolution);
    long voxel_z_min = std::round((pt_b(2) - max_dis) / voxel_resolution);
    long voxel_z_max = std::round((pt_b(2) + max_dis) / voxel_resolution);
    for (long i = voxel_x_min; i <= voxel_x_max; i++)
    {
        for (long j = voxel_y_min; j <= voxel_y_max; j++)
        {
            for (long k = voxel_z_min; k <= voxel_z_max; k++)
            {

                if (!voxel_pts.if_exist(i, j, k))
                {
                    continue;
                }

                for (auto it : voxel_pts.m_map_3d_hash_map[i][j][k]->all_pts_idx)
                {
                    if (it == idx) { continue; }
                    // 若不是三角形顶点则跳过
                    if (pts_state[it] != 2)
                    {
                        continue;
                    }
                    Eigen::Vector3d pt_a;
                    pt_a << pts_list->points[it].x, pts_list->points[it].y, pts_list->points[it].z;

                    double dis = (pt_a - pt_b).norm();
                    if (dis <= max_dis /*&& ang >= 0*/)
                    {
                        // 距离小于阈值，法向夹角不超过90度
                        vec.push_back(it);
                    }
                }
            }
        }
    }
}

void NonPlaneMesh::search_nearest_point_under_distance_r(const int& idx, std::vector<int>& vec)
{
    Eigen::Vector3d pt_b;
    pt_b << pts_list->points[idx].x, pts_list->points[idx].y, pts_list->points[idx].z;
    long voxel_x_min = std::round((pt_b(0) - r_dis) / voxel_resolution);
    long voxel_x_max = std::round((pt_b(0) + r_dis) / voxel_resolution);
    long voxel_y_min = std::round((pt_b(1) - r_dis) / voxel_resolution);
    long voxel_y_max = std::round((pt_b(1) + r_dis) / voxel_resolution);
    long voxel_z_min = std::round((pt_b(2) - r_dis) / voxel_resolution);
    long voxel_z_max = std::round((pt_b(2) + r_dis) / voxel_resolution);
    for (long i = voxel_x_min; i <= voxel_x_max; i++)
    {
        for (long j = voxel_y_min; j <= voxel_y_max; j++)
        {
            for (long k = voxel_z_min; k <= voxel_z_max; k++)
            {
                if (!voxel_pts.if_exist(i, j, k))
                {
                    continue;
                }

                for (auto it : voxel_pts.m_map_3d_hash_map[i][j][k]->all_pts_idx)
                {
                    if (it == idx) { continue; }
                    // 若不是孤立点则跳过
                    if (pts_state[it] != 0)
                    {
                        continue;
                    }
                    Eigen::Vector3d pt_a;
                    pt_a << pts_list->points[it].x, pts_list->points[it].y, pts_list->points[it].z;

                    double dis = (pt_a - pt_b).norm();
                    if (dis <= minimum_pt_dis /*&& ang >= 0*/)
                    {
                        // 距离小于阈值，法向夹角不超过90度
                        vec.push_back(it);
                    }
                }
            }
        }
    }
}


void NonPlaneMesh::search_nearest_valid_ptv(const Eigen::Vector3d& center_pt, double max_dis, std::vector<int>& pt_vec)
{
    long voxel_x_min = std::round((center_pt(0) - max_dis) / voxel_resolution);
    long voxel_x_max = std::round((center_pt(0) + max_dis) / voxel_resolution);
    long voxel_y_min = std::round((center_pt(1) - max_dis) / voxel_resolution);
    long voxel_y_max = std::round((center_pt(1) + max_dis) / voxel_resolution);
    long voxel_z_min = std::round((center_pt(2) - max_dis) / voxel_resolution);
    long voxel_z_max = std::round((center_pt(2) + max_dis) / voxel_resolution);
    for (long i = voxel_x_min; i <= voxel_x_max; i++)
    {
        for (long j = voxel_y_min; j <= voxel_y_max; j++)
        {
            for (long k = voxel_z_min; k <= voxel_z_max; k++)
            {
                if (!voxel_pts.if_exist(i, j, k))
                {
                    continue;
                }

                for (auto it : voxel_pts.m_map_3d_hash_map[i][j][k]->all_pts_idx)
                {
                    // 若不是有效点则跳过
                    if (pts_state[it] != 1 && pts_state[it] != 2 && pts_state[it] != 0)
                    {
                        continue;
                    }
                    Eigen::Vector3d pt_a;
                    pt_a << pts_list->points[it].x, pts_list->points[it].y, pts_list->points[it].z;

                    double dis = (pt_a - center_pt).norm();
                    if (dis <= max_dis /*&& ang >= 0*/)
                    {
                        // 距离小于阈值，法向夹角不超过90度
                        pt_vec.push_back(it);
                    }
                }
            }
        }
    }
}



bool NonPlaneMesh::if_projected_inner(const Eigen::Vector3d& Q, const Eigen::Vector3d& A,
                                      const Eigen::Vector3d& B, const Eigen::Vector3d& C,
                                      Eigen::Vector3d P )
{
    // P是lidar的位置, Q是点的位置
    Eigen::Vector3d ab = B - A;
    Eigen::Vector3d ac = C - A;
    Eigen::Vector3d qp = P - Q;
    Eigen::Vector3d norm = ab.cross(ac);   // 三角形顶点逆时针排序

    Eigen::Vector3d bc = C - B;
    Eigen::Vector3d aq = Q - A;
    Eigen::Vector3d bq = Q - B;
    Eigen::Vector3d norm_ab = norm.cross(ab);
    if (norm_ab.dot(aq) <= 0) { return false; }

    Eigen::Vector3d norm_bc = norm.cross(bc);
    if (norm_bc.dot(bq) <= 0) { return false; }

    Eigen::Vector3d norm_ca = norm.cross(-1 * ac);
    if (norm_ca.dot(aq) <= 0) { return false; }

    // TODO 判断一下点到三角形的距离，不能过大

    return true;
}

bool NonPlaneMesh::if_p_in_tri(const Eigen::Vector3d& Q, const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C)
{
    // ABC任意顺序
    Eigen::Vector3d ab = B - A;
    Eigen::Vector3d ac = C - A;
    Eigen::Vector3d norm = ab.cross(ac);

    Eigen::Vector3d bc = C - B;
    Eigen::Vector3d aq = Q - A;
    Eigen::Vector3d bq = Q - B;
    Eigen::Vector3d norm_ab = norm.cross(ab);
    if (norm_ab.dot(aq) <= 0) { return false; }

    Eigen::Vector3d norm_bc = norm.cross(bc);
    if (norm_bc.dot(bq) <= 0) { return false; }

    Eigen::Vector3d norm_ca = norm.cross(-1 * ac);
    if (norm_ca.dot(aq) <= 0) { return false; }

    // TODO 判断一下点到三角形的距离，不能过大

    return true;
}

bool NonPlaneMesh::if_edge_vertex(const int& pt_idx)
{
    // 顶点对应的三角形索引列表不能为空
    assert(pt2face_idx[pt_idx].size() != 0);

    std::unordered_map<int, int> vertex_idx_num;
    for (auto it : pt2face_idx[pt_idx])
    {
        for (int i = 0; i < 3; i++)
        {
            int key = faces_list_all[it][i];

            if (key == pt_idx) { continue; }

            if (vertex_idx_num.find(key) == vertex_idx_num.end())
            {
                vertex_idx_num[key] = 1;
            }
            else
            {
                vertex_idx_num[key]++;
            }
        }
    }

    for (auto it = vertex_idx_num.begin(); it != vertex_idx_num.end(); ++it) {
        if (it->second == 1)
        {
            return true;
        }
    }
    return false;
}

void NonPlaneMesh::find_edge_vertex(const int& pt_idx, std::vector<int> &adj_edge_v)
{
    if (pt2face_idx[pt_idx].size() > 0)
    {
        adj_edge_v.clear();
        std::unordered_map<int, int> vertex_idx_num;
        for (auto it : pt2face_idx[pt_idx])
        {
            for (int i = 0; i < 3; i++)
            {
                int key = faces_list_all[it][i];

                if (key == pt_idx) { continue; }

                if (vertex_idx_num.find(key) == vertex_idx_num.end())
                {
                    vertex_idx_num[key] = 1;
                }
                else
                {
                    vertex_idx_num[key]++;
                }
            }
        }

        for (auto it = vertex_idx_num.begin(); it != vertex_idx_num.end(); ++it) {
            if (it->second == 1)
            {
                adj_edge_v.push_back(it->first);
            }
        }
    }
}

void NonPlaneMesh::find_edge_and_abnormal_vertex(const int& pt_idx, std::vector<int> &adj_edge_v, std::vector<int> &abnormal_v)
{
    // 顶点对应的三角形索引列表不能为空
    assert(pt2face_idx[pt_idx].size() != 0);

    adj_edge_v.clear();
    abnormal_v.clear();
    std::unordered_map<int, int> vertex_idx_num;
    for (auto it : pt2face_idx[pt_idx])
    {
        for (int i = 0; i < 3; i++)
        {
            int key = faces_list_all[it][i];

            if (key == pt_idx) { continue; }

            if (vertex_idx_num.find(key) == vertex_idx_num.end())
            {
                vertex_idx_num[key] = 1;
            }
            else
            {
                vertex_idx_num[key]++;
            }
        }
    }

    for (auto it = vertex_idx_num.begin(); it != vertex_idx_num.end(); ++it) {
        if (it->second == 1)
        {
            adj_edge_v.push_back(it->first);
        }
        if (it->second > 2)
        {
            abnormal_v.push_back(it->first);
        }
    }
}

bool NonPlaneMesh::if_reasonable_normal(const int& A_idx, const int& B_idx, const int& C_idx)
{
    Eigen::Vector3d A, B, C, normA, normB, normC;
    A << pts_list->points[A_idx].x, pts_list->points[A_idx].y,
            pts_list->points[A_idx].z;
    B << pts_list->points[B_idx].x, pts_list->points[B_idx].y,
            pts_list->points[B_idx].z;
    C << pts_list->points[C_idx].x, pts_list->points[C_idx].y,
            pts_list->points[C_idx].z;
    normA << pts_list->points[A_idx].normal_x, pts_list->points[A_idx].normal_y,
            pts_list->points[A_idx].normal_z;
    normB << pts_list->points[B_idx].normal_x, pts_list->points[B_idx].normal_y,
            pts_list->points[B_idx].normal_z;
    normC << pts_list->points[C_idx].normal_x, pts_list->points[C_idx].normal_y,
            pts_list->points[C_idx].normal_z;

    Eigen::Vector3d ab = B - A;
    Eigen::Vector3d ac = C - A;
    Eigen::Vector3d norm = ab.cross(ac);   // 三角形顶点逆时针排序

    if (norm.dot(normA) <= 0)
    {
        return false;
    }
    if (norm.dot(normB) <= 0)
    {
        return false;
    }
    if (norm.dot(normC) <= 0)
    {
        return false;
    }
    return true;
}

bool NonPlaneMesh::if_different_side(const Eigen::Vector3d& edge_A, const Eigen::Vector3d& edge_B,
                                     const Eigen::Vector3d& P, const Eigen::Vector3d& Q)
{
    // 沿lidar-point投影还是沿法向投影
    Eigen::Vector3d ab = edge_B - edge_A;
    Eigen::Vector3d ac = P - edge_A;
//    Eigen::Vector3d n = ab.cross(lidar_beam);
    Eigen::Vector3d n = ab.cross(ac);
    Eigen::Vector3d n_ab = ab.cross(n);
    Eigen::Vector3d ap = P - edge_A;
    Eigen::Vector3d aq = Q - edge_A;
    double d1 = n_ab.dot(ap);
    double d2 = n_ab.dot(aq);
    if (d1 * d2 < 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}


void NonPlaneMesh::fit_2_tri(std::vector<int> &pts_idx_list, Eigen::Vector3d &avg_normal,
                             std::vector<std::array<int, 3>> &tri_pt_indices)
{
    int             pt_size = pts_idx_list.size();
    Eigen::MatrixXd pc_mat;
    pc_mat.resize( pt_size, 3 );
    // 计算坐标变换矩阵
    for ( int i = 0; i < pts_idx_list.size(); i++ )
    {
        int temp_idx = pts_idx_list[i];
        pc_mat(i, 0) = pts_list->points[temp_idx].x;
        pc_mat(i, 1) = pts_list->points[temp_idx].y;
        pc_mat(i, 2) = pts_list->points[temp_idx].z;
    }
    Eigen::Vector3d pc_center = pc_mat.colwise().mean().transpose();
    Eigen::MatrixXd pt_sub_center = pc_mat.rowwise() - pc_center.transpose();
    Eigen::Matrix3d cov = ( pt_sub_center.transpose() * pt_sub_center ) / double( pc_mat.rows() );

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver;
    eigensolver.compute( cov );

    Eigen::Vector3d short_axis = eigensolver.eigenvectors().col( 0 );
    Eigen::Vector3d mid_axis = eigensolver.eigenvectors().col( 1 );
    Eigen::Vector3d eigenvalues = eigensolver.eigenvalues();
    if ( avg_normal.dot( short_axis ) < 0 )
    {
        short_axis *= -1;
    }
    Eigen::Vector3d long_axis = short_axis.cross( mid_axis );

    Eigen::MatrixXd pro_pc_1(pt_size, 1);
    Eigen::MatrixXd pro_pc_2(pt_size, 1);
    pro_pc_1 = pt_sub_center * mid_axis;
    pro_pc_2 = pt_sub_center * long_axis;


    Eigen::MatrixXd new_pc_mat;
    new_pc_mat.resize( 4, 3 );

    new_pc_mat.row(0) = pc_center + pro_pc_1.minCoeff() * mid_axis;
    new_pc_mat.row(1) = pc_center + pro_pc_1.maxCoeff() * mid_axis;
    new_pc_mat.row(2) = pc_center + pro_pc_2.minCoeff() * long_axis;
    new_pc_mat.row(3) = pc_center + pro_pc_2.maxCoeff() * long_axis;

    // 将新点加入到点列表中
    std::vector<int> new_pc_idx;
    for (int i = 0; i < 4; i++)
    {
        pcl::PointXYZINormal temp_pt;
        temp_pt.x = new_pc_mat(i, 0);
        temp_pt.y = new_pc_mat(i, 1);
        temp_pt.z = new_pc_mat(i, 2);

        mtx_pts_list.lock();
        int temp_cur_idx = pts_list->size();
        new_pc_idx.push_back(temp_cur_idx);
        pts_list->push_back(temp_pt);
        mtx_pts_list.unlock();

        mtx_pts_state.lock();
        pts_state.push_back(2);
        mtx_pts_state.unlock();

        mtx_pts_process_num.lock();
        pts_process_num.push_back(1);
        hole_process_num.push_back(1);
        hole_v_process.push_back(1);
        mtx_pts_process_num.unlock();

        mtx_tri_edge_v_list_all.lock();
        tri_edge_v_list_all.push_back(temp_cur_idx);
        mtx_tri_edge_v_list_all.unlock();

        mtx_pt2face_idx.lock();
        pt2face_idx.resize(pts_list->size());
        mtx_pt2face_idx.unlock();
    }
    // 确定三角形和顶点顺序
    Eigen::Vector3d ab, ac, norm;
    ab = new_pc_mat.row(1) - new_pc_mat.row(0);
    ac = new_pc_mat.row(2) - new_pc_mat.row(0);
    norm = ab.cross(ac);
    if (norm.dot(short_axis) > 0)
    {
        std::array<int, 3> temp_tri;
        temp_tri[0] = new_pc_idx[0];
        temp_tri[1] = new_pc_idx[1];
        temp_tri[2] = new_pc_idx[2];
        tri_pt_indices.push_back(temp_tri);
    }
    else
    {
        std::array<int, 3> temp_tri;
        temp_tri[0] = new_pc_idx[0];
        temp_tri[1] = new_pc_idx[2];
        temp_tri[2] = new_pc_idx[1];
        tri_pt_indices.push_back(temp_tri);
    }

    ac = new_pc_mat.row(3) - new_pc_mat.row(0);
    norm = ab.cross(ac);
    if (norm.dot(short_axis) > 0)
    {
        std::array<int, 3> temp_tri;
        temp_tri[0] = new_pc_idx[0];
        temp_tri[1] = new_pc_idx[1];
        temp_tri[2] = new_pc_idx[3];
        tri_pt_indices.push_back(temp_tri);
    }
    else
    {
        std::array<int, 3> temp_tri;
        temp_tri[0] = new_pc_idx[0];
        temp_tri[1] = new_pc_idx[3];
        temp_tri[2] = new_pc_idx[1];
        tri_pt_indices.push_back(temp_tri);
    }
    // 将新点添加到体素索引中
    mtx_voxel.lock();
    for (int i = 0; i < 4; i++)
    {
        long voxel_x = std::round(new_pc_mat(i, 0) / voxel_resolution);
        long voxel_y = std::round(new_pc_mat(i, 1) / voxel_resolution);
        long voxel_z = std::round(new_pc_mat(i, 2) / voxel_resolution);
        if (voxel_pts.if_exist(voxel_x, voxel_y, voxel_z))
        {
            voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->all_pts_idx.push_back(new_pc_idx[i]);
        }
        else
        {
            std::shared_ptr<m_voxel> temp_voxel = std::make_shared<m_voxel>();
            voxel_pts.insert(voxel_x, voxel_y, voxel_z, temp_voxel);
            voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->all_pts_idx.push_back(new_pc_idx[i]);
        }
    }
    mtx_voxel.unlock();
}

/*
 * 要排序的点索引列表、中心点索引、三角形除中心点外的另外两个点、排序后的点索引列表、排序后的点对应的像限索引
 */
void NonPlaneMesh::anticlockwise_sort(std::vector<int> &pts_idx_list, int ref_pt_idx,
                                      std::vector<std::vector<int>> &tri_v_idx,
                                      std::vector<int> &pts_idx_list_sorted,
                                      std::vector<int> &pts_axis_idx)
{
    pts_idx_list_sorted.clear();
    pts_axis_idx.clear();
    std::vector<int> pts_idx_list_sorted_initial;
    std::vector<int> pts_axis_idx_initial;
    int             pt_size = pts_idx_list.size() + 1;
    Eigen::MatrixXd pc_mat;
    pc_mat.resize( pt_size, 3 );
    // 计算坐标变换矩阵
    for ( int i = 0; i < pts_idx_list.size(); i++ )
    {
        int temp_idx = pts_idx_list[i];
        pc_mat(i, 0) = pts_list->points[temp_idx].x;
        pc_mat(i, 1) = pts_list->points[temp_idx].y;
        pc_mat(i, 2) = pts_list->points[temp_idx].z;
    }
    pc_mat(pt_size-1, 0) = pts_list->points[ref_pt_idx].x;
    pc_mat(pt_size-1, 1) = pts_list->points[ref_pt_idx].y;
    pc_mat(pt_size-1, 2) = pts_list->points[ref_pt_idx].z;

    Eigen::Vector3d pc_center = pc_mat.colwise().mean().transpose();
    Eigen::MatrixXd pt_sub_center = pc_mat.rowwise() - pc_center.transpose();
    Eigen::Matrix3d cov = ( pt_sub_center.transpose() * pt_sub_center ) / double( pc_mat.rows() );

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver;
    eigensolver.compute( cov );

    Eigen::Vector3d short_axis = eigensolver.eigenvectors().col( 0 );
    Eigen::Vector3d mid_axis = eigensolver.eigenvectors().col( 1 );
    Eigen::Vector3d long_axis = short_axis.cross( mid_axis );

    // 投影坐标
    Eigen::MatrixXd pro_pc(pt_size, 2);
    pro_pc.col(0) = pt_sub_center * mid_axis;
    pro_pc.col(1) = pt_sub_center * long_axis;
    pro_pc.rowwise() -= pro_pc.row(pt_size-1);  // 以参考点为中心

    // 分像限
    std::vector<std::vector<int>> axis_pts_idx;
    axis_pts_idx.resize(6);
    for (int i = 0; i < pt_size; i++)
    {
        if (pro_pc(i, 0) == 0)
        {
            if (pro_pc(i, 1) > 0)
            {
                axis_pts_idx[1].push_back(i);
            }
            else if (pro_pc(i, 1) < 0)
            {
                axis_pts_idx[4].push_back(i);
            }
            else
            {
                //
            }
        }
        else if (pro_pc(i, 0) > 0)
        {
            if (pro_pc(i, 1) >= 0)
            {
                axis_pts_idx[0].push_back(i);
            }
            else
            {
                axis_pts_idx[5].push_back(i);
            }
        }
        else
        {
            // < 0
            if (pro_pc(i, 1) >= 0)
            {
                axis_pts_idx[2].push_back(i);
            }
            else
            {
                axis_pts_idx[3].push_back(i);
            }
        }
    }

    // 分别排序 tan
    for (int i = 0; i < axis_pts_idx.size(); i++)
    {
        std::unordered_map<int, double> idx_tan_dic;
        if (axis_pts_idx[i].size() == 0) { continue; }
        for (int j = 0; j < axis_pts_idx[i].size(); j++)
        {
            if (i == 1 || i == 4)
            {
                pts_idx_list_sorted_initial.push_back(pts_idx_list[axis_pts_idx[i][j]]);
                pts_axis_idx_initial.push_back(i);
                continue;
            }
            int temp_idx = axis_pts_idx[i][j];
            idx_tan_dic[pts_idx_list[temp_idx]] = pro_pc(temp_idx, 1) / pro_pc(temp_idx, 0);
        }
        if (i == 1 || i == 4) { continue; }
        std::vector<std::pair<int, double>> indexVec(idx_tan_dic.begin(), idx_tan_dic.end());
        std::sort(indexVec.begin(), indexVec.end(), compare);
        for (const auto& pair : indexVec)
        {
            pts_idx_list_sorted_initial.push_back(pair.first);
            pts_axis_idx_initial.push_back(i);
        }
    }

    // 调整一下初始点位置，保证初始轴不穿过三角形
    bool if_first_is_tri = false;
    int first_v = pts_idx_list_sorted_initial[0];
    int second_v = pts_idx_list_sorted_initial[1];
    std::vector<int> ref_tri_list;
    for (int i = 0; i < tri_v_idx.size(); i++)
    {
        for (int j = 0; j < tri_v_idx[i].size(); j++)
        {
            if (tri_v_idx[i][j] == first_v)
            {
                int j_next = (j+1)%tri_v_idx[i].size();
                if (tri_v_idx[i][j_next] == second_v)
                {
                    if_first_is_tri = true;
                    break;
                }
            }
        }
        if (if_first_is_tri)
        {
            break;
        }
    }

    if (if_first_is_tri)
    {
        for (int i = 0; i < pts_idx_list_sorted_initial.size(); i++)
        {
            pts_idx_list_sorted.push_back(pts_idx_list_sorted_initial[i]);
            pts_axis_idx.push_back(pts_axis_idx_initial[i]);
        }
    }
    else
    {
        for (int i = 1; i < pts_idx_list_sorted_initial.size(); i++)
        {
            pts_idx_list_sorted.push_back(pts_idx_list_sorted_initial[i]);
            pts_axis_idx.push_back(pts_axis_idx_initial[i]);
        }
        pts_idx_list_sorted.push_back(pts_idx_list_sorted_initial[0]);
        pts_axis_idx.push_back(pts_axis_idx_initial[0]);
    }

}


bool NonPlaneMesh::if_lay_in_hole(int pt_idx, std::vector<int>& adj_v_list, std::vector<int>& all_edge_v)
{
    all_edge_v.clear();

    int initial_adj_v1, initial_adj_v2;
    initial_adj_v1 = adj_v_list[0];
    if (adj_v_list.size() == 2)
    {
        initial_adj_v2 = adj_v_list[1];
    }
    // 两个缺口，选一个处理
    else if (adj_v_list.size() == 4)
    {
        // 目标点相连的所有顶点
        std::set<int> circle_tri_v;
        std::vector<std::vector<int>> tri_v_idx;
        tri_v_idx.resize(pt2face_idx[pt_idx].size());
        for (int j = 0; j < pt2face_idx[pt_idx].size(); j++)
        {
            int f_idx = pt2face_idx[pt_idx][j];
            for (int l = 0; l < 3; l++)
            {
                if (faces_list_all[f_idx][l] == pt_idx) { continue; }
                circle_tri_v.insert(faces_list_all[f_idx][l]);
                tri_v_idx[j].push_back(faces_list_all[f_idx][l]);
            }
        }

        // 点排序
        std::vector<int> cir_tri_v, cir_tri_v_sorted, pts_to_axis;
        for (auto it : circle_tri_v)
        {
            cir_tri_v.push_back(it);
        }
        anticlockwise_sort(cir_tri_v, pt_idx, tri_v_idx, cir_tri_v_sorted, pts_to_axis);

        // 判断next相邻的前后点是否一个边缘点一个不是边缘点
        int sorted_adj_v1, sorted_adj_v2;
        sorted_adj_v1 = -1;
        for (int j = 0; j < cir_tri_v_sorted.size(); j++)
        {
            if (cir_tri_v_sorted[j] == initial_adj_v1)
            {
                int n = cir_tri_v_sorted.size();
                sorted_adj_v1 = cir_tri_v_sorted[(j-1+n)%n];
                sorted_adj_v2 = cir_tri_v_sorted[(j+1+n)%n];
            }
        }

        // 某些地方有问题
        if (sorted_adj_v1 == -1)
        {
            return false;
        }

        // 相邻的边缘点就是要找的点
        auto sorted_adj_v1_idx = std::find(adj_v_list.begin(), adj_v_list.end(), sorted_adj_v1);
        if (sorted_adj_v1_idx == adj_v_list.end())
        {
            auto sorted_adj_v2_idx = std::find(adj_v_list.begin(), adj_v_list.end(), sorted_adj_v2);
            if (sorted_adj_v2_idx == adj_v_list.end())
            {
                return false;
            }
            else
            {
                initial_adj_v2 = sorted_adj_v2;
            }
        }
        else
        {
            auto sorted_adj_v2_idx = std::find(adj_v_list.begin(), adj_v_list.end(), sorted_adj_v2);
            if (sorted_adj_v2_idx == adj_v_list.end())
            {
                initial_adj_v2 = sorted_adj_v1;
            }
            else
            {
                return false;
            }
        }

    }
    else
    {
        return false;
    }

    // 查找边缘点的相邻边缘点
    all_edge_v.push_back(initial_adj_v1);
    all_edge_v.push_back(pt_idx);
    int last_edge_v = pt_idx;
    int next_edge_v = initial_adj_v2;

    while(next_edge_v != all_edge_v[0] && all_edge_v.size() <= 10)
    {
        all_edge_v.push_back(next_edge_v);
        std::vector<int> temp_adj_v_list;
        find_edge_vertex(next_edge_v, temp_adj_v_list);
        if (temp_adj_v_list.size() == 2)
        {
            if (temp_adj_v_list[0] == last_edge_v)
            {
                last_edge_v = next_edge_v;
                next_edge_v = temp_adj_v_list[1];
            }
            else
            {
                last_edge_v = next_edge_v;
                next_edge_v = temp_adj_v_list[0];
            }
        }
        else if (temp_adj_v_list.size() == 4)
        {
            // 目标点相连的所有顶点
            std::set<int> circle_tri_v;
            std::vector<std::vector<int>> tri_v_idx;
            tri_v_idx.resize(pt2face_idx[next_edge_v].size());
            for (int j = 0; j < pt2face_idx[next_edge_v].size(); j++)
            {
                int f_idx = pt2face_idx[next_edge_v][j];
                for (int l = 0; l < 3; l++)
                {
                    if (faces_list_all[f_idx][l] == next_edge_v) { continue; }
                    circle_tri_v.insert(faces_list_all[f_idx][l]);
                    tri_v_idx[j].push_back(faces_list_all[f_idx][l]);
                }
            }

            // 点排序
            std::vector<int> cir_tri_v, cir_tri_v_sorted, pts_to_axis;
            for (auto it : circle_tri_v)
            {
                cir_tri_v.push_back(it);
            }
            anticlockwise_sort(cir_tri_v, next_edge_v, tri_v_idx, cir_tri_v_sorted, pts_to_axis);

            // 判断next相邻的前后点是否一个边缘点一个不是边缘点
            int sorted_adj_v1, sorted_adj_v2;
            sorted_adj_v1 = -1;
            for (int j = 0; j < cir_tri_v_sorted.size(); j++)
            {
                if (cir_tri_v_sorted[j] == last_edge_v)
                {
                    int n = cir_tri_v_sorted.size();
                    sorted_adj_v1 = cir_tri_v_sorted[(j-1+n)%n];
                    sorted_adj_v2 = cir_tri_v_sorted[(j+1+n)%n];
                }
            }

            // 某些地方有问题
            if (sorted_adj_v1 == -1) { break; }

            // 相邻的边缘点就是要找的点
            auto sorted_adj_v1_idx = std::find(temp_adj_v_list.begin(), temp_adj_v_list.end(), sorted_adj_v1);
            if (sorted_adj_v1_idx == temp_adj_v_list.end())
            {
                auto sorted_adj_v2_idx = std::find(temp_adj_v_list.begin(), temp_adj_v_list.end(), sorted_adj_v2);
                if (sorted_adj_v2_idx == temp_adj_v_list.end())
                {
                    break;
                }
                else
                {
                    last_edge_v = next_edge_v;
                    next_edge_v = sorted_adj_v2;
                }
            }
            else
            {
                auto sorted_adj_v2_idx = std::find(temp_adj_v_list.begin(), temp_adj_v_list.end(), sorted_adj_v2);
                if (sorted_adj_v2_idx == temp_adj_v_list.end())
                {
                    last_edge_v = next_edge_v;
                    next_edge_v = sorted_adj_v1;
                }
                else
                {
                    break;
                }
            }
        }
        else
        {
            break;
        }
    }

    if (next_edge_v != all_edge_v[0])
    {
        return false;
    }
    else
    {
        return true;
    }

}

// 递归处理
void NonPlaneMesh::fill_hole(std::vector<int>& hole_pts_list)
{
    std::vector<int> added_tri_list;

    // 如果只有三个点
    if (hole_pts_list.size() == 3)
    {
        // 根据邻接三角形确定法向
        int adj_tri_idx = -1;

        if (adj_tri_idx == -1)
        {
            adj_tri_idx = pt2face_idx[hole_pts_list[1]][0];
        }

        assert(adj_tri_idx != -1);
        int tri_A_idx = faces_list_all[adj_tri_idx][0];
        int tri_B_idx = faces_list_all[adj_tri_idx][1];
        int tri_C_idx = faces_list_all[adj_tri_idx][2];
        Eigen::Vector3d tri_A, tri_B, tri_C, tri_ab, tri_ac, tri_n;
        tri_A << pts_list->points[tri_A_idx].x, pts_list->points[tri_A_idx].y,
                pts_list->points[tri_A_idx].z;
        tri_B << pts_list->points[tri_B_idx].x, pts_list->points[tri_B_idx].y,
                pts_list->points[tri_B_idx].z;
        tri_C << pts_list->points[tri_C_idx].x, pts_list->points[tri_C_idx].y,
                pts_list->points[tri_C_idx].z;
        tri_ab = tri_B - tri_A;
        tri_ac = tri_C - tri_A;
        tri_n = tri_ab.cross(tri_ac);

        // 添加三角形
        int A_idx = hole_pts_list[0];
        int B_idx = hole_pts_list[1];
        int C_idx = hole_pts_list[2];
        Eigen::Vector3d A, B, C, ab, ac, n;
        A << pts_list->points[A_idx].x, pts_list->points[A_idx].y, pts_list->points[A_idx].z;
        B << pts_list->points[B_idx].x, pts_list->points[B_idx].y, pts_list->points[B_idx].z;
        C << pts_list->points[C_idx].x, pts_list->points[C_idx].y, pts_list->points[C_idx].z;
        ab = B - A;
        ac = C - A;
        n = ab.cross(ac);

        bool if_acw = true;
        if (n.dot(tri_n) < 0)
        {
            if_acw = false;
        }

        // 判断法向是否符合要求  TODO
        if (if_acw)
        {
            // 添加三角形
            std::array<int, 3> temp_add_t;
            temp_add_t[0] = A_idx;
            temp_add_t[1] = B_idx;
            temp_add_t[2] = C_idx;
            added_tri_list.push_back(faces_list_all.size());
            faces_list_all.push_back(temp_add_t);
        }
        else
        {
            // 添加三角形
            std::array<int, 3> temp_add_t;
            temp_add_t[0] = A_idx;
            temp_add_t[1] = C_idx;
            temp_add_t[2] = B_idx;
            added_tri_list.push_back(faces_list_all.size());
            faces_list_all.push_back(temp_add_t);
        }
    }
    // 如果只有四个点
    else if (hole_pts_list.size() == 4)
    {
        // 确定一个较小角
        Eigen::Vector3d v1, v2, v3, v4, v1_4, v1_2, v2_1, v2_3;
        v1 << pts_list->points[hole_pts_list[0]].x, pts_list->points[hole_pts_list[0]].y,
                pts_list->points[hole_pts_list[0]].z;
        v2 << pts_list->points[hole_pts_list[1]].x, pts_list->points[hole_pts_list[1]].y,
                pts_list->points[hole_pts_list[1]].z;
        v3 << pts_list->points[hole_pts_list[2]].x, pts_list->points[hole_pts_list[2]].y,
                pts_list->points[hole_pts_list[2]].z;
        v4 << pts_list->points[hole_pts_list[3]].x, pts_list->points[hole_pts_list[3]].y,
                pts_list->points[hole_pts_list[3]].z;
        v1_4 = v4 - v1;
        v1_2 = v2 - v1;
        v2_1 = v1 - v2;
        v2_3 = v3 - v2;
        double cos_ang_v1 = v1_4.dot(v1_2) / v1_4.norm() / v1_2.norm();
        double cos_ang_v2 = v2_1.dot(v2_3) / v2_1.norm() / v2_3.norm();
        int A_idx, B_idx, C_idx, D_idx;
        if (cos_ang_v1 > cos_ang_v2)
        {
            A_idx = hole_pts_list[3];
            B_idx = hole_pts_list[0];
            C_idx = hole_pts_list[1];
            D_idx = hole_pts_list[2];
        }
        else
        {
            A_idx = hole_pts_list[0];
            B_idx = hole_pts_list[1];
            C_idx = hole_pts_list[2];
            D_idx = hole_pts_list[3];
        }

        // 根据邻接三角形确定法向
        int adj_tri_idx = -1;
        if (adj_tri_idx == -1)
        {
            adj_tri_idx = pt2face_idx[hole_pts_list[1]][0];
        }

        assert(adj_tri_idx != -1);
        int tri_A_idx = faces_list_all[adj_tri_idx][0];
        int tri_B_idx = faces_list_all[adj_tri_idx][1];
        int tri_C_idx = faces_list_all[adj_tri_idx][2];
        Eigen::Vector3d tri_A, tri_B, tri_C, tri_ab, tri_ac, tri_n;
        tri_A << pts_list->points[tri_A_idx].x, pts_list->points[tri_A_idx].y,
                pts_list->points[tri_A_idx].z;
        tri_B << pts_list->points[tri_B_idx].x, pts_list->points[tri_B_idx].y,
                pts_list->points[tri_B_idx].z;
        tri_C << pts_list->points[tri_C_idx].x, pts_list->points[tri_C_idx].y,
                pts_list->points[tri_C_idx].z;
        tri_ab = tri_B - tri_A;
        tri_ac = tri_C - tri_A;
        tri_n = tri_ab.cross(tri_ac);

        // 添加三角形
        Eigen::Vector3d A, B, C, D, ab, ac, ad, n;
        A << pts_list->points[A_idx].x, pts_list->points[A_idx].y, pts_list->points[A_idx].z;
        B << pts_list->points[B_idx].x, pts_list->points[B_idx].y, pts_list->points[B_idx].z;
        C << pts_list->points[C_idx].x, pts_list->points[C_idx].y, pts_list->points[C_idx].z;
        D << pts_list->points[D_idx].x, pts_list->points[D_idx].y, pts_list->points[D_idx].z;
        ab = B - A;
        ac = C - A;
        ad = D - A;

        // 第一个三角形
        n = ab.cross(ac);
        bool if_acw = true;
        if (n.dot(tri_n) < 0)
        {
            if_acw = false;
        }

        // 判断法向是否符合要求  TODO
        if (if_acw)
        {
            // 添加三角形
            std::array<int, 3> temp_add_t;
            temp_add_t[0] = A_idx;
            temp_add_t[1] = B_idx;
            temp_add_t[2] = C_idx;
            added_tri_list.push_back(faces_list_all.size());
            faces_list_all.push_back(temp_add_t);
        }
        else
        {
            // 添加三角形
            std::array<int, 3> temp_add_t;
            temp_add_t[0] = A_idx;
            temp_add_t[1] = C_idx;
            temp_add_t[2] = B_idx;
            added_tri_list.push_back(faces_list_all.size());
            faces_list_all.push_back(temp_add_t);
        }

        // 第二个三角形
        n = ad.cross(ac);
        if_acw = true;
        if (n.dot(tri_n) < 0)
        {
            if_acw = false;
        }

        // 判断法向是否符合要求  TODO
        if (if_acw)
        {
            // 添加三角形
            std::array<int, 3> temp_add_t;
            temp_add_t[0] = A_idx;
            temp_add_t[1] = D_idx;
            temp_add_t[2] = C_idx;
            added_tri_list.push_back(faces_list_all.size());
            faces_list_all.push_back(temp_add_t);
        }
        else
        {
            // 添加三角形
            std::array<int, 3> temp_add_t;
            temp_add_t[0] = A_idx;
            temp_add_t[1] = C_idx;
            temp_add_t[2] = D_idx;
            added_tri_list.push_back(faces_list_all.size());
            faces_list_all.push_back(temp_add_t);
        }
    }
    // 多于四个点
    else
    {
        std::vector<int> new_hole_pts_list;
        std::vector<double> pts_cos_ang_list;
        for (int i = 0; i < hole_pts_list.size(); i++)
        {
            // 判断是否有大于180度的顶点
            int cur_idx = hole_pts_list[i];
            int v1_idx = hole_pts_list[(i-1+hole_pts_list.size()) % hole_pts_list.size()];
            int v2_idx = hole_pts_list[(i+1+hole_pts_list.size()) % hole_pts_list.size()];
            int v3_idx = -1;
            for (auto it : pt2face_idx[cur_idx])
            {
                for (int l = 0; l < 3; l++)
                {
                    if (faces_list_all[it][l] == v1_idx)
                    {
                        for (int ll = 0; ll < 3; ll++)
                        {
                            if (faces_list_all[it][ll] != cur_idx && faces_list_all[it][ll] != v1_idx)
                            {
                                v3_idx = faces_list_all[it][ll];
                                break;
                            }
                        }
                        break;
                    }
                }
                if (v3_idx != -1)
                {
                    break;
                }
            }
//            assert(v3_idx != -1);

            Eigen::Vector3d cur_pt, v1, v2, v3, cur_v1, cur_v2, cur_v3;
            cur_pt << pts_list->points[cur_idx].x, pts_list->points[cur_idx].y,
                    pts_list->points[cur_idx].z;
            v1 << pts_list->points[v1_idx].x, pts_list->points[v1_idx].y,
                    pts_list->points[v1_idx].z;
            v2 << pts_list->points[v2_idx].x, pts_list->points[v2_idx].y,
                    pts_list->points[v2_idx].z;
            cur_v1 = v1 - cur_pt;
            cur_v2 = v2 - cur_pt;

            if (v3_idx == -1)
            {
                // 递归中间断点 视为小于180度 TODO
                new_hole_pts_list.push_back(hole_pts_list[i]);
                double cur_cos = cur_v1.dot(cur_v2) / cur_v1.norm() / cur_v2.norm();
                pts_cos_ang_list.push_back(cur_cos);
                continue;
            }


            v3 << pts_list->points[v3_idx].x, pts_list->points[v3_idx].y,
                    pts_list->points[v3_idx].z;
            cur_v3 = v3 - cur_pt;
            Eigen::Vector3d cv1_cross_cv2 = cur_v1.cross(cur_v2);
            Eigen::Vector3d cv1_cross_cv3 = cur_v1.cross(cur_v3);

            // 当前顶点夹角大于180度
            if (cv1_cross_cv2.dot(cv1_cross_cv3) > 0)
            {
                // 搜索其最近顶点，除去自身和相邻顶点外的最近
                Eigen::Vector3d temp_vv;
                int temp_vv_idx = hole_pts_list[(i-2+hole_pts_list.size()) % hole_pts_list.size()];
                temp_vv << pts_list->points[temp_vv_idx].x, pts_list->points[temp_vv_idx].y,
                        pts_list->points[temp_vv_idx].z;
                double min_dis = (cur_pt - temp_vv).norm();
                int min_dis_v_idx = (i-2+hole_pts_list.size()) % hole_pts_list.size();
                for (int j = 0; j < hole_pts_list.size(); j++)
                {
                    if (j == i) { continue; }
                    if (j == (i-1+hole_pts_list.size()) % hole_pts_list.size()) { continue; }
                    if (j == (i+1+hole_pts_list.size()) % hole_pts_list.size()) { continue; }
                    if (j == (i-2+hole_pts_list.size()) % hole_pts_list.size()) { continue; }

                    Eigen::Vector3d temp_v;
                    temp_v << pts_list->points[hole_pts_list[j]].x, pts_list->points[hole_pts_list[j]].y,
                            pts_list->points[hole_pts_list[j]].z;
                    double v_dis = (cur_pt - temp_v).norm();
                    if (v_dis < min_dis)
                    {
                        min_dis = v_dis;
                        min_dis_v_idx = j;
                    }
                }
                // 相连将一个洞分为两个，分别处理  TODO 可能一次处理后角还是大于180度
                std::vector<int> one_hole_pts_list, another_hole_pts_list;
                if (min_dis_v_idx < i)
                {
                    // 一个洞的所有点索引
                    for (int j = 0; j <= min_dis_v_idx; j++)
                    {
                        one_hole_pts_list.push_back(hole_pts_list[j]);
                    }
                    for (int j = i; j < hole_pts_list.size(); j++)
                    {
                        one_hole_pts_list.push_back(hole_pts_list[j]);
                    }

                    // 另一个洞的所有点索引
                    for (int j = min_dis_v_idx; j <= i; j++)
                    {
                        another_hole_pts_list.push_back(hole_pts_list[j]);
                    }
                }
                else
                {
                    // 一个洞的所有点索引
                    for (int j = 0; j <= i; j++)
                    {
                        one_hole_pts_list.push_back(hole_pts_list[j]);
                    }
                    for (int j = min_dis_v_idx; j < hole_pts_list.size(); j++)
                    {
                        one_hole_pts_list.push_back(hole_pts_list[j]);
                    }

                    // 另一个洞的所有点索引
                    for (int j = i; j <= min_dis_v_idx; j++)
                    {
                        another_hole_pts_list.push_back(hole_pts_list[j]);
                    }
                }

                // 递归处理两个洞
                fill_hole(one_hole_pts_list);
                fill_hole(another_hole_pts_list);
                new_hole_pts_list.clear();
                break;

            }
            // 当前顶点夹角小于180度
            else
            {
                new_hole_pts_list.push_back(hole_pts_list[i]);
                double cur_cos = cur_v1.dot(cur_v2) / cur_v1.norm() / cur_v2.norm();
                pts_cos_ang_list.push_back(cur_cos);
            }

        }

        // 计算相邻三角形的法向
        Eigen::Vector3d adj_tri_n;
        if (new_hole_pts_list.size() > 2)
        {
            int cur_idx, v1_idx, f_idx;
            cur_idx = new_hole_pts_list[1];
            v1_idx = new_hole_pts_list[2];
            f_idx = -1;

            if (f_idx == -1)
            {
                f_idx = pt2face_idx[cur_idx][0];
            }

            assert(f_idx != -1);

            int A_idx = faces_list_all[f_idx][0];
            int B_idx = faces_list_all[f_idx][1];
            int C_idx = faces_list_all[f_idx][2];
            Eigen::Vector3d A, B, C, ab, ac;
            A << pts_list->points[A_idx].x, pts_list->points[A_idx].y,
                    pts_list->points[A_idx].z;
            B << pts_list->points[B_idx].x, pts_list->points[B_idx].y,
                    pts_list->points[B_idx].z;
            C << pts_list->points[C_idx].x, pts_list->points[C_idx].y,
                    pts_list->points[C_idx].z;
            ab = B - A;
            ac = C - A;
            adj_tri_n = ab.cross(ac);
        }

        // 对内角全小于180度的多边形内部构建三角剖分
        while (new_hole_pts_list.size() > 2)
        {
            // 内角全小于180度，按角度排序
            std::unordered_map<int, double> ptsidx_to_cos;
            for (int i = 0; i < new_hole_pts_list.size(); i++)
            {
                ptsidx_to_cos[i] = pts_cos_ang_list[i];
            }

            std::vector<int> new_hole_pts_idx_sorted;
            std::vector<std::pair<int, double>> indexVec(ptsidx_to_cos.begin(), ptsidx_to_cos.end());
            std::sort(indexVec.begin(), indexVec.end(), compare_1);  // cos从大到小，角度从小到大
            for (const auto& pair : indexVec)
            {
                new_hole_pts_idx_sorted.push_back(pair.first);
            }

            // 找最小角+两个相邻点组成三角形
            int min_cos_pt_idx, cur_idx, v1_idx, v2_idx;
            Eigen::Vector3d cur_pt, v1, v2, cur_v1, cur_v2, n;

            min_cos_pt_idx = new_hole_pts_idx_sorted[0];
            cur_idx = new_hole_pts_list[min_cos_pt_idx];
            v1_idx = new_hole_pts_list[(min_cos_pt_idx-1+new_hole_pts_list.size()) % new_hole_pts_list.size()];
            v2_idx = new_hole_pts_list[(min_cos_pt_idx+1+new_hole_pts_list.size()) % new_hole_pts_list.size()];

            // 添加三角形
            cur_pt << pts_list->points[cur_idx].x, pts_list->points[cur_idx].y,
                    pts_list->points[cur_idx].z;
            v1 << pts_list->points[v1_idx].x, pts_list->points[v1_idx].y,
                    pts_list->points[v1_idx].z;
            v2 << pts_list->points[v2_idx].x, pts_list->points[v2_idx].y,
                    pts_list->points[v2_idx].z;
            cur_v1 = v1 - cur_pt;
            cur_v2 = v2 - cur_pt;
            n = cur_v1.cross(cur_v2);

            bool if_acw = true;
            if (n.dot(adj_tri_n) < 0)
            {
                if_acw = false;
            }

            // 判断法向是否符合要求  TODO
            if (if_acw)
            {
                // 添加三角形
                std::array<int, 3> temp_add_t;
                temp_add_t[0] = cur_idx;
                temp_add_t[1] = v1_idx;
                temp_add_t[2] = v2_idx;
                added_tri_list.push_back(faces_list_all.size());
                faces_list_all.push_back(temp_add_t);
            }
            else
            {
                // 添加三角形
                std::array<int, 3> temp_add_t;
                temp_add_t[0] = cur_idx;
                temp_add_t[1] = v2_idx;
                temp_add_t[2] = v1_idx;
                added_tri_list.push_back(faces_list_all.size());
                faces_list_all.push_back(temp_add_t);
            }

            if (new_hole_pts_list.size() > 3)
            {
                // 更新角度列表，删除 + 更新
                // 第一个角cos值
                cur_idx = new_hole_pts_list[(min_cos_pt_idx-1+new_hole_pts_list.size()) % new_hole_pts_list.size()];
                v1_idx = new_hole_pts_list[(min_cos_pt_idx-2+new_hole_pts_list.size()) % new_hole_pts_list.size()];
                v2_idx = new_hole_pts_list[(min_cos_pt_idx+1+new_hole_pts_list.size()) % new_hole_pts_list.size()];

                cur_pt << pts_list->points[cur_idx].x, pts_list->points[cur_idx].y,
                        pts_list->points[cur_idx].z;
                v1 << pts_list->points[v1_idx].x, pts_list->points[v1_idx].y,
                        pts_list->points[v1_idx].z;
                v2 << pts_list->points[v2_idx].x, pts_list->points[v2_idx].y,
                        pts_list->points[v2_idx].z;
                cur_v1 = v1 - cur_pt;
                cur_v2 = v2 - cur_pt;
                pts_cos_ang_list[(min_cos_pt_idx-1+new_hole_pts_list.size()) % new_hole_pts_list.size()] =
                        cur_v1.dot(cur_v2) / cur_v1.norm() / cur_v2.norm();
                // 第二个角cos值
                cur_idx = new_hole_pts_list[(min_cos_pt_idx+1+new_hole_pts_list.size()) % new_hole_pts_list.size()];
                v1_idx = new_hole_pts_list[(min_cos_pt_idx-1+new_hole_pts_list.size()) % new_hole_pts_list.size()];
                v2_idx = new_hole_pts_list[(min_cos_pt_idx+2+new_hole_pts_list.size()) % new_hole_pts_list.size()];

                cur_pt << pts_list->points[cur_idx].x, pts_list->points[cur_idx].y,
                        pts_list->points[cur_idx].z;
                v1 << pts_list->points[v1_idx].x, pts_list->points[v1_idx].y,
                        pts_list->points[v1_idx].z;
                v2 << pts_list->points[v2_idx].x, pts_list->points[v2_idx].y,
                        pts_list->points[v2_idx].z;
                cur_v1 = v1 - cur_pt;
                cur_v2 = v2 - cur_pt;
                pts_cos_ang_list[(min_cos_pt_idx+1+new_hole_pts_list.size()) % new_hole_pts_list.size()] =
                        cur_v1.dot(cur_v2) / cur_v1.norm() / cur_v2.norm();
                // 第三个角cos值
                pts_cos_ang_list.erase(pts_cos_ang_list.begin() + min_cos_pt_idx);
            }
            // 更新边缘点列表，删除
            new_hole_pts_list.erase(new_hole_pts_list.begin() + min_cos_pt_idx);
        }
    }

    // 更新点-三角形索引列表
    for (int j = 0; j < added_tri_list.size(); j++)
    {
        for (int l = 0; l < 3; l++)
        {
            int temp_id = faces_list_all[added_tri_list[j]][l];
            pt2face_idx[temp_id].push_back(added_tri_list[j]);
        }
    }


}

void NonPlaneMesh::mesh_hole_process(FILE *f)
{
    // 初始化所有点的状态为 未处理
    int pt_list_size = tri_edge_v_list_all.size();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, pt_list_size), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t j = r.begin(); j != r.end(); ++j)
        {
            int i = tri_edge_v_list_all[j];
            hole_v_process[i] = 0;

        }
    });

    // 并行检测hole
    std::vector<int> delete_edge_idx_list;
    std::vector<std::vector<int>> hole_list;
    std::vector<int> hole_v_idx_list;
    std::mutex mtx_delete_edge_idx_list, mtx_edge_v_state, mtx_hole;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, pt_list_size), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t j = r.begin(); j != r.end(); ++j)
        {
            int v_idx = tri_edge_v_list_all[j];
            // 已处理过
            if (hole_v_process[v_idx] == 1)
            {
                continue;
            }

            // 检查状态
            if (pts_state[v_idx] != 2)
            {
                mtx_delete_edge_idx_list.lock();
                delete_edge_idx_list.push_back(j);
                mtx_delete_edge_idx_list.unlock();
                continue;
            }
            // 只处理 多次操作仍然是边缘点的 点
            if (pts_process_num[v_idx] < 10)
            {
                continue;
            }

            // 检测该点是否处于闭合空洞
            std::vector<int> edge_v;
            std::vector<int> adj_edge_vt;
            find_edge_vertex(v_idx, adj_edge_vt);

            if (adj_edge_vt.size() == 0)
            {
                pts_state[v_idx] = 1;
                mtx_delete_edge_idx_list.lock();
                delete_edge_idx_list.push_back(j);
                mtx_delete_edge_idx_list.unlock();
                continue;
            }

            bool if_hole = if_lay_in_hole(v_idx, adj_edge_vt, edge_v);

            mtx_edge_v_state.lock();
            for (auto it : edge_v)
            {
                hole_v_process[it] = 1;
            }
            mtx_edge_v_state.unlock();

            if (if_hole)
            {
                mtx_hole.lock();
                bool if_new = false;
                for (auto it : edge_v)
                {
                    auto temp_it = std::find(hole_v_idx_list.begin(), hole_v_idx_list.end(), it);
                    if (temp_it == hole_v_idx_list.end())
                    {
                        if_new = true;
                        break;
                    }
                }

                if (if_new)
                {
                    for (auto it : edge_v)
                    {
                        hole_v_idx_list.push_back(it);
                    }
                    hole_list.push_back(edge_v);
                }
                mtx_hole.unlock();
            }
            else
            {
                hole_process_num[v_idx]++;
                if (hole_process_num[v_idx] > 1000)
                {
                    mtx_delete_edge_idx_list.lock();
                    delete_edge_idx_list.push_back(j);
                    mtx_delete_edge_idx_list.unlock();
                }
            }
        }
    });

    // 补洞
    for(auto edge_v : hole_list)
    {
        // 添加三角形补洞
        fill_hole(edge_v);

        // 更新添加三角形相关联点的状态
        for (auto it : edge_v)
        {
            if (!if_edge_vertex(it))
            {
                if (pts_state[it] == 2)
                {
                    pts_state[it] = 1;
                }

                if (pts_state[it] == 5)
                {
                    pts_state[it] = 6;
                }
            }
        }
    }

    // 删去非边缘点
    std::sort(delete_edge_idx_list.begin(), delete_edge_idx_list.end());
    size_t i = 0;
    for (auto it : delete_edge_idx_list)
    {
        size_t it_it = it;
        tri_edge_v_list_all.erase(tri_edge_v_list_all.begin() + it_it - i);
        i++;
    }
}

// 思想参考IDTMM，但是由于其没有开源，实际上存在差异
void NonPlaneMesh::mesh_update(const std::vector<std::array<long, 3>>& voxel_to_update,
                               Eigen::Vector3d& body_pos/*, FILE *f*/)
{
    // 依次处理每个体素
    for (const auto& voxel_idx : voxel_to_update)
    {
        long voxel_x = voxel_idx[0];
        long voxel_y = voxel_idx[1];
        long voxel_z = voxel_idx[2];

        // 第一次遍历，处理自由点，连接自由点与三角形顶点形成新的mesh
        for (int i = 0; i < voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->all_pts_idx.size(); i++)
        {
            int current_idx = voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->all_pts_idx[i];
            // 若不是孤立点则跳过
            if (pts_state[current_idx] != 0)
            {
                continue;
            }

            std::vector<int> nearest_vertexs;
            std::vector<double> nearest_vertexs_dis;
            search_nearest_vertex_under_distance_r(current_idx, nearest_vertexs, nearest_vertexs_dis);

            if (nearest_vertexs.size() == 0)
            {
                std::vector<int> nearest_pts;
                search_nearest_point_under_distance_r(current_idx, nearest_pts);
                /*** 情况 1 ***/
                if (nearest_pts.size() < N_pts)
                {
                    continue;
                }
                /*** 情况 2 ***/
                else
                {
                    // 按协方差矩阵计算的方向投影，在二维新建delaunay三角剖分
                    nearest_pts.push_back(current_idx);
                    // 确定表面的大致法向
                    Eigen::Vector3d avg_normal;
                    avg_normal(0) = body_pos(0) - pts_list->points[current_idx].x;
                    avg_normal(1) = body_pos(1) - pts_list->points[current_idx].y;
                    avg_normal(2) = body_pos(2) - pts_list->points[current_idx].z;

                    // 拟合两个三角形
                    std::vector<std::array<int, 3>> tri_pt_indices;
                    fit_2_tri(nearest_pts, avg_normal, tri_pt_indices);

                    // 更新三角形列表、点-三角形索引
                    for (int j = 0; j < tri_pt_indices.size(); j++)
                    {
                        mtx_faces_list_all.lock();
                        int cur_tri_idx = faces_list_all.size();
                        faces_list_all.push_back(tri_pt_indices[j]);
                        mtx_faces_list_all.unlock();

                        mtx_pt2face_idx.lock();
                        for (int k = 0; k < 3; k++)
                        {
                            pt2face_idx[tri_pt_indices[j][k]].push_back(cur_tri_idx);
                        }
                        mtx_pt2face_idx.unlock();
                    }
                    // 更新点的状态：初始点为删除状态
                    mtx_pts_state.lock();
                    for (const auto& idx : nearest_pts)
                    {
                        pts_state[idx] = 3;
                    }
                    mtx_pts_state.unlock();
                }
            }
            else
            {
                // 判断沿法向投影在三角形内部还是外部
                std::vector<int> tri_tested_idx;
                bool if_inner = false;
                int projected_tri_idx;
                for (int j = 0; j < nearest_vertexs.size(); j++)
                {
                    int vertex_idx = nearest_vertexs[j];
                    for (int k = 0; k < pt2face_idx[vertex_idx].size(); k++)
                    {
                        int face_idx = pt2face_idx[vertex_idx][k];
                        auto it = std::find(tri_tested_idx.begin(), tri_tested_idx.end(), face_idx);
                        if (it != tri_tested_idx.end())
                        {
                            // 该三角形已经检查过了
                            continue;
                        }
                        tri_tested_idx.push_back(face_idx);

                        int A_idx = faces_list_all[face_idx][0];
                        int B_idx = faces_list_all[face_idx][1];
                        int C_idx = faces_list_all[face_idx][2];
                        Eigen::Vector3d Q, A, B, C;
                        Q << pts_list->points[current_idx].x, pts_list->points[current_idx].y,
                                pts_list->points[current_idx].z;
                        A << pts_list->points[A_idx].x, pts_list->points[A_idx].y,
                                pts_list->points[A_idx].z;
                        B << pts_list->points[B_idx].x, pts_list->points[B_idx].y,
                                pts_list->points[B_idx].z;
                        C << pts_list->points[C_idx].x, pts_list->points[C_idx].y,
                                pts_list->points[C_idx].z;
                        // 沿三角形法向投影在三角形内部
                        if ( if_projected_inner(Q, A, B, C, body_pos) )
                        {
                            if_inner = true;
                            projected_tri_idx = face_idx;
                            break;                             // TODO 是否只投影在一个三角形内部？多个判断离哪个最近？
                        }
                    }
                    if (if_inner) {break;}
                }
                /*** 情况 3 ***/
                if (if_inner)
                {
                    // 卡尔曼滤波器优化点的位置 TODO

                    // 更新当前点的状态
                    mtx_pts_state.lock();
                    pts_state[current_idx] = 3;
                    mtx_pts_state.unlock();
                }
                /*** 情况 4 ***/
                else
                {
                    // 保持任意两个点之间的最小距离
                    if (nearest_vertexs_dis[0] < minimum_pt_dis)
                    {
                        mtx_pts_state.lock();
                        pts_state[current_idx] = 3;
                        mtx_pts_state.unlock();
                    }
                    else
                    {
                        // 确定距离r内顶点关联的三角形列表
                        std::set<int> voxel_faces;
                        for (auto it_vertex : nearest_vertexs)  // TODO 如何利用上顶点的已完成/未完成状态
                        {
                            for (auto it_face : pt2face_idx[it_vertex])
                            {
                                voxel_faces.insert(it_face);
                            }
                        }

                        // 建立边-三角形索引
                        std::unordered_map<std::string, int> edge_tri_idx;
                        std::vector<std::vector<std::array<int, 3>>> edge_tri_list;
                        std::vector<std::vector<int>> edge_tri_id_list;
                        for (auto it : voxel_faces)
                        {
                            for (int j = 0; j < 3; j++)
                            {
                                int a, b, c;
                                if (j ==0 )
                                {
                                    a = faces_list_all[it][0];
                                    b = faces_list_all[it][1];
                                    c = faces_list_all[it][2];
                                }
                                else if (j == 1)
                                {
                                    a = faces_list_all[it][1];
                                    b = faces_list_all[it][2];
                                    c = faces_list_all[it][0];
                                }
                                else
                                {
                                    a = faces_list_all[it][2];
                                    b = faces_list_all[it][0];
                                    c = faces_list_all[it][1];
                                }
                                if (a > b)
                                {
                                    int temp_int = a;
                                    a = b;
                                    b = temp_int;
                                }
                                std::string edge_key = std::to_string(a) + ' ' + std::to_string(b);
                                if (edge_tri_idx.find(edge_key) == edge_tri_idx.end())
                                {
                                    edge_tri_idx[edge_key] = edge_tri_list.size();
                                    std::vector<std::array<int, 3>> temp_a;
                                    std::array<int, 3> temp_b;
                                    temp_b[0] = a;
                                    temp_b[1] = b;
                                    temp_b[2] = c;
                                    temp_a.push_back(temp_b);
                                    edge_tri_list.push_back(temp_a);

                                    std::vector<int> temp_c;
                                    temp_c.push_back(it);
                                    edge_tri_id_list.push_back(temp_c);
                                }
                                else
                                {
                                    int temp_idx = edge_tri_idx[edge_key];
                                    std::array<int, 3> temp_b;
                                    temp_b[0] = a;
                                    temp_b[1] = b;
                                    temp_b[2] = c;
                                    edge_tri_list[temp_idx].push_back(temp_b);
                                }
                            }
                        }
                        // 确定边缘边，格式 A-B-C
                        std::vector<std::array<int, 3>> edge_edge_inner;
                        std::vector<int> edge_edge_inner_tri_id;
                        for (int j = 0; j < edge_tri_list.size(); j++)
                        {
                            if (edge_tri_list[j].size() == 1)
                            {
                                // 进一步判断是否真正边缘边，由于后面会判断是否异侧，这里不处理也可以！
                                int edge_v1 = edge_tri_list[j][0][0];
                                int edge_v2 = edge_tri_list[j][0][1];
                                int common_tri_num = 0;
                                for (auto it : pt2face_idx[edge_v1])
                                {
                                    for (int l = 0; l < 3; l++)
                                    {
                                        if (faces_list_all[it][l] == edge_v2)
                                        {
                                            common_tri_num++;
                                        }
                                    }
                                }
                                if (common_tri_num == 1)
                                {
                                    edge_edge_inner.push_back(edge_tri_list[j][0]);
                                    edge_edge_inner_tri_id.push_back(edge_tri_id_list[j][0]);
                                }
                            }
                        }

                        // 依次处理每个边缘边
                        Eigen::Vector3d cur_pt_eigen;
                        cur_pt_eigen << pts_list->points[current_idx].x, pts_list->points[current_idx].y,
                                pts_list->points[current_idx].z;
                        Eigen::Vector3d lidar_beam = body_pos - cur_pt_eigen;
                        std::vector<int> added_tri_list, deleted_tri_list;
                        for (int j = 0; j < edge_edge_inner.size(); j++)
                        {
                            // 判断相邻三角形法向与lidar-point方向  TODO 改为判断两个三角形的法向夹角
                            int tri_A_idx = faces_list_all[edge_edge_inner_tri_id[j]][0];
                            int tri_B_idx = faces_list_all[edge_edge_inner_tri_id[j]][1];
                            int tri_C_idx = faces_list_all[edge_edge_inner_tri_id[j]][2];
                            Eigen::Vector3d tri_A, tri_B, tri_C;
                            tri_A << pts_list->points[tri_A_idx].x, pts_list->points[tri_A_idx].y,
                                    pts_list->points[tri_A_idx].z;
                            tri_B << pts_list->points[tri_B_idx].x, pts_list->points[tri_B_idx].y,
                                    pts_list->points[tri_B_idx].z;
                            tri_C << pts_list->points[tri_C_idx].x, pts_list->points[tri_C_idx].y,
                                    pts_list->points[tri_C_idx].z;
                            Eigen::Vector3d tri_ab = tri_B - tri_A;
                            Eigen::Vector3d tri_ac = tri_C - tri_A;
                            Eigen::Vector3d tri_n = tri_ab.cross(tri_ac);

                            // 新三角形的相关点
                            int A_idx = edge_edge_inner[j][0];
                            int B_idx = edge_edge_inner[j][1];
                            int C_idx = edge_edge_inner[j][2];
                            Eigen::Vector3d A, B, C;
                            A << pts_list->points[A_idx].x, pts_list->points[A_idx].y,
                                    pts_list->points[A_idx].z;
                            B << pts_list->points[B_idx].x, pts_list->points[B_idx].y,
                                    pts_list->points[B_idx].z;
                            C << pts_list->points[C_idx].x, pts_list->points[C_idx].y,
                                    pts_list->points[C_idx].z;
                            Eigen::Vector3d cur_A = A - cur_pt_eigen;
                            Eigen::Vector3d cur_B = B - cur_pt_eigen;
                            Eigen::Vector3d A_B = B - A;
                            Eigen::Vector3d A_cur = cur_pt_eigen - A;
                            Eigen::Vector3d B_cur = cur_pt_eigen - B;
                            Eigen::Vector3d B_A = A - B;

                            // 判断边长是否符合要求
                            if ( cur_A.norm() > r_dis * N_r_dis) { continue; }
                            if ( cur_B.norm() > r_dis * N_r_dis) { continue; }

                            // 判断三角形法向投影是否在边的两侧
                            if (!if_different_side(A, B, C, cur_pt_eigen))
                            {
                                continue;
                            }

                            // 若夹角过大则，删除再添加
                            double cos_cur_AB = cur_A.dot(cur_B) / cur_A.norm() / cur_B.norm();
                            if (cos_cur_AB < -0.5)  // 140-180  -0.766
                            {
                                // 删除三角形
                                mtx_faces_to_delete.lock();
                                faces_to_delete.push_back(edge_edge_inner_tri_id[j]);
                                mtx_faces_to_delete.unlock();
                                deleted_tri_list.push_back(edge_edge_inner_tri_id[j]);
                                // 确定第三个点的索引
                                int v3_idx;
                                if (tri_A_idx != A_idx && tri_A_idx != B_idx)
                                {
                                    v3_idx = tri_A_idx;
                                }
                                else if (tri_B_idx != A_idx && tri_B_idx != B_idx)
                                {
                                    v3_idx = tri_B_idx;
                                }
                                else
                                {
                                    v3_idx = tri_C_idx;
                                }
                                // 添加三角形
                                Eigen::Vector3d v3, cur_v3, n;
                                v3 << pts_list->points[v3_idx].x, pts_list->points[v3_idx].y,
                                        pts_list->points[v3_idx].z;
                                cur_v3 = v3 - cur_pt_eigen;
                                // 第一个三角形
                                n = cur_v3.cross(cur_A);
                                bool if_acw_nc = true;
                                if (n.dot(tri_n) < 0)
                                {
                                    if_acw_nc = false;
                                }
                                if (if_acw_nc)
                                {
                                    std::array<int, 3> temp_add_t;
                                    temp_add_t[0] = current_idx;
                                    temp_add_t[1] = v3_idx;
                                    temp_add_t[2] = A_idx;

                                    mtx_faces_list_all.lock();
                                    added_tri_list.push_back(faces_list_all.size());
                                    faces_list_all.push_back(temp_add_t);
                                    mtx_faces_list_all.unlock();
                                }
                                else
                                {
                                    std::array<int, 3> temp_add_t;
                                    temp_add_t[0] = current_idx;
                                    temp_add_t[1] = A_idx;
                                    temp_add_t[2] = v3_idx;

                                    mtx_faces_list_all.lock();
                                    added_tri_list.push_back(faces_list_all.size());
                                    faces_list_all.push_back(temp_add_t);
                                    mtx_faces_list_all.unlock();
                                }
                                // 第二个三角形
                                n = cur_v3.cross(cur_B);
                                if_acw_nc = true;
                                if (n.dot(tri_n) < 0)
                                {
                                    if_acw_nc = false;
                                }
                                if (if_acw_nc)
                                {
                                    std::array<int, 3> temp_add_t;
                                    temp_add_t[0] = current_idx;
                                    temp_add_t[1] = v3_idx;
                                    temp_add_t[2] = B_idx;

                                    mtx_faces_list_all.lock();
                                    added_tri_list.push_back(faces_list_all.size());
                                    faces_list_all.push_back(temp_add_t);
                                    mtx_faces_list_all.unlock();
                                }
                                else
                                {
                                    std::array<int, 3> temp_add_t;
                                    temp_add_t[0] = current_idx;
                                    temp_add_t[1] = B_idx;
                                    temp_add_t[2] = v3_idx;

                                    mtx_faces_list_all.lock();
                                    added_tri_list.push_back(faces_list_all.size());
                                    faces_list_all.push_back(temp_add_t);
                                    mtx_faces_list_all.unlock();
                                }
                                continue;
                            }

                            // 判断角度是否符合要求
                            if (cos_cur_AB > ang_thre) { continue; }
                            if (A_B.dot(A_cur) / A_B.norm() / A_cur.norm() > ang_thre) { continue; }
                            if (B_A.dot(B_cur) / B_A.norm() / B_cur.norm() > ang_thre) { continue; }

                            // TODO 判断可见性

                            // 根据lidar射线/相邻三角形法向确定三个点的逆时针顺序
                            Eigen::Vector3d norm = cur_A.cross(cur_B);
                            bool if_acw = true;
                            if (norm.dot(lidar_beam) < 0)
                            {
                                if_acw = false;
                            }

                            // 判断法向是否符合要求  TODO
                            if (if_acw)
                            {
                                //if (if_reasonable_normal(current_idx, A_idx, B_idx))
                                if (1)
                                {
                                    // 添加三角形
                                    std::array<int, 3> temp_add_t;
                                    temp_add_t[0] = current_idx;
                                    temp_add_t[1] = A_idx;
                                    temp_add_t[2] = B_idx;

                                    mtx_faces_list_all.lock();
                                    added_tri_list.push_back(faces_list_all.size());
                                    faces_list_all.push_back(temp_add_t);
                                    mtx_faces_list_all.unlock();
                                }
                            }
                            else
                            {
                                //if (if_reasonable_normal(current_idx, B_idx, A_idx))
                                if (1)
                                {
                                    // 添加三角形
                                    std::array<int, 3> temp_add_t;
                                    temp_add_t[0] = current_idx;
                                    temp_add_t[1] = B_idx;
                                    temp_add_t[2] = A_idx;

                                    mtx_faces_list_all.lock();
                                    added_tri_list.push_back(faces_list_all.size());
                                    faces_list_all.push_back(temp_add_t);
                                    mtx_faces_list_all.unlock();
                                }
                            }
                        }

                        // 添加不连续开放边之间的三角形
                        std::vector<int> added_tri_list_nc;
                        if (added_tri_list.size() > 1)
                        {
                            // 三角形顶点沿法向投影，逆(顺)时针排序
                            std::vector<int> added_tri_list_sorted;
                            std::unordered_map<int, int> ptidx_to_acidx;

                            std::set<int> candidate_vertex;
                            std::vector<std::vector<int>> tri_v_idx;
                            tri_v_idx.resize(added_tri_list.size());
                            for (int j = 0; j < added_tri_list.size(); j++)
                            {
                                int idx = added_tri_list[j];
                                for (int k = 0; k < 3; k++)
                                {
                                    if (faces_list_all[idx][k] == current_idx) { continue; }
                                    candidate_vertex.insert(faces_list_all[idx][k]);
                                    tri_v_idx[j].push_back(faces_list_all[idx][k]);
                                }
                            }
                            std::vector<int> can_v, can_v_sorted, pts_to_axis;  // 最后一个目前没用到
                            for(auto it = candidate_vertex.begin(); it != candidate_vertex.end(); it++)
                            {
                                can_v.push_back(*it);
                            }
                            anticlockwise_sort(can_v, current_idx, tri_v_idx, can_v_sorted, pts_to_axis);

                            // 三角形排序
                            for (int j = 0; j < can_v_sorted.size(); j++)
                            {
                                ptidx_to_acidx[can_v_sorted[j]] = j;
                            }
                            std::vector<std::pair<int, double>> tri_acidx;
                            for (int j = 0; j < added_tri_list.size(); j++)
                            {
                                std::pair<int, double> temp_pair;
                                temp_pair.first = j;
                                if (ptidx_to_acidx[tri_v_idx[j][0]] == 0 || ptidx_to_acidx[tri_v_idx[j][1]] == 0)
                                {
                                    if (ptidx_to_acidx[tri_v_idx[j][0]] == 1 || ptidx_to_acidx[tri_v_idx[j][1]] == 1)
                                    {
                                        temp_pair.second = 1.0;
                                    }
                                    else
                                    {
                                        temp_pair.second = double(2 * can_v.size() - 1);
                                    }
                                }
                                else
                                {
                                    temp_pair.second = double(ptidx_to_acidx[tri_v_idx[j][0]] + ptidx_to_acidx[tri_v_idx[j][1]]);
                                }
                                tri_acidx.push_back(temp_pair);
                            }
                            std::sort(tri_acidx.begin(), tri_acidx.end(), compare);
                            for (const auto& pair : tri_acidx)
                            {
                                added_tri_list_sorted.push_back(pair.first);
                            }

                            // 判断条件，添加三角形
                            for (int j = 0; j < added_tri_list.size(); j++)
                            {
                                int j_sorted = added_tri_list_sorted[j];
                                int j_next = added_tri_list_sorted[(j+1)%added_tri_list.size()];
                                bool if_continuous = false;
                                int last_A, last_B, next_A, next_B;
                                last_A = tri_v_idx[j_sorted][0];
                                last_B = tri_v_idx[j_sorted][1];
                                next_A = tri_v_idx[j_next][0];
                                next_B = tri_v_idx[j_next][1];
                                if (last_A == next_A || last_A == next_B || last_B == next_A || last_B == next_B)
                                {
                                    if_continuous = true;
                                }

                                if (!if_continuous)
                                {
                                    // 确定三个顶点的索引
                                    int idx_1, idx_2, idx_1_pre;
                                    if (ptidx_to_acidx[last_A] > ptidx_to_acidx[last_B])
                                    {
                                        idx_1 = last_A;
                                        idx_1_pre = last_B;
                                    }
                                    else
                                    {
                                        idx_1 = last_B;
                                        idx_1_pre = last_A;
                                    }
                                    if (ptidx_to_acidx[next_A] == 0)
                                    {
                                        if (ptidx_to_acidx[next_B] == 1)
                                        {
                                            idx_2 = next_A;
                                        }
                                        else
                                        {
                                            idx_2 = next_B;
                                        }
                                    }
                                    else if (ptidx_to_acidx[next_B] == 0)
                                    {
                                        if (ptidx_to_acidx[next_A] == 1)
                                        {
                                            idx_2 = next_B;
                                        }
                                        else
                                        {
                                            idx_2 = next_A;
                                        }
                                    }
                                    else
                                    {
                                        if (ptidx_to_acidx[next_A] > ptidx_to_acidx[next_B])
                                        {
                                            idx_2 = next_B;
                                        }
                                        else
                                        {
                                            idx_2 = next_A;
                                        }
                                    }
                                    // 判断两个顶点的顺序是否符合要求
                                    if (ptidx_to_acidx[idx_2] < ptidx_to_acidx[idx_1] && ptidx_to_acidx[idx_2] != 0)
                                    {
                                        continue;
                                    }

                                    Eigen::Vector3d v_1, v_2, v_1_pre;
                                    v_1 << pts_list->points[idx_1].x, pts_list->points[idx_1].y,
                                            pts_list->points[idx_1].z;
                                    v_2 << pts_list->points[idx_2].x, pts_list->points[idx_2].y,
                                            pts_list->points[idx_2].z;
                                    v_1_pre << pts_list->points[idx_1_pre].x, pts_list->points[idx_1_pre].y,
                                            pts_list->points[idx_1_pre].z;
                                    Eigen::Vector3d cur_v1 = v_1 - cur_pt_eigen;
                                    Eigen::Vector3d cur_v2 = v_2 - cur_pt_eigen;
                                    Eigen::Vector3d v1_v2 = v_2 - v_1;
                                    Eigen::Vector3d v1_cur = cur_pt_eigen - v_1;
                                    Eigen::Vector3d v2_cur = cur_pt_eigen - v_2;
                                    Eigen::Vector3d v2_v1 = v_1 - v_2;
                                    Eigen::Vector3d cur_v1pre = v_1_pre - cur_pt_eigen;

                                    // 判断夹角是否小于180度
                                    Eigen::Vector3d cross_1 = cur_v1pre.cross(cur_v1);
                                    Eigen::Vector3d cross_2 = cur_v1.cross(cur_v2);
                                    if (cross_1.dot(cross_2) < 0)
                                    {
                                        continue;
                                    }

                                    // 判断边长是否符合要求
//                                    if ( cur_v1.norm() > r_dis) { continue; }
//                                    if ( cur_v2.norm() > r_dis) { continue; }
//                                    if ( v1_v2.norm() > r_dis) { continue; }    // TODO 调节

                                    // 判断角度是否符合要求 TODO 10-90度之间
                                    if (cur_v1.dot(cur_v2) / cur_v1.norm() / cur_v2.norm() > ang_thre) { continue; }
                                    if (v1_v2.dot(v1_cur) / v1_v2.norm() / v1_cur.norm() > ang_thre) { continue; }
                                    if (v2_v1.dot(v2_cur) / v2_v1.norm() / v2_cur.norm() > ang_thre) { continue; }

                                    // 判断有无三角形顶点位于三角形内部
                                    bool if_any_v_inner = false;
                                    for (int nvi = 0; nvi < nearest_vertexs.size(); nvi++)
                                    {
                                        if (nearest_vertexs[nvi] == idx_1 || nearest_vertexs[nvi] == idx_2)
                                        {
                                            continue;
                                        }
                                        Eigen::Vector3d temp_v;
                                        temp_v << pts_list->points[nearest_vertexs[nvi]].x, pts_list->points[nearest_vertexs[nvi]].y,
                                                pts_list->points[nearest_vertexs[nvi]].z;
                                        if (if_p_in_tri(temp_v, cur_pt_eigen, v_1, v_2))
                                        {
                                            if_any_v_inner = true;
                                            break;
                                        }
                                    }
                                    if (if_any_v_inner) { continue; }


                                    // 确定顶点逆时针顺序
                                    int tri_A_idx = faces_list_all[added_tri_list[j_sorted]][0];
                                    int tri_B_idx = faces_list_all[added_tri_list[j_sorted]][1];
                                    int tri_C_idx = faces_list_all[added_tri_list[j_sorted]][2];
                                    Eigen::Vector3d tri_A, tri_B, tri_C;
                                    tri_A << pts_list->points[tri_A_idx].x, pts_list->points[tri_A_idx].y,
                                            pts_list->points[tri_A_idx].z;
                                    tri_B << pts_list->points[tri_B_idx].x, pts_list->points[tri_B_idx].y,
                                            pts_list->points[tri_B_idx].z;
                                    tri_C << pts_list->points[tri_C_idx].x, pts_list->points[tri_C_idx].y,
                                            pts_list->points[tri_C_idx].z;
                                    Eigen::Vector3d tri_ab = tri_B - tri_A;
                                    Eigen::Vector3d tri_ac = tri_C - tri_A;
                                    Eigen::Vector3d tri_n = tri_ab.cross(tri_ac);

                                    Eigen::Vector3d n_nc = cur_v1.cross(cur_v2);

                                    bool if_acw_nc = true;
                                    if (n_nc.dot(tri_n) < 0)
                                    {
                                        if_acw_nc = false;
                                    }

                                    // 判断法向是否符合要求  TODO
                                    if (if_acw_nc)
                                    {
                                        //if (if_reasonable_normal(current_idx, A_idx, B_idx))
                                        if (1)
                                        {
                                            // 添加三角形
                                            std::array<int, 3> temp_add_t;
                                            temp_add_t[0] = current_idx;
                                            temp_add_t[1] = idx_1;
                                            temp_add_t[2] = idx_2;

                                            mtx_faces_list_all.lock();
                                            added_tri_list_nc.push_back(faces_list_all.size());
                                            faces_list_all.push_back(temp_add_t);
                                            mtx_faces_list_all.unlock();
                                        }
                                    }
                                    else
                                    {
                                        //if (if_reasonable_normal(current_idx, B_idx, A_idx))
                                        if (1)
                                        {
                                            // 添加三角形
                                            std::array<int, 3> temp_add_t;
                                            temp_add_t[0] = current_idx;
                                            temp_add_t[1] = idx_2;
                                            temp_add_t[2] = idx_1;

                                            mtx_faces_list_all.lock();
                                            added_tri_list_nc.push_back(faces_list_all.size());
                                            faces_list_all.push_back(temp_add_t);
                                            mtx_faces_list_all.unlock();
                                        }
                                    }
                                }
                            }
                        }

                        // 更新点-三角形索引
                        for (int j = 0; j < added_tri_list_nc.size(); j++)
                        {
                            added_tri_list.push_back(added_tri_list_nc[j]);
                        }

                        mtx_pt2face_idx.lock();
                        std::set<int> vertex_to_update_state;
                        for (int j = 0; j < added_tri_list.size(); j++)
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                int temp_id = faces_list_all[added_tri_list[j]][l];
                                pt2face_idx[temp_id].push_back(added_tri_list[j]);
                                vertex_to_update_state.insert(temp_id);
                            }
                        }
                        for (int j = 0; j < deleted_tri_list.size(); j++)
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                int temp_id = faces_list_all[deleted_tri_list[j]][l];
                                pt2face_idx[temp_id].erase(std::remove(pt2face_idx[temp_id].begin(), pt2face_idx[temp_id].end(), deleted_tri_list[j]), pt2face_idx[temp_id].end());
                            }
                        }
                        mtx_pt2face_idx.unlock();

                        // 更新点状态
                        if (added_tri_list.size() == 0)
                        {
                            mtx_pts_state.lock();
                            pts_state[current_idx] = 0;
                            mtx_pts_state.unlock();
                        }
                        else
                        {
                            if (if_edge_vertex(current_idx))  // TODO 为什么之前没有判断？
                            {
                                mtx_pts_state.lock();
                                pts_state[current_idx] = 2;   // 边缘点
                                mtx_pts_state.unlock();

                                mtx_tri_edge_v_list_all.lock();
                                tri_edge_v_list_all.push_back(current_idx);
                                mtx_tri_edge_v_list_all.unlock();
                            }
                            else
                            {
                                mtx_pts_state.lock();
                                pts_state[current_idx] = 1;   // 边缘点
                                mtx_pts_state.unlock();
                            }
                        }
                    }
                }

            }
        }

        // 第二次遍历，处理边缘三角形顶点，尽可能使其成为完成状态，同时调整其连接的三角形使尽可能规则
        for (int i = 0; i < voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->all_pts_idx.size(); i++)
        {
            int current_idx = voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->all_pts_idx[i];
            // 若不是边缘点则跳过
            if (pts_state[current_idx] != 2)
            {
                continue;
            }
            mtx_pts_process_num.lock();
            pts_process_num[current_idx]++;
            mtx_pts_process_num.unlock();

            std::vector<int> adj_edge_v_list, abnormal_v;
            find_edge_and_abnormal_vertex(current_idx, adj_edge_v_list, abnormal_v);

            if (abnormal_v.size() > 0)
            {
                // 删除该顶点
                mtx_pts_state.lock();
                pts_state[current_idx] = 3;
                mtx_pts_state.unlock();

                // 删除错误连接的三角形
                for (auto it : pt2face_idx[current_idx])
                {
                    mtx_faces_to_delete.lock();
                    faces_to_delete.push_back(it);
                    mtx_faces_to_delete.unlock();

                    for (int l = 0; l < 3; l++)
                    {
                        int temp_id = faces_list_all[it][l];
                        if (temp_id == current_idx)
                        {
                            continue;
                        }
                        mtx_pt2face_idx.lock();
                        pt2face_idx[temp_id].erase(std::remove(pt2face_idx[temp_id].begin(), pt2face_idx[temp_id].end(), it), pt2face_idx[temp_id].end());
                        mtx_pt2face_idx.unlock();

                        if (pt2face_idx[temp_id].size() == 0)
                        {
                            mtx_pts_state.lock();
                            pts_state[temp_id] = 0;
                            mtx_pts_state.unlock();
                        }
                        else
                        {
                            mtx_pts_state.lock();
                            pts_state[temp_id] = 2;
                            mtx_pts_state.unlock();

                            mtx_tri_edge_v_list_all.lock();
                            tri_edge_v_list_all.push_back(temp_id);
                            mtx_tri_edge_v_list_all.unlock();

                        }
                    }
                }
                mtx_pt2face_idx.lock();
                pt2face_idx[current_idx].clear();
                mtx_pt2face_idx.unlock();
                continue;
            }

            // 两个相同的三角形叠加状态
            if (adj_edge_v_list.size() == 0)
            {
                mtx_pts_state.lock();
                pts_state[current_idx] = 1;
                mtx_pts_state.unlock();
                if (pt2face_idx[current_idx].size() == 2)
                {
                    // 删除该顶点
                    mtx_pts_state.lock();
                    pts_state[current_idx] = 3;
                    mtx_pts_state.unlock();
                    // 删除错误连接的三角形
                    for (auto it : pt2face_idx[current_idx])
                    {
                        mtx_faces_to_delete.lock();
                        faces_to_delete.push_back(it);
                        mtx_faces_to_delete.unlock();

                        for (int l = 0; l < 3; l++)
                        {
                            int temp_id = faces_list_all[it][l];
                            if (temp_id == current_idx)
                            {
                                continue;
                            }
                            mtx_pt2face_idx.lock();
                            pt2face_idx[temp_id].erase(std::remove(pt2face_idx[temp_id].begin(), pt2face_idx[temp_id].end(), it), pt2face_idx[temp_id].end());
                            mtx_pt2face_idx.unlock();
                            if (pt2face_idx[temp_id].size() == 0)
                            {
                                mtx_pts_state.lock();
                                pts_state[temp_id] = 0;
                                mtx_pts_state.unlock();
                            }
                            else
                            {
                                mtx_pts_state.lock();
                                pts_state[temp_id] = 2;
                                mtx_pts_state.unlock();
                                mtx_tri_edge_v_list_all.lock();
                                tri_edge_v_list_all.push_back(temp_id);
                                mtx_tri_edge_v_list_all.unlock();
                            }
                        }
                    }
                    mtx_pt2face_idx.lock();
                    pt2face_idx[current_idx].clear();
                    mtx_pt2face_idx.unlock();
                }
            }

            // 连接两条边缘边的边缘顶点，判断该顶点+两条边缘边是否可以组成新的三角形
            if (adj_edge_v_list.size() == 2)
            {
                // 判断当前点与两个边缘点组成的向量之间的夹角
                Eigen::Vector3d v1, v2, c_v1, c_v2, v1_2, cur_pt_eigen;
                cur_pt_eigen << pts_list->points[current_idx].x, pts_list->points[current_idx].y,
                        pts_list->points[current_idx].z;
                v1 << pts_list->points[adj_edge_v_list[0]].x, pts_list->points[adj_edge_v_list[0]].y,
                        pts_list->points[adj_edge_v_list[0]].z;
                v2 << pts_list->points[adj_edge_v_list[1]].x, pts_list->points[adj_edge_v_list[1]].y,
                        pts_list->points[adj_edge_v_list[1]].z;
                c_v1 = v1 - cur_pt_eigen;
                c_v2 = v2 - cur_pt_eigen;
                v1_2 = v2 - v1;

                bool if_over_180 = false;
                if (pt2face_idx[current_idx].size() == 1)
                {
                    if_over_180 = true;
                }
                else
                {
                    // 根据某个自由边与两侧边的叉乘的点积，判断夹角是否大于180度
                    int v3_with_v1_idx = -1;
                    for (auto it : pt2face_idx[current_idx])
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            if (faces_list_all[it][l] == adj_edge_v_list[0])
                            {
                                for (int ll = 0; ll < 3; ll++)
                                {
                                    if (faces_list_all[it][ll] != current_idx && faces_list_all[it][ll] != adj_edge_v_list[0])
                                    {
                                        v3_with_v1_idx = faces_list_all[it][ll];
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                        if (v3_with_v1_idx != -1) { break; }
                    }
                    Eigen::Vector3d v3_with_v1, v1_cross_v3, v1_cross_v2, c_v3_with_v1;
                    v3_with_v1 << pts_list->points[v3_with_v1_idx].x, pts_list->points[v3_with_v1_idx].y,
                            pts_list->points[v3_with_v1_idx].z;
                    c_v3_with_v1 = v3_with_v1 - cur_pt_eigen;
                    v1_cross_v3 = c_v1.cross(c_v3_with_v1);
                    v1_cross_v2 = c_v1.cross(c_v2);
                    if (v1_cross_v2.dot(v1_cross_v3) > 0)
                    {
                        if_over_180 = true;
                    }
                }

                // 夹角大于180度
                if ( if_over_180 )
                {
                    std::vector<int> nearest_vertexs;
                    search_edge_vertex_under_distance_r(current_idx, nearest_vertexs);

                    // 确定距离r内顶点关联的三角形列表
                    std::set<int> voxel_faces;
                    for (auto it_vertex : nearest_vertexs)  // TODO 如何利用上顶点的已完成/未完成状态
                    {
                        for (auto it_face : pt2face_idx[it_vertex])
                        {
                            voxel_faces.insert(it_face);
                        }
                    }

                    // 建立边-三角形索引
                    std::unordered_map<std::string, int> edge_tri_idx;
                    std::vector<std::vector<std::array<int, 3>>> edge_tri_list;
                    std::vector<std::vector<int>> edge_tri_id_list;
                    for (auto it : voxel_faces)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            int a, b, c;
                            if (j ==0 )
                            {
                                a = faces_list_all[it][0];
                                b = faces_list_all[it][1];
                                c = faces_list_all[it][2];
                            }
                            else if (j == 1)
                            {
                                a = faces_list_all[it][1];
                                b = faces_list_all[it][2];
                                c = faces_list_all[it][0];
                            }
                            else
                            {
                                a = faces_list_all[it][2];
                                b = faces_list_all[it][0];
                                c = faces_list_all[it][1];
                            }
                            if (a > b)
                            {
                                int temp_int = a;
                                a = b;
                                b = temp_int;
                            }
                            std::string edge_key = std::to_string(a) + ' ' + std::to_string(b);
                            if (edge_tri_idx.find(edge_key) == edge_tri_idx.end())
                            {
                                edge_tri_idx[edge_key] = edge_tri_list.size();
                                std::vector<std::array<int, 3>> temp_a;
                                std::array<int, 3> temp_b;
                                temp_b[0] = a;
                                temp_b[1] = b;
                                temp_b[2] = c;
                                temp_a.push_back(temp_b);
                                edge_tri_list.push_back(temp_a);

                                std::vector<int> temp_c;
                                temp_c.push_back(it);
                                edge_tri_id_list.push_back(temp_c);

                            }
                            else
                            {
                                int temp_idx = edge_tri_idx[edge_key];
                                std::array<int, 3> temp_b;
                                temp_b[0] = a;
                                temp_b[1] = b;
                                temp_b[2] = c;
                                edge_tri_list[temp_idx].push_back(temp_b);
                            }
                        }
                    }
                    // 确定边缘边，格式 A-B-C
                    std::vector<std::array<int, 3>> edge_edge_inner;
                    std::vector<int> edge_edge_inner_tri_id;
                    for (int j = 0; j < edge_tri_list.size(); j++)
                    {
                        if (edge_tri_list[j].size() == 1)
                        {
                            // 进一步判断是否真正边缘边，由于后面会判断是否异侧，这里不处理也可以！
                            int edge_v1 = edge_tri_list[j][0][0];
                            int edge_v2 = edge_tri_list[j][0][1];
                            // 去除包含当前点的边
                            if (edge_v1 == current_idx || edge_v2 == current_idx || edge_tri_list[j][0][2] == current_idx)
                            {
                                continue;
                            }
                            int common_tri_num = 0;
                            for (auto it : pt2face_idx[edge_v1])
                            {
                                for (int l = 0; l < 3; l++)
                                {
                                    if (faces_list_all[it][l] == edge_v2)
                                    {
                                        common_tri_num++;
                                    }
                                }
                            }
                            if (common_tri_num == 1)
                            {
                                edge_edge_inner.push_back(edge_tri_list[j][0]);
                                edge_edge_inner_tri_id.push_back(edge_tri_id_list[j][0]);
                            }
                        }
                    }

                    // 依次处理每个边缘边
                    std::vector<int> added_tri_list, deleted_tri_list;
                    for (int j = 0; j < edge_edge_inner.size(); j++)
                    {
                        // 判断相邻三角形法向与lidar-point方向  TODO 改为判断两个三角形的法向夹角
                        int tri_A_idx = faces_list_all[edge_edge_inner_tri_id[j]][0];
                        int tri_B_idx = faces_list_all[edge_edge_inner_tri_id[j]][1];
                        int tri_C_idx = faces_list_all[edge_edge_inner_tri_id[j]][2];
                        Eigen::Vector3d tri_A, tri_B, tri_C;
                        tri_A << pts_list->points[tri_A_idx].x, pts_list->points[tri_A_idx].y,
                                pts_list->points[tri_A_idx].z;
                        tri_B << pts_list->points[tri_B_idx].x, pts_list->points[tri_B_idx].y,
                                pts_list->points[tri_B_idx].z;
                        tri_C << pts_list->points[tri_C_idx].x, pts_list->points[tri_C_idx].y,
                                pts_list->points[tri_C_idx].z;
                        Eigen::Vector3d tri_ab = tri_B - tri_A;
                        Eigen::Vector3d tri_ac = tri_C - tri_A;
                        Eigen::Vector3d tri_n = tri_ab.cross(tri_ac);

                        // 新三角形的相关点
                        int A_idx = edge_edge_inner[j][0];
                        int B_idx = edge_edge_inner[j][1];
                        int C_idx = edge_edge_inner[j][2];
                        Eigen::Vector3d A, B, C;
                        A << pts_list->points[A_idx].x, pts_list->points[A_idx].y,
                                pts_list->points[A_idx].z;
                        B << pts_list->points[B_idx].x, pts_list->points[B_idx].y,
                                pts_list->points[B_idx].z;
                        C << pts_list->points[C_idx].x, pts_list->points[C_idx].y,
                                pts_list->points[C_idx].z;
                        Eigen::Vector3d cur_A = A - cur_pt_eigen;
                        Eigen::Vector3d cur_B = B - cur_pt_eigen;
                        Eigen::Vector3d A_B = B - A;
                        Eigen::Vector3d A_cur = cur_pt_eigen - A;
                        Eigen::Vector3d B_cur = cur_pt_eigen - B;
                        Eigen::Vector3d B_A = A - B;

                        // 判断边长是否符合要求
                        if ( cur_A.norm() > r_dis * N_r_dis) { continue; }
                        if ( cur_B.norm() > r_dis * N_r_dis) { continue; }

                        // 判断三角形法向投影是否在边的两侧
                        if (!if_different_side(A, B, C, cur_pt_eigen))
                        {
                            continue;
                        }

                        // 若夹角过大则，删除再添加
                        double cos_cur_AB = cur_A.dot(cur_B) / cur_A.norm() / cur_B.norm();
                        if (cos_cur_AB < -0.5)  // 120-180  -0.766
                        {
                            // 删除三角形
                            mtx_faces_to_delete.lock();
                            faces_to_delete.push_back(edge_edge_inner_tri_id[j]);
                            mtx_faces_to_delete.unlock();
                            deleted_tri_list.push_back(edge_edge_inner_tri_id[j]);
                            // 确定第三个点的索引
                            int v3_idx;
                            if (tri_A_idx != A_idx && tri_A_idx != B_idx)
                            {
                                v3_idx = tri_A_idx;
                            }
                            else if (tri_B_idx != A_idx && tri_B_idx != B_idx)
                            {
                                v3_idx = tri_B_idx;
                            }
                            else
                            {
                                v3_idx = tri_C_idx;
                            }
                            // 添加三角形
                            Eigen::Vector3d v3, cur_v3, n;
                            v3 << pts_list->points[v3_idx].x, pts_list->points[v3_idx].y,
                                    pts_list->points[v3_idx].z;
                            cur_v3 = v3 - cur_pt_eigen;
                            // 第一个三角形
                            n = cur_v3.cross(cur_A);
                            bool if_acw_nc = true;
                            if (n.dot(tri_n) < 0)
                            {
                                if_acw_nc = false;
                            }
                            if (if_acw_nc)
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = current_idx;
                                temp_add_t[1] = v3_idx;
                                temp_add_t[2] = A_idx;

                                mtx_faces_list_all.lock();
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                                mtx_faces_list_all.unlock();
                            }
                            else
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = current_idx;
                                temp_add_t[1] = A_idx;
                                temp_add_t[2] = v3_idx;

                                mtx_faces_list_all.lock();
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                                mtx_faces_list_all.unlock();
                            }
                            // 第二个三角形
                            n = cur_v3.cross(cur_B);
                            if_acw_nc = true;
                            if (n.dot(tri_n) < 0)
                            {
                                if_acw_nc = false;
                            }
                            if (if_acw_nc)
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = current_idx;
                                temp_add_t[1] = v3_idx;
                                temp_add_t[2] = B_idx;

                                mtx_faces_list_all.lock();
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                                mtx_faces_list_all.unlock();
                            }
                            else
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = current_idx;
                                temp_add_t[1] = B_idx;
                                temp_add_t[2] = v3_idx;

                                mtx_faces_list_all.lock();
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                                mtx_faces_list_all.unlock();
                            }
                            continue;
                        }


                        // 判断角度是否符合要求
                        if (cos_cur_AB > ang_thre) { continue; }
                        if (A_B.dot(A_cur) / A_B.norm() / A_cur.norm() > ang_thre) { continue; }
                        if (B_A.dot(B_cur) / B_A.norm() / B_cur.norm() > ang_thre) { continue; }

                        // TODO 判断可见性

                        // 根据lidar射线/相邻三角形法向确定三个点的逆时针顺序
                        Eigen::Vector3d norm = cur_A.cross(cur_B);
                        bool if_acw = true;
                        if (norm.dot(tri_n) < 0)
                        {
                            if_acw = false;
                        }

                        // 判断法向是否符合要求  TODO
                        if (if_acw)
                        {
                            //if (if_reasonable_normal(current_idx, A_idx, B_idx))
                            if (1)
                            {
                                // 添加三角形
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = current_idx;
                                temp_add_t[1] = A_idx;
                                temp_add_t[2] = B_idx;

                                mtx_faces_list_all.lock();
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                                mtx_faces_list_all.unlock();
                            }
                        }
                        else
                        {
                            //if (if_reasonable_normal(current_idx, B_idx, A_idx))
                            if (1)
                            {
                                // 添加三角形
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = current_idx;
                                temp_add_t[1] = B_idx;
                                temp_add_t[2] = A_idx;

                                mtx_faces_list_all.lock();
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                                mtx_faces_list_all.unlock();
                            }
                        }
                    }

                    mtx_pt2face_idx.lock();
                    for (int j = 0; j < added_tri_list.size(); j++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            int temp_id = faces_list_all[added_tri_list[j]][l];
                            pt2face_idx[temp_id].push_back(added_tri_list[j]);
                        }
                    }
                    for (int j = 0; j < deleted_tri_list.size(); j++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            int temp_id = faces_list_all[deleted_tri_list[j]][l];
                            pt2face_idx[temp_id].erase(std::remove(pt2face_idx[temp_id].begin(), pt2face_idx[temp_id].end(), deleted_tri_list[j]), pt2face_idx[temp_id].end());
                        }
                    }
                    mtx_pt2face_idx.unlock();


                }
                // 夹角小于180度
                else
                {

                    std::vector<int> tri_idx_add, tri_idx_delete;
                    double cos_c_v1_v2 = c_v1.dot(c_v2) / c_v1.norm() / c_v2.norm();
                    // 夹角过大 120 - 180
                    if (cos_c_v1_v2 < -0.5)
                    {
                        // 判断是否存在包含另外两个边缘点的三角形
                        int target_tri = -1;
                        for (auto it : pt2face_idx[adj_edge_v_list[0]])
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                if (faces_list_all[it][l] == adj_edge_v_list[1])
                                {
                                    target_tri = it;
                                    break;
                                }
                            }
                            if (target_tri != -1) { break; }
                        }
                        if (target_tri != -1)
                        {
                            mtx_faces_to_delete.lock();
                            faces_to_delete.push_back(target_tri);
                            mtx_faces_to_delete.unlock();
                            tri_idx_delete.push_back(target_tri);
                            // 要删除的三角形的三个点，计算法向
                            int A_idx = faces_list_all[target_tri][0];
                            int B_idx = faces_list_all[target_tri][1];
                            int C_idx = faces_list_all[target_tri][2];
                            Eigen::Vector3d A, B, C;
                            A << pts_list->points[A_idx].x, pts_list->points[A_idx].y,
                                    pts_list->points[A_idx].z;
                            B << pts_list->points[B_idx].x, pts_list->points[B_idx].y,
                                    pts_list->points[B_idx].z;
                            C << pts_list->points[C_idx].x, pts_list->points[C_idx].y,
                                    pts_list->points[C_idx].z;
                            Eigen::Vector3d tri_ab = B - A;
                            Eigen::Vector3d tri_ac = C - A;
                            Eigen::Vector3d tri_n = tri_ab.cross(tri_ac);
                            // 确定第四个点
                            int v3_idx;
                            if (A_idx != adj_edge_v_list[0] && A_idx != adj_edge_v_list[1])
                            {
                                v3_idx = A_idx;
                            }
                            else if (B_idx != adj_edge_v_list[0] && B_idx != adj_edge_v_list[1])
                            {
                                v3_idx = B_idx;
                            }
                            else
                            {
                                v3_idx = C_idx;
                            }
                            // 添加三角形
                            Eigen::Vector3d v3, c_v3, n;
                            v3 << pts_list->points[v3_idx].x, pts_list->points[v3_idx].y,
                                    pts_list->points[v3_idx].z;
                            c_v3 = v3 - cur_pt_eigen;
                            // 第一个三角形
                            n = c_v1.cross(c_v3);
                            bool if_acw_nc = true;
                            if (n.dot(tri_n) < 0)
                            {
                                if_acw_nc = false;
                            }
                            if (if_acw_nc)
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = current_idx;
                                temp_add_t[1] = adj_edge_v_list[0];
                                temp_add_t[2] = v3_idx;

                                mtx_faces_list_all.lock();
                                tri_idx_add.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                                mtx_faces_list_all.unlock();
                            }
                            else
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = current_idx;
                                temp_add_t[1] = v3_idx;
                                temp_add_t[2] = adj_edge_v_list[0];

                                mtx_faces_list_all.lock();
                                tri_idx_add.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                                mtx_faces_list_all.unlock();
                            }
                            // 第二个三角形
                            n = c_v2.cross(c_v3);
                            if_acw_nc = true;
                            if (n.dot(tri_n) < 0)
                            {
                                if_acw_nc = false;
                            }
                            if (if_acw_nc)
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = current_idx;
                                temp_add_t[1] = adj_edge_v_list[1];
                                temp_add_t[2] = v3_idx;

                                mtx_faces_list_all.lock();
                                tri_idx_add.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                                mtx_faces_list_all.unlock();
                            }
                            else
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = current_idx;
                                temp_add_t[1] = v3_idx;
                                temp_add_t[2] = adj_edge_v_list[1];

                                mtx_faces_list_all.lock();
                                tri_idx_add.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                                mtx_faces_list_all.unlock();
                            }
                        }
                    }
                        // 夹角过小 0 - 30
                    else if (cos_c_v1_v2 > 0.866)
                    {
                        // 找v1、v2中的较小角、较大角
                        int min_ang_pt_idx, max_ang_pt_idx;
                        if (c_v2.dot(v1_2) < 0)
                        {
                            min_ang_pt_idx = adj_edge_v_list[0];
                            max_ang_pt_idx = adj_edge_v_list[1];
                        }
                        else
                        {
                            min_ang_pt_idx = adj_edge_v_list[1];
                            max_ang_pt_idx = adj_edge_v_list[0];
                        }
                        // 找要删除的三角形索引
                        int target_tri = -1;
                        for (auto it : pt2face_idx[current_idx])
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                if (faces_list_all[it][l] == min_ang_pt_idx)
                                {
                                    target_tri = it;
                                    break;
                                }
                            }
                            if (target_tri != -1) { break; }
                        }
                        assert(target_tri != -1);
                        // 添加删除三角形
                        mtx_faces_to_delete.lock();
                        faces_to_delete.push_back(target_tri);
                        mtx_faces_to_delete.unlock();
                        tri_idx_delete.push_back(target_tri);
                        // 要删除的三角形的三个点，计算法向
                        int A_idx = faces_list_all[target_tri][0];
                        int B_idx = faces_list_all[target_tri][1];
                        int C_idx = faces_list_all[target_tri][2];
                        Eigen::Vector3d A, B, C;
                        A << pts_list->points[A_idx].x, pts_list->points[A_idx].y,
                                pts_list->points[A_idx].z;
                        B << pts_list->points[B_idx].x, pts_list->points[B_idx].y,
                                pts_list->points[B_idx].z;
                        C << pts_list->points[C_idx].x, pts_list->points[C_idx].y,
                                pts_list->points[C_idx].z;
                        Eigen::Vector3d tri_ab = B - A;
                        Eigen::Vector3d tri_ac = C - A;
                        Eigen::Vector3d tri_n = tri_ab.cross(tri_ac);
                        // 确定第四个点
                        int v3_idx;
                        if (A_idx != current_idx && A_idx != min_ang_pt_idx)
                        {
                            v3_idx = A_idx;
                        }
                        else if (B_idx != current_idx && B_idx != min_ang_pt_idx)
                        {
                            v3_idx = B_idx;
                        }
                        else
                        {
                            v3_idx = C_idx;
                        }
                        // 添加三角形
                        Eigen::Vector3d v3, min_pt, max_pt, v3_c, v3_min, v3_max, n;
                        v3 << pts_list->points[v3_idx].x, pts_list->points[v3_idx].y,
                                pts_list->points[v3_idx].z;
                        min_pt << pts_list->points[min_ang_pt_idx].x, pts_list->points[min_ang_pt_idx].y,
                                pts_list->points[min_ang_pt_idx].z;
                        max_pt << pts_list->points[max_ang_pt_idx].x, pts_list->points[max_ang_pt_idx].y,
                                pts_list->points[max_ang_pt_idx].z;
                        v3_c = cur_pt_eigen - v3;
                        v3_min = min_pt - v3;
                        v3_max = max_pt - v3;
                        // 第一个三角形
                        n = v3_c.cross(v3_max);
                        bool if_acw_nc = true;
                        if (n.dot(tri_n) < 0)
                        {
                            if_acw_nc = false;
                        }
                        if (if_acw_nc)
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = v3_idx;
                            temp_add_t[1] = current_idx;
                            temp_add_t[2] = max_ang_pt_idx;

                            mtx_faces_list_all.lock();
                            tri_idx_add.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                            mtx_faces_list_all.unlock();
                        }
                        else
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = v3_idx;
                            temp_add_t[1] = max_ang_pt_idx;
                            temp_add_t[2] = current_idx;

                            mtx_faces_list_all.lock();
                            tri_idx_add.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                            mtx_faces_list_all.unlock();
                        }
                        // 第二个三角形
                        n = v3_min.cross(v3_max);
                        if_acw_nc = true;
                        if (n.dot(tri_n) < 0)
                        {
                            if_acw_nc = false;
                        }
                        if (if_acw_nc)
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = v3_idx;
                            temp_add_t[1] = min_ang_pt_idx;
                            temp_add_t[2] = max_ang_pt_idx;

                            mtx_faces_list_all.lock();
                            tri_idx_add.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                            mtx_faces_list_all.unlock();
                        }
                        else
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = v3_idx;
                            temp_add_t[1] = max_ang_pt_idx;
                            temp_add_t[2] = min_ang_pt_idx;

                            mtx_faces_list_all.lock();
                            tri_idx_add.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                            mtx_faces_list_all.unlock();
                        }

                    }
                        // 30 - 120
                    else
                    {
                        // 判断三角形边长
                        if (v1_2.norm() > r_dis * N_r_dis) { continue; }     // TODO 调节
//                        if (v1_2.norm() > r_dis) { continue; }
//                    if (c_v1.norm() > r_dis) { continue; }
//                    if (c_v2.norm() > r_dis) { continue; }

                        // 判断三角形角度
                        Eigen::Vector3d v1_c, v2_c, v2_1;
                        v1_c = cur_pt_eigen - v1;
                        v2_c = cur_pt_eigen - v2;
                        v2_1 = v1 - v2;
//                        double ang_thre = 0.866; // 0.866 0.93969;
//                        if (v1_c.dot(v1_2) / v1_c.norm() / v1_2.norm() > ang_thre) { continue; }  // TODO 调节
//                        if (v2_c.dot(v2_1) / v2_c.norm() / v2_1.norm() > ang_thre) { continue; }

                        // 判断有无三角形顶点位于三角形内部
                        std::vector<int> nearest_vertexs;
                        std::vector<double> nearest_vertexs_dis;
                        search_nearest_vertex_under_distance_r(current_idx, nearest_vertexs, nearest_vertexs_dis);
                        bool if_any_v_inner = false;
                        for (int nvi = 0; nvi < nearest_vertexs.size(); nvi++)
                        {
                            if (nearest_vertexs[nvi] == adj_edge_v_list[0] || nearest_vertexs[nvi] == adj_edge_v_list[1])
                            {
                                continue;
                            }
                            Eigen::Vector3d temp_v;
                            temp_v << pts_list->points[nearest_vertexs[nvi]].x, pts_list->points[nearest_vertexs[nvi]].y,
                                    pts_list->points[nearest_vertexs[nvi]].z;
                            if (if_p_in_tri(temp_v, cur_pt_eigen, v1, v2))
                            {
                                if_any_v_inner = true;
                                break;
                            }
                        }
                        if (if_any_v_inner)
                        {
                            continue;
                        }

                        // 确定相邻三角形id
                        int target_tri = -1;
                        for (auto it : pt2face_idx[current_idx])
                        {
                            for (int l = 0; l < 3; l++)
                            {
                                if (faces_list_all[it][l] == adj_edge_v_list[0])
                                {
                                    target_tri = it;
                                    break;
                                }
                            }
                            if (target_tri != -1) { break; }
                        }
                        assert(target_tri != -1);
                        // 计算法向
                        int A_idx = faces_list_all[target_tri][0];
                        int B_idx = faces_list_all[target_tri][1];
                        int C_idx = faces_list_all[target_tri][2];
                        Eigen::Vector3d A, B, C;
                        A << pts_list->points[A_idx].x, pts_list->points[A_idx].y,
                                pts_list->points[A_idx].z;
                        B << pts_list->points[B_idx].x, pts_list->points[B_idx].y,
                                pts_list->points[B_idx].z;
                        C << pts_list->points[C_idx].x, pts_list->points[C_idx].y,
                                pts_list->points[C_idx].z;
                        Eigen::Vector3d tri_ab = B - A;
                        Eigen::Vector3d tri_ac = C - A;
                        Eigen::Vector3d tri_n = tri_ab.cross(tri_ac);
                        // 添加三角形
                        Eigen::Vector3d n;
                        n = c_v1.cross(c_v2);
                        bool if_acw_nc = true;
                        if (n.dot(tri_n) < 0)
                        {
                            if_acw_nc = false;
                        }
                        if (if_acw_nc)
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = current_idx;
                            temp_add_t[1] = adj_edge_v_list[0];
                            temp_add_t[2] = adj_edge_v_list[1];

                            mtx_faces_list_all.lock();
                            tri_idx_add.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                            mtx_faces_list_all.unlock();
                        }
                        else
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = current_idx;
                            temp_add_t[1] = adj_edge_v_list[1];
                            temp_add_t[2] = adj_edge_v_list[0];

                            mtx_faces_list_all.lock();
                            tri_idx_add.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                            mtx_faces_list_all.unlock();
                        }
                    }

                    // 更新点-三角形索引列表
                    mtx_pt2face_idx.lock();
                    for (int j = 0; j < tri_idx_delete.size(); j++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            int temp_id = faces_list_all[tri_idx_delete[j]][l];
                            pt2face_idx[temp_id].erase(std::remove(pt2face_idx[temp_id].begin(), pt2face_idx[temp_id].end(), tri_idx_delete[j]), pt2face_idx[temp_id].end());
                        }
                    }
                    for (int j = 0; j < tri_idx_add.size(); j++)
                    {
                        for (int l = 0; l < 3; l++)
                        {
                            int temp_id = faces_list_all[tri_idx_add[j]][l];
                            pt2face_idx[temp_id].push_back(tri_idx_add[j]);
                        }
                    }
                    mtx_pt2face_idx.unlock();
                }
            }
        }

        voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->if_need_update = false;
//        voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->new_pt_num = 0;
    }
}

// 判断某个世界坐标系下的点是否投影在某个平面的内部
bool NonPlaneMesh::if_project_inner_plane(const Eigen::Vector3d& pt, const std::shared_ptr<MeshFragment>& plane)
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
        double y_start = plane->y0 - plane->dis_resolution * plane->leftW;
        double z_start = plane->z0 + plane->dis_resolution * plane->topH;
        int u = int((z_start - pt_pro_plane(2)) / plane->dis_resolution);
        int v = int((pt_pro_plane(1) - y_start) / plane->dis_resolution);
        int H = plane->bottomH + plane->topH;
        int W = plane->leftW + plane->rightW;
        if (u < 0 || u > H - 1 || v < 0 || v > W - 1)
        {
            // 投影在平面外
            return false;
        }
        if (plane->uniform_grid[u][v].node_type == 0)
        {
            // 投影网格还不是内部点
            return false;
        }
        if (std::abs(pt_pro_plane(0) - plane->uniform_grid[u][v].x_avg) > 0.2)  // TODO 阈值
        {
            // 虽投影在平面内，但与平面法向距离太大
            return false;
        }
        if (plane->uniform_grid[u][v].node_type == 2)
        {
            for (auto it : plane->uniform_grid[u][v].tri_idx_list)
            {
                int A_idx = plane->faces_list[it][0];
                int B_idx = plane->faces_list[it][1];
                int C_idx = plane->faces_list[it][2];
                Eigen::Vector3d A, B, C;
                A << plane->ptcl_grid->points[A_idx].x, plane->ptcl_grid->points[A_idx].y,
                        plane->ptcl_grid->points[A_idx].z;
                B << plane->ptcl_grid->points[B_idx].x, plane->ptcl_grid->points[B_idx].y,
                        plane->ptcl_grid->points[B_idx].z;
                C << plane->ptcl_grid->points[C_idx].x, plane->ptcl_grid->points[C_idx].y,
                        plane->ptcl_grid->points[C_idx].z;

                if ( if_p_in_tri(pt, A, B, C) )
                {
                    return true;
                }
            }
            return false;
        }
    }
    return true;

}

void NonPlaneMesh::mesh_update_from_plane_map(const std::shared_ptr<MeshFragment>& plane)
{
    // 删除边缘点相关的三角形，点-三角形索引，修改边缘点状态
    std::vector<int> relevant_v_idx;
    for (auto pt_idx : plane->pts_edge_delete)
    {
        if (plane->plane_vidx_to_noplnae_vidx.find(pt_idx) == plane->plane_vidx_to_noplnae_vidx.end())
        {
            std::cout << "\033[31m Warning: No edge vertex index!  \033[0m" << std::endl;
            continue;
        }
        int pt_idx_np = plane->plane_vidx_to_noplnae_vidx[pt_idx];  // np中点的索引

        for (auto tri_idx : pt2face_idx[pt_idx_np])
        {
            // 添加到删除三角形列表
            faces_to_delete.push_back(tri_idx);

            for (int l = 0; l < 3; l++)
            {
                // 当前点、平面内部点不处理
                int temp_p_id = faces_list_all[tri_idx][l];
                if (temp_p_id == pt_idx_np || pts_state[temp_p_id] == 4)
                {
                    continue;
                }
                // 删除点-三角形索引，可能是边缘点、非平面点 TODO 若是边缘点有没有影响？？？
                pt2face_idx[temp_p_id].erase(std::remove(pt2face_idx[temp_p_id].begin(), pt2face_idx[temp_p_id].end(),
                                                         tri_idx), pt2face_idx[temp_p_id].end());

                // 添加相连接的非平面点
                if (pts_state[temp_p_id] == 5 || pts_state[temp_p_id] == 6)
                {
                    pts_state[temp_p_id] = 5;
                    continue;
                }
                auto temp_v_idx = std::find(relevant_v_idx.begin(), relevant_v_idx.end(), temp_p_id);
                if (temp_v_idx != relevant_v_idx.end())
                {
                    continue; // 已经添加过了
                }
                relevant_v_idx.push_back(temp_p_id);
            }
        }
        // 删除边缘点-三角形索引
        pt2face_idx[pt_idx_np].clear();
        // 修改边缘点的状态为 3-删除
        pts_state[pt_idx_np] = 3;
    }

    // 处理边缘点连接的非平面点、与之相连的三角形
    for (auto v_idx : relevant_v_idx)
    {
        // 判断是否投影在平面内部
        Eigen::Vector3d temp_pt;
        temp_pt << pts_list->points[v_idx].x, pts_list->points[v_idx].y, pts_list->points[v_idx].z;
        // 若投影不在平面内部，则只修改点状态
        if ( !if_project_inner_plane(temp_pt, plane) )
        {
            if (pt2face_idx[v_idx].size() > 0)
            {
                pts_state[v_idx] = 2;
                tri_edge_v_list_all.push_back(v_idx);
            }
            else
            {
                pts_state[v_idx] = 0;  // TODO 0 or 3?
            }
        }
        // 若在平面内部，则删除点、三角形、点-三角形索引 TODO 是否要广度优先扩散去删除？
        else
        {
            for (auto tri_idx : pt2face_idx[v_idx])
            {
                // 添加到删除三角形列表
                faces_to_delete.push_back(tri_idx);

                // 删除其他点-三角形索引
                for (int l = 0; l < 3; l++)
                {
                    // 当前点不处理
                    int temp_p_id = faces_list_all[tri_idx][l];
                    if (temp_p_id == v_idx)
                    {
                        continue;
                    }
                    // 删除点-三角形索引，可能是边缘点、非平面点
                    pt2face_idx[temp_p_id].erase(std::remove(pt2face_idx[temp_p_id].begin(), pt2face_idx[temp_p_id].end(),
                                                             tri_idx), pt2face_idx[temp_p_id].end());
                    if (pts_state[temp_p_id] == 6)
                    {
                        pts_state[temp_p_id] = 5;
                    }
                    else if (pts_state[temp_p_id] == 1 || pts_state[temp_p_id] == 2)
                    {
                        if (pt2face_idx[temp_p_id].size() > 0)
                        {
                            pts_state[temp_p_id] = 2;
                            tri_edge_v_list_all.push_back(temp_p_id);
                        }
                        else
                        {
                            pts_state[temp_p_id] = 0;
                        }
                    }
                    else
                    {
                        // 其余不处理
                    }
                }
            }
            // 删除邻近非平面点-三角形索引
            pt2face_idx[v_idx].clear();
            // 修改邻近非平面点的状态为 3-删除
            pts_state[v_idx] = 3;
        }
    }

    // 根据新加入的平面内部点，删除投影在平面内部的点和三角形
    for (auto pt_idx : plane->pts_inner_add)
    {
        // 平面内部点在世界坐标系下的坐标
        Eigen::Vector3d temp_pt;
        temp_pt << plane->ptcl_grid->points[pt_idx].x,
                plane->ptcl_grid->points[pt_idx].y,
                plane->ptcl_grid->points[pt_idx].z;

        // 搜索平面内部点r距离内的三角形顶点，r=plane->dis_resolution  TODO 阈值
        std::vector<int> near_v_list;
        search_nearest_valid_ptv(temp_pt, plane->dis_resolution, near_v_list);

        // 依次处理每个最近邻三角形顶点
        for (auto v_idx : near_v_list)
        {
            // 判断三角形顶点是否投影在平面内部
            Eigen::Vector3d temp_v;
            temp_v << pts_list->points[v_idx].x,
                    pts_list->points[v_idx].y,
                    pts_list->points[v_idx].z;
            // 若投影在平面内部
            if ( if_project_inner_plane(temp_v, plane) )
            {
                for (auto tri_idx : pt2face_idx[v_idx])
                {
                    // 添加到删除三角形列表
                    faces_to_delete.push_back(tri_idx);

                    // 删除其他点-三角形索引
                    for (int l = 0; l < 3; l++)
                    {
                        // 当前点、平面内部点不处理
                        int temp_p_id = faces_list_all[tri_idx][l];
                        if (temp_p_id == v_idx)
                        {
                            continue;
                        }
                        // 删除点-三角形索引，可能是边缘点、非平面点
                        pt2face_idx[temp_p_id].erase(std::remove(pt2face_idx[temp_p_id].begin(), pt2face_idx[temp_p_id].end(),
                                                                 tri_idx), pt2face_idx[temp_p_id].end());
                        if (pts_state[temp_p_id] == 6)
                        {
                            pts_state[temp_p_id] = 5;
                        }
                        else if (pts_state[temp_p_id] == 1 || pts_state[temp_p_id] == 2)
                        {
                            if (pt2face_idx[temp_p_id].size() > 0)
                            {
                                pts_state[temp_p_id] = 2;
                                tri_edge_v_list_all.push_back(temp_p_id);
                            }
                            else
                            {
                                pts_state[temp_p_id] = 0;
                            }
                        }
                        else
                        {
                            // 其余不处理
                        }
                    }
                }
                // 删除边缘点-三角形索引
                pt2face_idx[v_idx].clear();
                // 修改边缘点的状态为 3-删除
                pts_state[v_idx] = 3;
            }
        }

    }

    // 更新 需更新边缘点的位置
    for (auto pt_idx : plane->pts_edge_update)
    {
        if (plane->plane_vidx_to_noplnae_vidx.find(pt_idx) == plane->plane_vidx_to_noplnae_vidx.end())
        {
            std::cout << "\033[31m Warning: No edge vertex index!  \033[0m" << std::endl;
            continue;
        }
        int pt_idx_np = plane->plane_vidx_to_noplnae_vidx[pt_idx];  // np中点的索引

        pts_list->points[pt_idx_np].x = plane->ptcl_grid->points[pt_idx].x;
        pts_list->points[pt_idx_np].y = plane->ptcl_grid->points[pt_idx].y;
        pts_list->points[pt_idx_np].z = plane->ptcl_grid->points[pt_idx].z;

        // 判断更新点对应的三角形是否需要删除
        std::vector<int> temp_delete_tri;
        for (auto tri_idx : pt2face_idx[pt_idx_np])
        {
            std::vector<int> temp_tri_v_idx;
            temp_tri_v_idx.push_back(faces_list_all[tri_idx][0]);
            temp_tri_v_idx.push_back(faces_list_all[tri_idx][1]);
            temp_tri_v_idx.push_back(faces_list_all[tri_idx][2]);
            std::sort(temp_tri_v_idx.begin(), temp_tri_v_idx.end());
            std::string tri_key = std::to_string(temp_tri_v_idx[0]) + ' '
                                  + std::to_string(temp_tri_v_idx[1]) + ' '
                                  + std::to_string(temp_tri_v_idx[2]);
            if (plane->tri_delete_string_idx.find(tri_key) != plane->tri_delete_string_idx.end())
            {
                temp_delete_tri.push_back(tri_idx);
            }
        }
        // 删除顶点对该三角形的索引
        for (auto tri_idx : temp_delete_tri)
        {
            for (int l = 0; l < 3; l++)
            {
                // 平面内部点不处理
                int temp_p_id = faces_list_all[tri_idx][l];
                if (pts_state[temp_p_id] == 4) {
                    continue;
                }
                // 删除点-三角形索引，可能是边缘点、非平面点 TODO 若是边缘点有没有影响？？？
                pt2face_idx[temp_p_id].erase(std::remove(pt2face_idx[temp_p_id].begin(), pt2face_idx[temp_p_id].end(),
                                                         tri_idx), pt2face_idx[temp_p_id].end());
                if (pts_state[temp_p_id] == 6)
                {
                    pts_state[temp_p_id] = 5;
                }
                else if (pts_state[temp_p_id] == 1 || pts_state[temp_p_id] == 2)
                {
                    if (pt2face_idx[temp_p_id].size() > 0)
                    {
                        pts_state[temp_p_id] = 2;
                        tri_edge_v_list_all.push_back(temp_p_id);
                    }
                    else
                    {
                        pts_state[temp_p_id] = 0;
                    }
                }
                else
                {
                    // 其余不处理
                }
            }
        }
    }

    // 添加边缘三角形相关点，并设置状态，平面内部点设置一个特殊的状态
    for (auto it : plane->pts_edge_add)
    {
        if (plane->plane_vidx_to_noplnae_vidx.find(it) == plane->plane_vidx_to_noplnae_vidx.end())
        {
            pcl::PointXYZINormal temp_pt;
            temp_pt.x = plane->ptcl_grid->points[it].x;
            temp_pt.y = plane->ptcl_grid->points[it].y;
            temp_pt.z = plane->ptcl_grid->points[it].z;

            pts_list->push_back(temp_pt);
            if (plane->pts_state[it] == 1)
            {
                // 平面内部点
                pts_state.push_back(4);
            }
            if (plane->pts_state[it] == 2)
            {
                pts_state.push_back(5);   // TODO 边缘2？ 平面边缘点5？
            }
            pts_process_num.push_back(0);
            hole_process_num.push_back(0);
            hole_v_process.push_back(0);
            pt2face_idx.resize(pts_list->size());

            // 平面点索引 -》 非平面点索引
            int curpt_idx = pts_list->size() - 1;
            plane->plane_vidx_to_noplnae_vidx[it] = curpt_idx;

            // 记录边缘点的索引
            if (plane->pts_state[it] == 2)
            {
                plane_edge_v_list_all.push_back(curpt_idx);
            }

            // 边缘点添加体素索引
            if (plane->pts_state[it] == 2)
            {
                long voxel_x = std::round(temp_pt.x / voxel_resolution);
                long voxel_y = std::round(temp_pt.y / voxel_resolution);
                long voxel_z = std::round(temp_pt.z / voxel_resolution);
                if (voxel_pts.if_exist(voxel_x, voxel_y, voxel_z))
                {
                    voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->all_pts_idx.push_back(curpt_idx);
                }
                else
                {
                    std::shared_ptr<m_voxel> temp_voxel = std::make_shared<m_voxel>();
                    voxel_pts.insert(voxel_x, voxel_y, voxel_z, temp_voxel);
                    voxel_pts.m_map_3d_hash_map[voxel_x][voxel_y][voxel_z]->all_pts_idx.push_back(curpt_idx);
                }
            }
        }
        else
        {
            int pt_idx_np = plane->plane_vidx_to_noplnae_vidx[it];  // np中点的索引
            pts_list->points[pt_idx_np].x = plane->ptcl_grid->points[it].x;
            pts_list->points[pt_idx_np].y = plane->ptcl_grid->points[it].y;
            pts_list->points[pt_idx_np].z = plane->ptcl_grid->points[it].z;
            if (plane->pts_state[it] == 1)
            {
                // 平面内部点
                pts_state[pt_idx_np] = 4;
            }
            if (plane->pts_state[it] == 2)
            {
                pts_state[pt_idx_np] = 5;   // TODO 边缘2？ 平面边缘点5？
            }
        }
    }

    // 添加边缘三角形，添加点-三角形索引
    std::set<int> state_update_v_list;
    for (auto tri_idx : plane->faces_with_2_edge_vertex)
    {
        std::array<int, 3> temp_tri_idx;
        temp_tri_idx[0] = plane->plane_vidx_to_noplnae_vidx[plane->faces_list[tri_idx][0]];
        temp_tri_idx[1] = plane->plane_vidx_to_noplnae_vidx[plane->faces_list[tri_idx][1]];
        temp_tri_idx[2] = plane->plane_vidx_to_noplnae_vidx[plane->faces_list[tri_idx][2]];

        faces_list_all.push_back(temp_tri_idx);
        // 最终不会保存平面三角形，避免重复
        faces_to_delete.push_back(faces_list_all.size() - 1);

        for (int l = 0; l < 3; l++)
        {
            // 只有平面边缘点才添加三角形索引 TODO 对非平面点处理有没有影响？
            if (plane->pts_state[plane->faces_list[tri_idx][l]] == 2)
            {
                pt2face_idx[temp_tri_idx[l]].push_back(faces_list_all.size() - 1);
                state_update_v_list.insert(temp_tri_idx[l]);
            }
        }
    }
    // 更新新添加三角形相关顶点的状态
    for (auto it : state_update_v_list)
    {
        if (pts_state[it] != 5)
        {
            continue;
        }
        if (!if_edge_vertex(it))
        {
            pts_state[it] = 6;
        }
    }
}

void NonPlaneMesh::mesh_delete_plane_map(const std::shared_ptr<MeshFragment>& plane)
{
    // 删除边缘点相关的三角形，点-三角形索引，修改边缘点状态
    std::vector<int> relevant_v_idx;
    for (auto it : plane->plane_vidx_to_noplnae_vidx)
    {
        int pt_idx_np = it.second;  // np中点的索引

        for (auto tri_idx : pt2face_idx[pt_idx_np])
        {
            // 添加到删除三角形列表
            faces_to_delete.push_back(tri_idx);

            for (int l = 0; l < 3; l++)
            {
                // 当前点、平面内部点不处理
                int temp_p_id = faces_list_all[tri_idx][l];
                if (temp_p_id == pt_idx_np || pts_state[temp_p_id] == 4)
                {
                    continue;
                }
                // 删除点-三角形索引，可能是边缘点、非平面点 TODO 若是边缘点有没有影响？？？
                pt2face_idx[temp_p_id].erase(std::remove(pt2face_idx[temp_p_id].begin(), pt2face_idx[temp_p_id].end(),
                                                         tri_idx), pt2face_idx[temp_p_id].end());

                // 添加相连接的非平面点
                if (pts_state[temp_p_id] == 5 || pts_state[temp_p_id] == 6)
                {
                    pts_state[temp_p_id] = 5;
                    continue;
                }
                auto temp_v_idx = std::find(relevant_v_idx.begin(), relevant_v_idx.end(), temp_p_id);
                if (temp_v_idx != relevant_v_idx.end())
                {
                    continue; // 已经添加过了
                }
                relevant_v_idx.push_back(temp_p_id);
            }
        }
        // 删除边缘点-三角形索引
        pt2face_idx[pt_idx_np].clear();
        // 修改边缘点的状态为 3-删除
        pts_state[pt_idx_np] = 3;
    }

    // 处理边缘点连接的非平面点、与之相连的三角形
    for (auto v_idx : relevant_v_idx)
    {
        // 判断是否投影在平面内部
        Eigen::Vector3d temp_pt;
        temp_pt << pts_list->points[v_idx].x, pts_list->points[v_idx].y, pts_list->points[v_idx].z;
        // 若投影不在平面内部，则只修改点状态
        if ( !if_project_inner_plane(temp_pt, plane) )
        {
            if (pt2face_idx[v_idx].size() > 0)
            {
                pts_state[v_idx] = 2;
                tri_edge_v_list_all.push_back(v_idx);
            }
            else
            {
                pts_state[v_idx] = 0;  // TODO 0 or 3?
            }
        }
            // 若在平面内部，则删除点、三角形、点-三角形索引 TODO 是否要广度优先扩散去删除？
        else
        {
            for (auto tri_idx : pt2face_idx[v_idx])
            {
                // 添加到删除三角形列表
                faces_to_delete.push_back(tri_idx);

                // 删除其他点-三角形索引
                for (int l = 0; l < 3; l++)
                {
                    // 当前点不处理
                    int temp_p_id = faces_list_all[tri_idx][l];
                    if (temp_p_id == v_idx)
                    {
                        continue;
                    }
                    // 删除点-三角形索引，可能是边缘点、非平面点
                    pt2face_idx[temp_p_id].erase(std::remove(pt2face_idx[temp_p_id].begin(), pt2face_idx[temp_p_id].end(),
                                                             tri_idx), pt2face_idx[temp_p_id].end());
                    if (pts_state[temp_p_id] == 6)
                    {
                        pts_state[temp_p_id] = 5;
                    }
                    else if (pts_state[temp_p_id] == 1 || pts_state[temp_p_id] == 2)
                    {
                        if (pt2face_idx[temp_p_id].size() > 0)
                        {
                            pts_state[temp_p_id] = 2;
                            tri_edge_v_list_all.push_back(temp_p_id);
                        }
                        else
                        {
                            pts_state[temp_p_id] = 0;
                        }
                    }
                    else
                    {
                        // 其余不处理
                    }
                }
            }
            // 删除邻近非平面点-三角形索引
            pt2face_idx[v_idx].clear();
            // 修改邻近非平面点的状态为 3-删除
            pts_state[v_idx] = 3;
        }
    }
}

void NonPlaneMesh::plane_and_noplane_connect()
{
    double r_dis_for_edge_v = r_dis * N_r_dis_for_edge_v;
    for (int i = 0; i < plane_edge_v_list_all.size(); i++)
    {
        int cur_edge_v_idx = plane_edge_v_list_all[i];
        // 跳过非边缘点 5-新加入/未连接 3-删除 1-已完成  TODO 1 or 6
        if (pts_state[cur_edge_v_idx] != 5)
        {
            continue;
        }

        // 判断两个邻接边缘点是否都是平面边缘点
        bool if_two_side = false;
        bool if_passable = false;

        // 履带式连接 1-1，哪一边点多，改成 1-n，一般是非平面一侧点多
        // 前一个点连多的点 or 后一个点连多的点 or 一边一半
        while (1)
        {
            // 确定平面边缘点的两个邻接平面边缘点
            std::vector<int> adj_edge_v_list;
            find_edge_vertex(cur_edge_v_idx, adj_edge_v_list);
            if (adj_edge_v_list.size() != 2)
            {
                if (adj_edge_v_list.size() == 0)
                {
                    pts_state[cur_edge_v_idx] = 6;
                }
                break;
            }

            // 判断是否需要继续延展
            int if_need_extend = true;
            // 第一次扩展，判断扩展类型，1-1、1-n、n-1
            bool if_p1_npn = true;

            std::vector<int> plane_edge_v_list, noplane_edge_v_list;
            Eigen::Vector3d cur_eigen_pt;
            cur_eigen_pt << pts_list->points[cur_edge_v_idx].x,
                    pts_list->points[cur_edge_v_idx].y,
                    pts_list->points[cur_edge_v_idx].z;
            // 若两侧都是平面边缘点，则随机选择一侧延展处理  TODO 两侧都处理
            if (pts_state[adj_edge_v_list[0]] == 5 && pts_state[adj_edge_v_list[1]] == 5)
            {
                if_passable = true;
                // 确定一个平面边缘点，r距离内存在一个非平面边缘点
                std::vector<int> nearest_v_list;
                search_edge_vertex_under_distance_kr(cur_edge_v_idx, r_dis_for_edge_v, nearest_v_list);
                if (nearest_v_list.size() == 0)
                {
                    break;
                }
                // 找距离当前平面边缘点最近的点
                Eigen::Vector3d temp_v;
                int nearest_edge_v_idx = nearest_v_list[0];
                temp_v << pts_list->points[nearest_edge_v_idx].x,
                        pts_list->points[nearest_edge_v_idx].y,
                        pts_list->points[nearest_edge_v_idx].z;
                double min_dis = (cur_eigen_pt - temp_v).norm();
                for (int j = 1; j < nearest_v_list.size(); j++)
                {
                    int temp_idx = nearest_v_list[j];
                    temp_v << pts_list->points[temp_idx].x,
                            pts_list->points[temp_idx].y,
                            pts_list->points[temp_idx].z;
                    double temp_dis = (cur_eigen_pt - temp_v).norm();
                    if (temp_dis < min_dis)
                    {
                        min_dis = temp_dis;
                        nearest_edge_v_idx = temp_idx;
                    }
                }
                // 查找距离最近的非平面点的两个邻接非平面点
                std::vector<int> adj_edge_v_list_1;
                find_edge_vertex(nearest_edge_v_idx, adj_edge_v_list_1);
                if (adj_edge_v_list_1.size() != 2)
                {
                    break;
                }
                // 确定沿哪个方向延展
                Eigen::Vector3d p_v1, np_v0, np_v1, np_v2, cur_pv1, npv0_npv1, npv0_npv2;
                p_v1 << pts_list->points[adj_edge_v_list[0]].x,
                        pts_list->points[adj_edge_v_list[0]].y,
                        pts_list->points[adj_edge_v_list[0]].z;
                np_v0 << pts_list->points[nearest_edge_v_idx].x,
                        pts_list->points[nearest_edge_v_idx].y,
                        pts_list->points[nearest_edge_v_idx].z;
                np_v1 << pts_list->points[adj_edge_v_list_1[0]].x,
                        pts_list->points[adj_edge_v_list_1[0]].y,
                        pts_list->points[adj_edge_v_list_1[0]].z;
                np_v2 << pts_list->points[adj_edge_v_list_1[1]].x,
                        pts_list->points[adj_edge_v_list_1[1]].y,
                        pts_list->points[adj_edge_v_list_1[1]].z;
                cur_pv1 = p_v1 - cur_eigen_pt;
                npv0_npv1 = np_v1 - np_v0;
                npv0_npv2 = np_v2 - np_v0;

                // 添加非平面边缘点
                noplane_edge_v_list.push_back(nearest_edge_v_idx);
                double dot_1 = cur_pv1.dot(npv0_npv1);
                double dot_2 = cur_pv1.dot(npv0_npv2);

                // 无法区分同侧还是异侧，则跳过
                if (dot_1 * dot_2 > 0)
                {
                    break;
                }

                if (dot_1 >= 0)
                {
                    noplane_edge_v_list.push_back(adj_edge_v_list_1[0]);
                }
                else
                {
                    noplane_edge_v_list.push_back(adj_edge_v_list_1[1]);
                }

                // 添加平面边缘点
                plane_edge_v_list.push_back(cur_edge_v_idx);
                plane_edge_v_list.push_back(adj_edge_v_list[0]);

            }
            // 一侧或没有，选择一侧延展处理
            else if (pts_state[adj_edge_v_list[0]] == 5 && pts_state[adj_edge_v_list[1]] == 2)
            {
                // 添加平面边缘点
                plane_edge_v_list.push_back(cur_edge_v_idx);
                plane_edge_v_list.push_back(adj_edge_v_list[0]);

                // 添加非平面边缘点
                noplane_edge_v_list.push_back(adj_edge_v_list[1]);

                std::vector<int> adj_edge_v_list_2;
                find_edge_vertex(adj_edge_v_list[1], adj_edge_v_list_2);
                if (adj_edge_v_list_2.size() != 2)
                {
                    break;
                }
                // 判断是否是闭合三角形情况
                if (adj_edge_v_list_2[0] != cur_edge_v_idx)
                {
                    if (pts_state[adj_edge_v_list_2[0]] == 2)
                    {
                        noplane_edge_v_list.push_back(adj_edge_v_list_2[0]);
                    }
                    else
                    {
                        // 补洞
                        std::vector<int> edge_v;
                        bool if_hole = if_lay_in_hole(cur_edge_v_idx, adj_edge_v_list, edge_v);
                        if (if_hole)
                        {
                            fill_hole(edge_v);
                            for (auto temp_v_idx : edge_v)
                            {
                                if (!if_edge_vertex(temp_v_idx))
                                {
                                    if (pts_state[temp_v_idx] == 2)
                                    {
                                        pts_state[temp_v_idx] = 1;
                                    }

                                    if (pts_state[temp_v_idx] == 5)
                                    {
                                        pts_state[temp_v_idx] = 6;
                                    }
                                }
                            }
                        }
                        break;
                    }
                }
                else
                {
                    if (pts_state[adj_edge_v_list_2[1]] == 2)
                    {
                        noplane_edge_v_list.push_back(adj_edge_v_list_2[1]);
                    }
                    else
                    {
                        // 补洞
                        std::vector<int> edge_v;
                        bool if_hole = if_lay_in_hole(cur_edge_v_idx, adj_edge_v_list, edge_v);
                        if (if_hole)
                        {
                            fill_hole(edge_v);
                            for (auto temp_v_idx : edge_v)
                            {
                                if (!if_edge_vertex(temp_v_idx))
                                {
                                    if (pts_state[temp_v_idx] == 2)
                                    {
                                        pts_state[temp_v_idx] = 1;
                                    }

                                    if (pts_state[temp_v_idx] == 5)
                                    {
                                        pts_state[temp_v_idx] = 6;
                                    }
                                }
                            }
                        }
                        break;
                    }
                }
            }
            else if (pts_state[adj_edge_v_list[0]] == 2 && pts_state[adj_edge_v_list[1]] == 5)
            {
                // 添加平面边缘点
                plane_edge_v_list.push_back(cur_edge_v_idx);
                plane_edge_v_list.push_back(adj_edge_v_list[1]);

                // 添加非平面边缘点
                noplane_edge_v_list.push_back(adj_edge_v_list[0]);

                std::vector<int> adj_edge_v_list_2;
                find_edge_vertex(adj_edge_v_list[0], adj_edge_v_list_2);
                if (adj_edge_v_list_2.size() != 2)
                {
                    break;
                }
                // 判断是否是闭合三角形情况
                if (adj_edge_v_list_2[0] != cur_edge_v_idx)
                {
                    if (pts_state[adj_edge_v_list_2[0]] == 2)
                    {
                        noplane_edge_v_list.push_back(adj_edge_v_list_2[0]);
                    }
                    else
                    {
                        // 补洞
                        std::vector<int> edge_v;
                        bool if_hole = if_lay_in_hole(cur_edge_v_idx, adj_edge_v_list, edge_v);
                        if (if_hole)
                        {
                            fill_hole(edge_v);
                            for (auto temp_v_idx : edge_v)
                            {
                                if (!if_edge_vertex(temp_v_idx))
                                {
                                    if (pts_state[temp_v_idx] == 2)
                                    {
                                        pts_state[temp_v_idx] = 1;
                                    }

                                    if (pts_state[temp_v_idx] == 5)
                                    {
                                        pts_state[temp_v_idx] = 6;
                                    }
                                }
                            }
                        }
                        break;
                    }
                }
                else
                {
                    if (pts_state[adj_edge_v_list_2[1]] == 2)
                    {
                        noplane_edge_v_list.push_back(adj_edge_v_list_2[1]);
                    }
                    else
                    {
                        // 补洞
                        std::vector<int> edge_v;
                        bool if_hole = if_lay_in_hole(cur_edge_v_idx, adj_edge_v_list, edge_v);
                        if (if_hole)
                        {
                            fill_hole(edge_v);
                            for (auto temp_v_idx : edge_v)
                            {
                                if (!if_edge_vertex(temp_v_idx))
                                {
                                    if (pts_state[temp_v_idx] == 2)
                                    {
                                        pts_state[temp_v_idx] = 1;
                                    }

                                    if (pts_state[temp_v_idx] == 5)
                                    {
                                        pts_state[temp_v_idx] = 6;
                                    }
                                }
                            }
                        }
                        break;
                    }
                }
            }
            // 平面只有一个点，非平面不确定有几个点
            else if (pts_state[adj_edge_v_list[0]] == 2 && pts_state[adj_edge_v_list[1]] == 2)
            {
                // TODO 添加补洞操作
                // 添加平面边缘点
                plane_edge_v_list.push_back(cur_edge_v_idx);
                // 添加非平面边缘点
                noplane_edge_v_list.push_back(adj_edge_v_list[0]);
                noplane_edge_v_list.push_back(adj_edge_v_list[1]);

                if_need_extend = false;
                if_p1_npn = false;  // 保证非平面两个点

                // 补洞
                std::vector<int> edge_v;
                bool if_hole = if_lay_in_hole(cur_edge_v_idx, adj_edge_v_list, edge_v);
                if (if_hole)
                {
                    fill_hole(edge_v);
                    for (auto temp_v_idx : edge_v)
                    {
                        if (!if_edge_vertex(temp_v_idx))
                        {
                            if (pts_state[temp_v_idx] == 2)
                            {
                                pts_state[temp_v_idx] = 1;
                            }

                            if (pts_state[temp_v_idx] == 5)
                            {
                                pts_state[temp_v_idx] = 6;
                            }
                        }
                    }
                }
                break;
            }
            else
            {
                break;
            }

            if (if_need_extend)
            {
                // 判断1-1 or 1-n or n-1
                Eigen::Vector3d p_cur, p_next, np_cur, np_next;
                int p_cur_idx = plane_edge_v_list[plane_edge_v_list.size() - 1];
                int np_cur_idx = noplane_edge_v_list[noplane_edge_v_list.size() - 1];
                int p_prior_idx = plane_edge_v_list[plane_edge_v_list.size() - 2];
                int np_prior_idx = noplane_edge_v_list[noplane_edge_v_list.size() - 2];
                int p_next_idx, np_next_idx;
                p_cur << pts_list->points[p_cur_idx].x,
                        pts_list->points[p_cur_idx].y,
                        pts_list->points[p_cur_idx].z;
                np_cur << pts_list->points[np_cur_idx].x,
                        pts_list->points[np_cur_idx].y,
                        pts_list->points[np_cur_idx].z;

                // 如果平面、非平面的方向反向则跳出循环
                Eigen::Vector3d p_prior, np_prior, p_pri_cur, np_pri_cur;
                p_prior << pts_list->points[p_prior_idx].x,
                        pts_list->points[p_prior_idx].y,
                        pts_list->points[p_prior_idx].z;
                np_prior << pts_list->points[np_prior_idx].x,
                        pts_list->points[np_prior_idx].y,
                        pts_list->points[np_prior_idx].z;
                p_pri_cur = p_cur - p_prior;
                np_pri_cur = np_cur - np_prior;
                if (p_pri_cur.dot(np_pri_cur) < 0)
                {
                    break;
                }

                // 如果闭合一侧距离过远则跳出循环
                if ((p_prior - np_prior).norm() > r_dis * 2)
                {
                    break;
                }

                // 确定平面、非平面list最后一个点的下一个延展点
                std::vector<int> pcur_adj_v_list, npcur_adj_v_list;
                find_edge_vertex(p_cur_idx, pcur_adj_v_list);
                find_edge_vertex(np_cur_idx, npcur_adj_v_list);
                if (pcur_adj_v_list.size() != 2 || npcur_adj_v_list.size() != 2)
                {
                    break;
                }
                if (pcur_adj_v_list[0] != p_prior_idx)
                {
                    p_next_idx = pcur_adj_v_list[0];
                }
                else
                {
                    p_next_idx = pcur_adj_v_list[1];
                }
                if (npcur_adj_v_list[0] != np_prior_idx)
                {
                    np_next_idx = npcur_adj_v_list[0];
                }
                else
                {
                    np_next_idx = npcur_adj_v_list[1];
                }
                p_next << pts_list->points[p_next_idx].x,
                        pts_list->points[p_next_idx].y,
                        pts_list->points[p_next_idx].z;
                np_next << pts_list->points[np_next_idx].x,
                        pts_list->points[np_next_idx].y,
                        pts_list->points[np_next_idx].z;


                // 判断下一个延展点是否闭合
                if (pts_state[p_next_idx] == 2 && pts_state[np_next_idx] == 5)
                {
                    // 闭合四边形情况 -> 不需要扩展，p1_npn 2-2处理
                    if_need_extend = false;
                }
                else if (pts_state[p_next_idx] == 2 && pts_state[np_next_idx] == 2)
                {
                    // 闭合，p-np 2-n情况 -> 不需要扩展，补洞处理 TODO
                    if_need_extend = false;
                }
                else if (pts_state[p_next_idx] == 5 && pts_state[np_next_idx] == 5)
                {
                    // 闭合，p-np n-2情况 -> 不需要扩展，p1_npn 2-2处理
                    if_need_extend = false;
                }
                else
                {
                    // 正常情况
                }

                // 第一次扩展，判断扩展类型，1-1、1-n、n-1
                if (if_need_extend)
                {
                    // 计算三个距离，p_cur->np_cur、p_cur->np_next、p_next->np_cur
                    double dis_1 = (p_cur - np_cur).norm();
                    double dis_2 = (p_cur - np_next).norm();
                    double dis_3 = (p_next - np_cur).norm();

                    // p_cur->np_cur距离最短则无需扩展
                    if (dis_1 < dis_2 && dis_1 < dis_3)
                    {
                        if_need_extend = false;
                    }
                    else
                    {
                        // p-np 1-n
                        if (dis_2 <= dis_3)
                        {
                            if_p1_npn = true;
                            noplane_edge_v_list.push_back(np_next_idx);
                            np_prior_idx = np_cur_idx;
                            np_cur_idx = np_next_idx;
                        }
                            // p-np n-1
                        else
                        {
                            if_p1_npn = false;
                            plane_edge_v_list.push_back(p_next_idx);
                            p_prior_idx = p_cur_idx;
                            p_cur_idx = p_next_idx;
                        }
                    }
                }

                // 补充完整 plane_edge_v_list, noplane_edge_v_list
                // 根据扩展类型，循环直至满足 p_cur->np_cur 距离最短
                while (if_need_extend)
                {
                    // p-np 1-n
                    if (if_p1_npn)
                    {
                        std::vector<int> temp_npcur_adj_v_list;
                        find_edge_vertex(np_cur_idx, temp_npcur_adj_v_list);
                        if (temp_npcur_adj_v_list.size() != 2)
                        {
                            break;
                        }
                        if (temp_npcur_adj_v_list[0] != np_prior_idx)
                        {
                            np_next_idx = temp_npcur_adj_v_list[0];
                        }
                        else
                        {
                            np_next_idx = temp_npcur_adj_v_list[1];
                        }
                        np_next << pts_list->points[np_next_idx].x,
                                pts_list->points[np_next_idx].y,
                                pts_list->points[np_next_idx].z;
                        // 判断下一个延展点是否闭合  TODO 判断条件指向性太强，有无问题？ 会造成空洞？
                        if (pts_state[np_next_idx] != 2)
                        {
                            break;
                        }
                        // 计算三个距离，p_cur->np_cur、p_cur->np_next、p_next->np_cur
                        double dis_cur = (p_cur - np_cur).norm();
                        double dis_next = (p_cur - np_next).norm();

                        if (dis_cur >= dis_next)
                        {
                            auto temp_ys_idx = std::find(noplane_edge_v_list.begin(), noplane_edge_v_list.end(), np_next_idx);
                            if (temp_ys_idx == noplane_edge_v_list.end())
                            {
                                noplane_edge_v_list.push_back(np_next_idx);
                                np_prior_idx = np_cur_idx;
                                np_cur_idx = np_next_idx;
                            }
                            else
                            {
                                // 避免出现死循环
                                break;
                            }
                        }
                        else
                        {
                            break;
                        }
                    }
                    // p-np n-1
                    else
                    {
                        std::vector<int> temp_pcur_adj_v_list;
                        find_edge_vertex(p_cur_idx, temp_pcur_adj_v_list);
                        if (temp_pcur_adj_v_list.size() != 2)
                        {
                            break;
                        }
                        if (temp_pcur_adj_v_list[0] != p_prior_idx)
                        {
                            p_next_idx = temp_pcur_adj_v_list[0];
                        }
                        else
                        {
                            p_next_idx = temp_pcur_adj_v_list[1];
                        }
                        p_next << pts_list->points[p_next_idx].x,
                                pts_list->points[p_next_idx].y,
                                pts_list->points[p_next_idx].z;
                        // 判断下一个延展点是否闭合  TODO 判断条件指向性太强，有无问题？
                        if (pts_state[p_next_idx] != 5)
                        {
                            break;
                        }
                        // 计算三个距离，p_cur->np_cur、p_cur->np_next、p_next->np_cur
                        double dis_cur = (np_cur - p_cur).norm();
                        double dis_next = (np_cur - p_next).norm();

                        if (dis_cur >= dis_next)
                        {
                            auto temp_ys_idx = std::find(plane_edge_v_list.begin(), plane_edge_v_list.end(), p_next_idx);
                            if (temp_ys_idx == plane_edge_v_list.end())
                            {
                                plane_edge_v_list.push_back(p_next_idx);
                                p_prior_idx = p_cur_idx;
                                p_cur_idx = p_next_idx;
                            }
                            else
                            {
                                // 避免出现死循环
                                break;
                            }
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }

            // 构建三角剖分
            std::vector<int> added_tri_list;
            // p-np 1-n
            if (if_p1_npn)
            {
                // 三个点
                if (noplane_edge_v_list.size() == 1)
                {
                    // 三角形法向
                    int plane_v0_idx = plane_edge_v_list[0];
                    int adj_tri_idx = pt2face_idx[plane_v0_idx][0];
                    int tri_A_idx = faces_list_all[adj_tri_idx][0];
                    int tri_B_idx = faces_list_all[adj_tri_idx][1];
                    int tri_C_idx = faces_list_all[adj_tri_idx][2];
                    Eigen::Vector3d tri_A, tri_B, tri_C;
                    tri_A << pts_list->points[tri_A_idx].x, pts_list->points[tri_A_idx].y,
                            pts_list->points[tri_A_idx].z;
                    tri_B << pts_list->points[tri_B_idx].x, pts_list->points[tri_B_idx].y,
                            pts_list->points[tri_B_idx].z;
                    tri_C << pts_list->points[tri_C_idx].x, pts_list->points[tri_C_idx].y,
                            pts_list->points[tri_C_idx].z;
                    Eigen::Vector3d tri_ab = tri_B - tri_A;
                    Eigen::Vector3d tri_ac = tri_C - tri_A;
                    Eigen::Vector3d tri_n = tri_ab.cross(tri_ac);

                    // 添加三角形
                    int plane_v1_idx = plane_edge_v_list[1];
                    int noplane_v0_idx = noplane_edge_v_list[0];
                    Eigen::Vector3d plane_v0, plane_v1, noplane_v0;
                    plane_v0 << pts_list->points[plane_v0_idx].x,
                            pts_list->points[plane_v0_idx].y,
                            pts_list->points[plane_v0_idx].z;
                    plane_v1 << pts_list->points[plane_v1_idx].x,
                            pts_list->points[plane_v1_idx].y,
                            pts_list->points[plane_v1_idx].z;
                    noplane_v0 << pts_list->points[noplane_v0_idx].x,
                            pts_list->points[noplane_v0_idx].y,
                            pts_list->points[noplane_v0_idx].z;
                    Eigen::Vector3d planev0_planev1 = plane_v1 - plane_v0;
                    Eigen::Vector3d planev0_noplanev0 = noplane_v0 - plane_v0;
                    Eigen::Vector3d n = planev0_planev1.cross(planev0_noplanev0);
                    if (tri_n.dot(n) > 0)
                    {
                        std::array<int, 3> temp_add_t;
                        temp_add_t[0] = plane_v0_idx;
                        temp_add_t[1] = plane_v1_idx;
                        temp_add_t[2] = noplane_v0_idx;
                        added_tri_list.push_back(faces_list_all.size());
                        faces_list_all.push_back(temp_add_t);
                    }
                    else
                    {
                        // 添加三角形
                        std::array<int, 3> temp_add_t;
                        temp_add_t[0] = plane_v0_idx;
                        temp_add_t[1] = noplane_v0_idx;
                        temp_add_t[2] = plane_v1_idx;
                        added_tri_list.push_back(faces_list_all.size());
                        faces_list_all.push_back(temp_add_t);
                    }
                }
                // 四个点，四边形
                else if (noplane_edge_v_list.size() == 2)
                {
                    // 判断哪个对角线较短
                    int plane_v0_idx = plane_edge_v_list[0];
                    int plane_v1_idx = plane_edge_v_list[1];
                    int noplane_v0_idx = noplane_edge_v_list[0];
                    int noplane_v1_idx = noplane_edge_v_list[1];
                    Eigen::Vector3d plane_v0, plane_v1, noplane_v0, noplane_v1;
                    plane_v0 << pts_list->points[plane_v0_idx].x,
                            pts_list->points[plane_v0_idx].y,
                            pts_list->points[plane_v0_idx].z;
                    plane_v1 << pts_list->points[plane_v1_idx].x,
                            pts_list->points[plane_v1_idx].y,
                            pts_list->points[plane_v1_idx].z;
                    noplane_v0 << pts_list->points[noplane_v0_idx].x,
                            pts_list->points[noplane_v0_idx].y,
                            pts_list->points[noplane_v0_idx].z;
                    noplane_v1 << pts_list->points[noplane_v1_idx].x,
                            pts_list->points[noplane_v1_idx].y,
                            pts_list->points[noplane_v1_idx].z;
                    double dis1 = (plane_v0 - noplane_v1).norm();
                    double dis2 = (plane_v1 - noplane_v0).norm();

                    // 三角形法向
                    int adj_tri_idx = pt2face_idx[plane_v0_idx][0];
                    int tri_A_idx = faces_list_all[adj_tri_idx][0];
                    int tri_B_idx = faces_list_all[adj_tri_idx][1];
                    int tri_C_idx = faces_list_all[adj_tri_idx][2];
                    Eigen::Vector3d tri_A, tri_B, tri_C;
                    tri_A << pts_list->points[tri_A_idx].x, pts_list->points[tri_A_idx].y,
                            pts_list->points[tri_A_idx].z;
                    tri_B << pts_list->points[tri_B_idx].x, pts_list->points[tri_B_idx].y,
                            pts_list->points[tri_B_idx].z;
                    tri_C << pts_list->points[tri_C_idx].x, pts_list->points[tri_C_idx].y,
                            pts_list->points[tri_C_idx].z;
                    Eigen::Vector3d tri_ab = tri_B - tri_A;
                    Eigen::Vector3d tri_ac = tri_C - tri_A;
                    Eigen::Vector3d tri_n = tri_ab.cross(tri_ac);

                    // 添加三角形
                    Eigen::Vector3d planev0_planev1 = plane_v1 - plane_v0;
                    Eigen::Vector3d noplanev0_noplanev1 = noplane_v1 - noplane_v0;
                    Eigen::Vector3d planev0_noplanev1 = noplane_v1 - plane_v0;
                    Eigen::Vector3d noplanev0_planev1 = plane_v1 - noplane_v0;

                    if (dis1 < dis2)
                    {
                        // 第一个三角形
                        Eigen::Vector3d n = planev0_planev1.cross(planev0_noplanev1);
                        if (tri_n.dot(n) > 0)
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = plane_v0_idx;
                            temp_add_t[1] = plane_v1_idx;
                            temp_add_t[2] = noplane_v1_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                        else
                        {
                            // 添加三角形
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = plane_v0_idx;
                            temp_add_t[1] = noplane_v1_idx;
                            temp_add_t[2] = plane_v1_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                        // 第二个三角形
                        n = noplanev0_noplanev1.cross(planev0_noplanev1);
                        if (tri_n.dot(n) > 0)
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = noplane_v1_idx;
                            temp_add_t[1] = noplane_v0_idx;
                            temp_add_t[2] = plane_v0_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                        else
                        {
                            // 添加三角形
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = noplane_v1_idx;
                            temp_add_t[1] = plane_v0_idx;
                            temp_add_t[2] = noplane_v0_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                    }
                    else
                    {
                        // 第一个三角形
                        Eigen::Vector3d n = planev0_planev1.cross(noplanev0_planev1);
                        if (tri_n.dot(n) > 0)
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = plane_v1_idx;
                            temp_add_t[1] = plane_v0_idx;
                            temp_add_t[2] = noplane_v0_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                        else
                        {
                            // 添加三角形
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = plane_v1_idx;
                            temp_add_t[1] = noplane_v0_idx;
                            temp_add_t[2] = plane_v0_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                        // 第二个三角形
                        n = noplanev0_noplanev1.cross(noplanev0_planev1);
                        if (tri_n.dot(n) > 0)
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = noplane_v0_idx;
                            temp_add_t[1] = noplane_v1_idx;
                            temp_add_t[2] = plane_v1_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                        else
                        {
                            // 添加三角形
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = noplane_v0_idx;
                            temp_add_t[1] = plane_v1_idx;
                            temp_add_t[2] = noplane_v1_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                    }
                }
                else
                {
                    int medium_idx = noplane_edge_v_list.size() / 2;

                    // 确定三角形法向
                    int plane_v0_idx = plane_edge_v_list[0];
                    int plane_v1_idx = plane_edge_v_list[1];
                    int adj_tri_idx = pt2face_idx[plane_v0_idx][0];
                    int tri_A_idx = faces_list_all[adj_tri_idx][0];
                    int tri_B_idx = faces_list_all[adj_tri_idx][1];
                    int tri_C_idx = faces_list_all[adj_tri_idx][2];
                    Eigen::Vector3d tri_A, tri_B, tri_C;
                    tri_A << pts_list->points[tri_A_idx].x, pts_list->points[tri_A_idx].y,
                            pts_list->points[tri_A_idx].z;
                    tri_B << pts_list->points[tri_B_idx].x, pts_list->points[tri_B_idx].y,
                            pts_list->points[tri_B_idx].z;
                    tri_C << pts_list->points[tri_C_idx].x, pts_list->points[tri_C_idx].y,
                            pts_list->points[tri_C_idx].z;
                    Eigen::Vector3d tri_ab = tri_B - tri_A;
                    Eigen::Vector3d tri_ac = tri_C - tri_A;
                    Eigen::Vector3d tri_n = tri_ab.cross(tri_ac);

                    // 添加位于中间的三角形
                    int noplane_vm_idx = noplane_edge_v_list[medium_idx];
                    Eigen::Vector3d plane_v0, plane_v1, noplane_vm;
                    plane_v0 << pts_list->points[plane_v0_idx].x,
                            pts_list->points[plane_v0_idx].y,
                            pts_list->points[plane_v0_idx].z;
                    plane_v1 << pts_list->points[plane_v1_idx].x,
                            pts_list->points[plane_v1_idx].y,
                            pts_list->points[plane_v1_idx].z;
                    noplane_vm << pts_list->points[noplane_vm_idx].x,
                            pts_list->points[noplane_vm_idx].y,
                            pts_list->points[noplane_vm_idx].z;
                    Eigen::Vector3d planev0_planev1 = plane_v1 - plane_v0;
                    Eigen::Vector3d planev0_noplanevm = noplane_vm - plane_v0;
                    Eigen::Vector3d n = planev0_planev1.cross(planev0_noplanevm);
                    if (tri_n.dot(n) > 0)
                    {
                        std::array<int, 3> temp_add_t;
                        temp_add_t[0] = plane_v0_idx;
                        temp_add_t[1] = plane_v1_idx;
                        temp_add_t[2] = noplane_vm_idx;
                        added_tri_list.push_back(faces_list_all.size());
                        faces_list_all.push_back(temp_add_t);
                    }
                    else
                    {
                        std::array<int, 3> temp_add_t;
                        temp_add_t[0] = plane_v0_idx;
                        temp_add_t[1] = noplane_vm_idx;
                        temp_add_t[2] = plane_v1_idx;
                        added_tri_list.push_back(faces_list_all.size());
                        faces_list_all.push_back(temp_add_t);
                    }

                    // 根据点多的一侧的中点分开，前面、中间、后面分别添加
                    for (int j = 0; j < noplane_edge_v_list.size()-1; j++)
                    {
                        int A_idx = noplane_edge_v_list[j];
                        int B_idx = noplane_edge_v_list[j+1];
                        Eigen::Vector3d A, B, n;
                        A << pts_list->points[A_idx].x, pts_list->points[A_idx].y, pts_list->points[A_idx].z;
                        B << pts_list->points[B_idx].x, pts_list->points[B_idx].y, pts_list->points[B_idx].z;

                        if (j < medium_idx)
                        {
                            Eigen::Vector3d planev0_A = A - plane_v0;
                            Eigen::Vector3d planev0_B = B - plane_v0;
                            n = planev0_A.cross(planev0_B);
                            if (tri_n.dot(n) > 0)
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = plane_v0_idx;
                                temp_add_t[1] = A_idx;
                                temp_add_t[2] = B_idx;
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                            }
                            else
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = plane_v0_idx;
                                temp_add_t[1] = B_idx;
                                temp_add_t[2] = A_idx;
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                            }
                        }
                        else
                        {
                            Eigen::Vector3d planev1_A = A - plane_v1;
                            Eigen::Vector3d planev1_B = B - plane_v1;
                            n = planev1_A.cross(planev1_B);
                            if (tri_n.dot(n) > 0)
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = plane_v1_idx;
                                temp_add_t[1] = A_idx;
                                temp_add_t[2] = B_idx;
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                            }
                            else
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = plane_v1_idx;
                                temp_add_t[1] = B_idx;
                                temp_add_t[2] = A_idx;
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                            }
                        }
                    }

                }
            }
            // p-np n-1
            else
            {
                // 三个点
                if (plane_edge_v_list.size() == 1)
                {
                    // 三角形法向
                    int plane_v0_idx = plane_edge_v_list[0];
                    int adj_tri_idx = pt2face_idx[plane_v0_idx][0];
                    int tri_A_idx = faces_list_all[adj_tri_idx][0];
                    int tri_B_idx = faces_list_all[adj_tri_idx][1];
                    int tri_C_idx = faces_list_all[adj_tri_idx][2];
                    Eigen::Vector3d tri_A, tri_B, tri_C;
                    tri_A << pts_list->points[tri_A_idx].x, pts_list->points[tri_A_idx].y,
                            pts_list->points[tri_A_idx].z;
                    tri_B << pts_list->points[tri_B_idx].x, pts_list->points[tri_B_idx].y,
                            pts_list->points[tri_B_idx].z;
                    tri_C << pts_list->points[tri_C_idx].x, pts_list->points[tri_C_idx].y,
                            pts_list->points[tri_C_idx].z;
                    Eigen::Vector3d tri_ab = tri_B - tri_A;
                    Eigen::Vector3d tri_ac = tri_C - tri_A;
                    Eigen::Vector3d tri_n = tri_ab.cross(tri_ac);

                    // 添加三角形
                    int noplane_v0_idx = noplane_edge_v_list[0];
                    int noplane_v1_idx = noplane_edge_v_list[1];
                    Eigen::Vector3d plane_v0, noplane_v1, noplane_v0;
                    plane_v0 << pts_list->points[plane_v0_idx].x,
                            pts_list->points[plane_v0_idx].y,
                            pts_list->points[plane_v0_idx].z;
                    noplane_v0 << pts_list->points[noplane_v0_idx].x,
                            pts_list->points[noplane_v0_idx].y,
                            pts_list->points[noplane_v0_idx].z;
                    noplane_v1 << pts_list->points[noplane_v1_idx].x,
                            pts_list->points[noplane_v1_idx].y,
                            pts_list->points[noplane_v1_idx].z;
                    Eigen::Vector3d planev0_noplanev0 = noplane_v0 - plane_v0;
                    Eigen::Vector3d planev0_noplanev1 = noplane_v1 - plane_v0;
                    Eigen::Vector3d n = planev0_noplanev0.cross(planev0_noplanev1);
                    if (tri_n.dot(n) > 0)
                    {
                        std::array<int, 3> temp_add_t;
                        temp_add_t[0] = plane_v0_idx;
                        temp_add_t[1] = noplane_v0_idx;
                        temp_add_t[2] = noplane_v1_idx;
                        added_tri_list.push_back(faces_list_all.size());
                        faces_list_all.push_back(temp_add_t);
                    }
                    else
                    {
                        // 添加三角形
                        std::array<int, 3> temp_add_t;
                        temp_add_t[0] = plane_v0_idx;
                        temp_add_t[1] = noplane_v1_idx;
                        temp_add_t[2] = noplane_v0_idx;
                        added_tri_list.push_back(faces_list_all.size());
                        faces_list_all.push_back(temp_add_t);
                    }
                }
                // 四个点，四边形
                else if (plane_edge_v_list.size() == 2)
                {
                    // 判断哪个对角线较短
                    int plane_v0_idx = plane_edge_v_list[0];
                    int plane_v1_idx = plane_edge_v_list[1];
                    int noplane_v0_idx = noplane_edge_v_list[0];
                    int noplane_v1_idx = noplane_edge_v_list[1];
                    Eigen::Vector3d plane_v0, plane_v1, noplane_v0, noplane_v1;
                    plane_v0 << pts_list->points[plane_v0_idx].x,
                            pts_list->points[plane_v0_idx].y,
                            pts_list->points[plane_v0_idx].z;
                    plane_v1 << pts_list->points[plane_v1_idx].x,
                            pts_list->points[plane_v1_idx].y,
                            pts_list->points[plane_v1_idx].z;
                    noplane_v0 << pts_list->points[noplane_v0_idx].x,
                            pts_list->points[noplane_v0_idx].y,
                            pts_list->points[noplane_v0_idx].z;
                    noplane_v1 << pts_list->points[noplane_v1_idx].x,
                            pts_list->points[noplane_v1_idx].y,
                            pts_list->points[noplane_v1_idx].z;
                    double dis1 = (plane_v0 - noplane_v1).norm();
                    double dis2 = (plane_v1 - noplane_v0).norm();

                    // 三角形法向
                    int adj_tri_idx = pt2face_idx[plane_v0_idx][0];
                    int tri_A_idx = faces_list_all[adj_tri_idx][0];
                    int tri_B_idx = faces_list_all[adj_tri_idx][1];
                    int tri_C_idx = faces_list_all[adj_tri_idx][2];
                    Eigen::Vector3d tri_A, tri_B, tri_C;
                    tri_A << pts_list->points[tri_A_idx].x, pts_list->points[tri_A_idx].y,
                            pts_list->points[tri_A_idx].z;
                    tri_B << pts_list->points[tri_B_idx].x, pts_list->points[tri_B_idx].y,
                            pts_list->points[tri_B_idx].z;
                    tri_C << pts_list->points[tri_C_idx].x, pts_list->points[tri_C_idx].y,
                            pts_list->points[tri_C_idx].z;
                    Eigen::Vector3d tri_ab = tri_B - tri_A;
                    Eigen::Vector3d tri_ac = tri_C - tri_A;
                    Eigen::Vector3d tri_n = tri_ab.cross(tri_ac);

                    // 添加三角形
                    Eigen::Vector3d planev0_planev1 = plane_v1 - plane_v0;
                    Eigen::Vector3d noplanev0_noplanev1 = noplane_v1 - noplane_v0;
                    Eigen::Vector3d planev0_noplanev1 = noplane_v1 - plane_v0;
                    Eigen::Vector3d noplanev0_planev1 = plane_v1 - noplane_v0;

                    if (dis1 < dis2)
                    {
                        // 第一个三角形
                        Eigen::Vector3d n = planev0_planev1.cross(planev0_noplanev1);
                        if (tri_n.dot(n) > 0)
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = plane_v0_idx;
                            temp_add_t[1] = plane_v1_idx;
                            temp_add_t[2] = noplane_v1_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                        else
                        {
                            // 添加三角形
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = plane_v0_idx;
                            temp_add_t[1] = noplane_v1_idx;
                            temp_add_t[2] = plane_v1_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                        // 第二个三角形
                        n = noplanev0_noplanev1.cross(planev0_noplanev1);
                        if (tri_n.dot(n) > 0)
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = noplane_v1_idx;
                            temp_add_t[1] = noplane_v0_idx;
                            temp_add_t[2] = plane_v0_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                        else
                        {
                            // 添加三角形
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = noplane_v1_idx;
                            temp_add_t[1] = plane_v0_idx;
                            temp_add_t[2] = noplane_v0_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                    }
                    else
                    {
                        // 第一个三角形
                        Eigen::Vector3d n = planev0_planev1.cross(noplanev0_planev1);
                        if (tri_n.dot(n) > 0)
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = plane_v1_idx;
                            temp_add_t[1] = plane_v0_idx;
                            temp_add_t[2] = noplane_v0_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                        else
                        {
                            // 添加三角形
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = plane_v1_idx;
                            temp_add_t[1] = noplane_v0_idx;
                            temp_add_t[2] = plane_v0_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                        // 第二个三角形
                        n = noplanev0_noplanev1.cross(noplanev0_planev1);
                        if (tri_n.dot(n) > 0)
                        {
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = noplane_v0_idx;
                            temp_add_t[1] = noplane_v1_idx;
                            temp_add_t[2] = plane_v1_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                        else
                        {
                            // 添加三角形
                            std::array<int, 3> temp_add_t;
                            temp_add_t[0] = noplane_v0_idx;
                            temp_add_t[1] = plane_v1_idx;
                            temp_add_t[2] = noplane_v1_idx;
                            added_tri_list.push_back(faces_list_all.size());
                            faces_list_all.push_back(temp_add_t);
                        }
                    }
                }
                else
                {
                    int medium_idx = plane_edge_v_list.size() / 2;

                    // 确定三角形法向
                    int plane_v0_idx = plane_edge_v_list[0];
                    int adj_tri_idx = pt2face_idx[plane_v0_idx][0];
                    int tri_A_idx = faces_list_all[adj_tri_idx][0];
                    int tri_B_idx = faces_list_all[adj_tri_idx][1];
                    int tri_C_idx = faces_list_all[adj_tri_idx][2];
                    Eigen::Vector3d tri_A, tri_B, tri_C;
                    tri_A << pts_list->points[tri_A_idx].x, pts_list->points[tri_A_idx].y,
                            pts_list->points[tri_A_idx].z;
                    tri_B << pts_list->points[tri_B_idx].x, pts_list->points[tri_B_idx].y,
                            pts_list->points[tri_B_idx].z;
                    tri_C << pts_list->points[tri_C_idx].x, pts_list->points[tri_C_idx].y,
                            pts_list->points[tri_C_idx].z;
                    Eigen::Vector3d tri_ab = tri_B - tri_A;
                    Eigen::Vector3d tri_ac = tri_C - tri_A;
                    Eigen::Vector3d tri_n = tri_ab.cross(tri_ac);

                    // 添加位于中间的三角形
                    int noplane_v0_idx = noplane_edge_v_list[0];
                    int noplane_v1_idx = noplane_edge_v_list[1];
                    int plane_vm_idx = plane_edge_v_list[medium_idx];
                    Eigen::Vector3d noplane_v0, noplane_v1, plane_vm;
                    noplane_v0 << pts_list->points[noplane_v0_idx].x,
                            pts_list->points[noplane_v0_idx].y,
                            pts_list->points[noplane_v0_idx].z;
                    noplane_v1 << pts_list->points[noplane_v1_idx].x,
                            pts_list->points[noplane_v1_idx].y,
                            pts_list->points[noplane_v1_idx].z;
                    plane_vm << pts_list->points[plane_vm_idx].x,
                            pts_list->points[plane_vm_idx].y,
                            pts_list->points[plane_vm_idx].z;
                    Eigen::Vector3d noplanev0_noplanev1 = noplane_v1 - noplane_v0;
                    Eigen::Vector3d noplanev0_planevm = plane_vm - noplane_v0;
                    Eigen::Vector3d n = noplanev0_noplanev1.cross(noplanev0_planevm);
                    if (tri_n.dot(n) > 0)
                    {
                        std::array<int, 3> temp_add_t;
                        temp_add_t[0] = noplane_v0_idx;
                        temp_add_t[1] = noplane_v1_idx;
                        temp_add_t[2] = plane_vm_idx;
                        added_tri_list.push_back(faces_list_all.size());
                        faces_list_all.push_back(temp_add_t);
                    }
                    else
                    {
                        std::array<int, 3> temp_add_t;
                        temp_add_t[0] = noplane_v0_idx;
                        temp_add_t[1] = plane_vm_idx;
                        temp_add_t[2] = noplane_v1_idx;
                        added_tri_list.push_back(faces_list_all.size());
                        faces_list_all.push_back(temp_add_t);
                    }

                    // 根据点多的一侧的中点分开，前面、中间、后面分别添加
                    for (int j = 0; j < plane_edge_v_list.size()-1; j++)
                    {
                        int A_idx = plane_edge_v_list[j];
                        int B_idx = plane_edge_v_list[j+1];
                        Eigen::Vector3d A, B, n;
                        A << pts_list->points[A_idx].x, pts_list->points[A_idx].y, pts_list->points[A_idx].z;
                        B << pts_list->points[B_idx].x, pts_list->points[B_idx].y, pts_list->points[B_idx].z;

                        if (j < medium_idx)
                        {
                            Eigen::Vector3d noplanev0_A = A - noplane_v0;
                            Eigen::Vector3d noplanev0_B = B - noplane_v0;
                            n = noplanev0_A.cross(noplanev0_B);
                            if (tri_n.dot(n) > 0)
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = noplane_v0_idx;
                                temp_add_t[1] = A_idx;
                                temp_add_t[2] = B_idx;
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                            }
                            else
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = noplane_v0_idx;
                                temp_add_t[1] = B_idx;
                                temp_add_t[2] = A_idx;
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                            }
                        }
                        else
                        {
                            Eigen::Vector3d noplanev1_A = A - noplane_v1;
                            Eigen::Vector3d noplanev1_B = B - noplane_v1;
                            n = noplanev1_A.cross(noplanev1_B);
                            if (tri_n.dot(n) > 0)
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = noplane_v1_idx;
                                temp_add_t[1] = A_idx;
                                temp_add_t[2] = B_idx;
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                            }
                            else
                            {
                                std::array<int, 3> temp_add_t;
                                temp_add_t[0] = noplane_v1_idx;
                                temp_add_t[1] = B_idx;
                                temp_add_t[2] = A_idx;
                                added_tri_list.push_back(faces_list_all.size());
                                faces_list_all.push_back(temp_add_t);
                            }
                        }
                    }

                }
            }

            // 更新点-三角形索引列表
            for (int j = 0; j < added_tri_list.size(); j++)
            {
                for (int l = 0; l < 3; l++)
                {
                    int temp_id = faces_list_all[added_tri_list[j]][l];
                    pt2face_idx[temp_id].push_back(added_tri_list[j]);
                }
            }
            // 更新平面点的状态
            for (int j = 0; j < plane_edge_v_list.size(); j++)
            {
                int temp_id = plane_edge_v_list[j];
                if (!if_edge_vertex(temp_id))
                {
                    pts_state[temp_id] = 6;
                }
            }
            // 更新非平面点的状态
            for (int j = 0; j < noplane_edge_v_list.size(); j++)
            {
                int temp_id = noplane_edge_v_list[j];
                if (!if_edge_vertex(temp_id))
                {
                    pts_state[temp_id] = 1;
                }
            }

            // 判断该点是否需要处理另一侧
            if (if_passable && added_tri_list.size() > 0)
            {
                if_two_side = true;
            }

            // 继续延展直到遇到闭合点
            if (pts_state[plane_edge_v_list[plane_edge_v_list.size() - 1]] == 6)
            {
                break;
            }

            cur_edge_v_idx = plane_edge_v_list[plane_edge_v_list.size() - 1];
        }

        if (if_two_side)
        {
            i--;
        }
    }
}

void NonPlaneMesh::save_to_ply(const std::string &fileName)
{
    pcl::PolygonMesh mesh;
    mesh.header = pcl::PCLHeader();
    mesh.cloud = pcl::PCLPointCloud2();
    mesh.polygons = std::vector<pcl::Vertices>();

    pcl::toPCLPointCloud2(*pts_list, mesh.cloud);


    std::sort(faces_to_delete.begin(), faces_to_delete.end());   // 从小到大排序
    int delete_i = 0;
    for (int i = 0; i < faces_list_all.size(); i++)
    {
        if (faces_to_delete.size() > 0)
        {
            if (i == faces_to_delete[delete_i])
            {
                while (delete_i != faces_to_delete.size() - 1 && faces_to_delete[delete_i] <= i)
                {
                    delete_i++;
                }
                continue;
            }
        }

        pcl::Vertices vertices;
        vertices.vertices.push_back(faces_list_all[i][0]);   // meshlab点序逆时针
        vertices.vertices.push_back(faces_list_all[i][1]);
        vertices.vertices.push_back(faces_list_all[i][2]);
        mesh.polygons.push_back(vertices);
    }

//    pcl::io::savePLYFile(fileName, mesh);
    pcl::io::savePLYFileBinary(fileName, mesh);

    std::cout << "Plane mesh saved to " << fileName << std::endl;
    int num_all = faces_list_all.size();
    int num_delete = faces_to_delete.size();
    int num = faces_list_all.size() - faces_to_delete.size();
    std::cout << "The number of all triangulation is " << num_all << std::endl;
    std::cout << "The number of delete triangulation is " << num_delete << std::endl;
    std::cout << "The number of triangulation is " << num << std::endl;
}

bool NonPlaneMesh::save_to_ply_without_redundancy(const std::string &fileName)
{
    pcl::PolygonMesh mesh;
    mesh.header = pcl::PCLHeader();
    mesh.cloud = pcl::PCLPointCloud2();
    mesh.polygons = std::vector<pcl::Vertices>();

    std::vector<int> face_exist;
    std::sort(faces_to_delete.begin(), faces_to_delete.end());   // 从小到大排序
    int delete_i = 0;
    for (int i = 0; i < faces_list_all.size(); i++)
    {
        if (faces_to_delete.size() > 0)
        {
            if (i == faces_to_delete[delete_i])
            {
                while (delete_i != faces_to_delete.size() - 1 && faces_to_delete[delete_i] <= i)
                {
                    delete_i++;
                }
                continue;
            }
        }

        face_exist.push_back(i);
    }


    std::set<int> vertex_exist_single;
    std::unordered_map<int, int> vertex_idx_old2new;
    PointCloudXYZ::Ptr ptcl_processed = boost::make_shared<PointCloudXYZ>();

    for (auto it : face_exist)
    {
        vertex_exist_single.insert(faces_list_all[it][0]);
        vertex_exist_single.insert(faces_list_all[it][1]);
        vertex_exist_single.insert(faces_list_all[it][2]);
    }

    int i_ = 0;
    for (auto it : vertex_exist_single)
    {
        pcl::PointXYZ temp_pt;
        temp_pt.x = pts_list->points[it].x;
        temp_pt.y = pts_list->points[it].y;
        temp_pt.z = pts_list->points[it].z;

        ptcl_processed->push_back(temp_pt);
        vertex_idx_old2new[it] = i_;
        i_++;
    }
    pcl::toPCLPointCloud2(*ptcl_processed, mesh.cloud);

    for (auto it : face_exist)
    {
        int x_ = vertex_idx_old2new[faces_list_all[it][0]];
        int y_ = vertex_idx_old2new[faces_list_all[it][1]];
        int z_ = vertex_idx_old2new[faces_list_all[it][2]];

        pcl::Vertices vertices;
        vertices.vertices.push_back(x_);   // meshlab点序逆时针
        vertices.vertices.push_back(y_);
        vertices.vertices.push_back(z_);
        mesh.polygons.push_back(vertices);
    }


//    pcl::io::savePLYFile(fileName, mesh);
    pcl::io::savePLYFileBinary(fileName, mesh);

    std::cout << "Non-Planar mesh saved to " << fileName << std::endl;
    int num_all = faces_list_all.size();
    int num_delete = faces_to_delete.size();
    int num = faces_list_all.size() - faces_to_delete.size();
    std::cout << "The number of all triangulation is " << num_all << std::endl;
    std::cout << "The number of delete triangulation is " << num_delete << std::endl;
    std::cout << "The number of triangulation is " << num << std::endl;

    return true;
}

void NonPlaneMesh::plane_vertex_adjust(const std::shared_ptr<MeshFragment>& plane)
{
    for (auto it : plane->plane_vidx_to_noplnae_vidx)
    {
        pts_list->points[it.second].x = plane->ptcl_grid->points[it.first].x;
        pts_list->points[it.second].y = plane->ptcl_grid->points[it.first].y;
        pts_list->points[it.second].z = plane->ptcl_grid->points[it.first].z;
    }
}















