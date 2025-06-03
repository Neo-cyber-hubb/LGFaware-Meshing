//
// Created by neo on 2024/2/29.
//

#include "MeshFragment.h"


void MeshFragment::give_plane_point(const Eigen::Tensor<double, 3>& project_image,
                                      const std::vector<Eigen::Vector2i>& point_uv,
                                      Eigen::Matrix3d rot_mat,
                                      Eigen::Vector3d pos_vec)
{
    ptcl_all->clear();
    for ( int i = 0; i < point_uv.size(); i++ )
    {
        Eigen::Vector3d p_body;
        p_body(0) = project_image(point_uv[i](0), point_uv[i](1), 0);
        p_body(1) = project_image(point_uv[i](0), point_uv[i](1), 1);
        p_body(2) = project_image(point_uv[i](0), point_uv[i](1), 2);
        Eigen::Vector3d p_global(rot_mat * (p_body) + pos_vec);
        PointType temp_pt;
        temp_pt.x = p_global(0);
        temp_pt.y = p_global(1);
        temp_pt.z = p_global(2);
        temp_pt.intensity = project_image(point_uv[i](0), point_uv[i](1), 10);
        ptcl_all->push_back(temp_pt);
    }
}

void MeshFragment::compute_plane_parameter(const Eigen::Vector3d& frame_pos)
{
    Eigen::MatrixXd eigen_cloud(ptcl_all->size(), 3);
    for (int i = 0; i < ptcl_all->size(); ++i)
    {
        eigen_cloud(i, 0) = ptcl_all->points[i].x;
        eigen_cloud(i, 1) = ptcl_all->points[i].y;
        eigen_cloud(i, 2) = ptcl_all->points[i].z;
    }

    center_point = eigen_cloud.colwise().mean();
    Eigen::Vector3d lidar_beam = center_point - frame_pos;
    Eigen::Matrix3d cov_matrix = (eigen_cloud.rowwise() - center_point.transpose()).transpose() *
            (eigen_cloud.rowwise() - center_point.transpose()) / eigen_cloud.rows();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver( cov_matrix );
    Eigen::Matrix3d eigenvectors = eigensolver.eigenvectors();
    Eigen::Vector3d eigenvalues = eigensolver.eigenvalues();

//    quadtree_axes = eigenvectors;
    plane_quality = eigenvalues(0);      // TODO 是否要除以三个特征值的和
    normal_vector = eigenvectors.col(0);
    if (normal_vector.dot(lidar_beam) > 0)
    {
        normal_vector *= -1;
    }

    rec_index.center = center_point;
    for (int i = 0; i < 3; ++i)
    {
        rec_index.max_pt(i) = eigen_cloud.col(i).maxCoeff();
        rec_index.min_pt(i) = eigen_cloud.col(i).minCoeff();
    }
}

void MeshFragment::update_plane(const std::shared_ptr<MeshFragment> new_plane)
{
    // TODO
    int ptsNum_1 = ptcl_all->size();
    int ptsNum_2 = new_plane->ptcl_all->size();
    *ptcl_all += *(new_plane->ptcl_all);

    if (ptcl_all->size() <= ptsNum_thres)    // TODO 阈值 有必要否？
    {
        Eigen::MatrixXd eigen_cloud(ptcl_all->size(), 3);
        for (int i = 0; i < ptcl_all->size(); ++i)
        {
            eigen_cloud(i, 0) = ptcl_all->points[i].x;
            eigen_cloud(i, 1) = ptcl_all->points[i].y;
            eigen_cloud(i, 2) = ptcl_all->points[i].z;
        }

        center_point = eigen_cloud.colwise().mean();
        Eigen::Matrix3d cov_matrix = (eigen_cloud.rowwise() - center_point.transpose()).transpose() *
                                     (eigen_cloud.rowwise() - center_point.transpose()) / eigen_cloud.rows();

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver( cov_matrix );
        Eigen::Matrix3d eigenvectors = eigensolver.eigenvectors();
//        quadtree_axes = eigenvectors;

        if (normal_vector.dot(eigenvectors.col(0)) > 0)
        {
            normal_vector = eigenvectors.col(0);
        }
        else
        {
            normal_vector = -1 * eigenvectors.col(0);
        }

        rec_index.center = center_point;
        for (int i = 0; i < 3; ++i)
        {
            rec_index.max_pt(i) = eigen_cloud.col(i).maxCoeff();
            rec_index.min_pt(i) = eigen_cloud.col(i).minCoeff();
        }
    }
    else
    {
        Eigen::Vector3d new_center_point, new_normal_vector;
        AABB new_rec_index;
        new_center_point = (ptsNum_1 * center_point + ptsNum_2 * new_plane->center_point) / (ptsNum_1 + ptsNum_2);
        new_normal_vector = (ptsNum_1 * normal_vector + ptsNum_2 * new_plane->normal_vector) / (ptsNum_1 + ptsNum_2);
        for (int i = 0; i < 3; ++i)
        {
            new_rec_index.max_pt(i) = std::max(rec_index.max_pt(i), new_plane->rec_index.max_pt(i));
            new_rec_index.min_pt(i) = std::min(rec_index.min_pt(i), new_plane->rec_index.min_pt(i));
        }
        center_point = new_center_point;
        normal_vector = new_normal_vector;
        rec_index = new_rec_index;
        rec_index.center = center_point;
    }
    plane_update_num++;
}

// 初始化网格
grid MeshFragment::build_new_grid(double y_start, double z_start, int i, int j)
{
    grid temp_grid;
    temp_grid.pts_num = 0;
    temp_grid.node_type = 0;
    temp_grid.x_avg = 0;
    temp_grid.if_grid_need_update = false;
    temp_grid.if_mesh_need_delete = false;
    temp_grid.if_mesh_need_add = false;
    temp_grid.grid_vertex_state << 1, 1, 1, 1;
    temp_grid.grid_center << y_start + (j + 0.5) * dis_resolution,
            z_start - (i + 0.5) * dis_resolution;
    temp_grid.grid_vertex_list << y_start + j * dis_resolution,
            z_start - i * dis_resolution,
            y_start + (j + 1) * dis_resolution,
            z_start - i * dis_resolution,
            y_start + (j + 1) * dis_resolution,
            z_start - (i + 1) * dis_resolution,
            y_start + j * dis_resolution,
            z_start - (i + 1) * dis_resolution;
    temp_grid.edge_vertex << y_start + j * dis_resolution,
            y_start + (j + 1) * dis_resolution,
            z_start - i * dis_resolution,
            z_start - (i + 1) * dis_resolution,
            y_start + (j + 1) * dis_resolution,
            y_start + j * dis_resolution,
            z_start - (i + 1) * dis_resolution,
            z_start - i * dis_resolution;
    temp_grid.cell_level = 1;
    return temp_grid;
}

void MeshFragment::grid_decimate(double dis_resol, FILE *fpTime_mesh)
{
    // 第一次降采样，确定投影中心和投影矩阵
    dis_resolution = dis_resol;
    grid_need_update_list.clear();
    if (ptcl_all->size() > ptsNum_thres)   // 只有点的数量大于阈值才执行降采样和构建mesh
    {
        // 该平面第一次处理
        if ( if_first_decimate )
        {
            // 所有平面点投影到二维平面，按特征向量定义的坐标变换
            Eigen::MatrixXd eigen_cloud(ptcl_all->size(), 3);
            for (int i = 0; i < ptcl_all->size(); ++i)
            {
                eigen_cloud(i, 0) = ptcl_all->points[i].x;
                eigen_cloud(i, 1) = ptcl_all->points[i].y;
                eigen_cloud(i, 2) = ptcl_all->points[i].z;
            }
            quadtree_center = eigen_cloud.colwise().mean();
            Eigen::Matrix3d cov_matrix = (eigen_cloud.rowwise() - quadtree_center.transpose()).transpose() *
                                         (eigen_cloud.rowwise() - quadtree_center.transpose()) / eigen_cloud.rows();

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver( cov_matrix );

            // 特征向量的方向是随机的
            // 只需保证投影平面yz轴与法向符合右手笛卡尔坐标系，不需要管quadtree_axes.col(0)的方向
            quadtree_axes = eigensolver.eigenvectors();
            Eigen::Vector3d cross_product_yz = quadtree_axes.col(1).cross(quadtree_axes.col(2));  // y叉乘z
            if (cross_product_yz.dot(normal_vector) < 0)
            {
                quadtree_axes.col(1) *= -1;
            }

            Eigen::MatrixXd projected_pts(ptcl_all->size(), 3);
            projected_pts = (eigen_cloud.rowwise() - quadtree_center.transpose()) * quadtree_axes;
            yMax = projected_pts.col(1).maxCoeff();
            yMin = projected_pts.col(1).minCoeff();
            zMax = projected_pts.col(2).maxCoeff();
            zMin = projected_pts.col(2).minCoeff();


            y0 = (yMin + yMax) / 2;
            z0 = (zMin + zMax) / 2;

            topH = int((zMax - z0) / dis_resolution) + 1;      // 任一方格范围 [0, 1)
            rightW = int((yMax - y0) / dis_resolution) + 1;
            bottomH = int((z0 - zMin) / dis_resolution) + 1;
            leftW = int((y0 - yMin) / dis_resolution) + 1;

            int H = topH + bottomH;
            int W = rightW + leftW;

            y_start = y0 - dis_resolution * leftW;
            z_start = z0 + dis_resolution * topH;

            for (int i = 0; i < H; i++)
            {
                std::deque<grid> temp_row;
                for (int j = 0; j < W; j++)
                {
                    temp_row.push_back(build_new_grid(y_start, z_start, i, j));
                }
                uniform_grid.push_back(temp_row);
            }

            // 记录网格顶点在点列表中的索引
            for (int i = 0; i < 2*H+1; i++)
            {
                std::deque<int> temp_row;
                for (int j = 0; j < 2*W+1; j++)
                {
                    temp_row.push_back(-1);
                }
                tri_vertex_idx_grid.push_back(temp_row);
            }

            for (int i = 0; i < projected_pts.rows(); i++)
            {
                double u_d = (z_start - projected_pts(i, 2)) / dis_resolution;
                double v_d = (projected_pts(i, 1) - y_start) / dis_resolution;

                int u = int(u_d); // [)
                int v = int(v_d);

                // 更新x方向均值、点数
                uniform_grid[u][v].x_avg =
                        (uniform_grid[u][v].pts_num * uniform_grid[u][v].x_avg + projected_pts(i, 0))
                        / (uniform_grid[u][v].pts_num + 1);
                uniform_grid[u][v].pts_num++;

                // 更新网格类型 外部->内部/边缘 --- 网格类型更改不完全正确，需后续处理-》MC_mesh()
                if (uniform_grid[u][v].node_type == 0 && uniform_grid[u][v].pts_num > grid_ptsNum_threshold)
                {
                    // 0->1
                    uniform_grid[u][v].node_type = 1;
                    if (!uniform_grid[u][v].if_grid_need_update)
                    {
                        uniform_grid[u][v].if_grid_need_update = true;
                        std::array<int, 2> temp_grid_uv;
                        temp_grid_uv[0] = u;
                        temp_grid_uv[1] = v;
                        grid_need_update_list.push_back(temp_grid_uv);
                    }
                    // 2->1?
                    for (int j = 0; j < 9; j++)
                    {
                        if (j == 4) {continue;}

                        int m = j / 3 - 1;
                        int n = j % 3 - 1;

                        if (u + m < 0 || u + m >= H || v + n < 0 || v + n >= W)
                        {
                            continue;
                        }
                        if (uniform_grid[u+m][v+n].node_type == 2)
                        {
                            if (!uniform_grid[u+m][v+n].if_grid_need_update)
                            {
                                uniform_grid[u+m][v+n].if_grid_need_update = true;
                                std::array<int, 2> temp_grid_uv;
                                temp_grid_uv[0] = u+m;
                                temp_grid_uv[1] = v+n;
                                grid_need_update_list.push_back(temp_grid_uv);
                            }
                        }
                    }
                }

                // 更新网格边界插值点的位置 (0, 1)      TODO 0.1是否合适
                if ((u_d - u) < 0.1)   // yMin_n yMax_n
                {
                    uniform_grid[u][v].edge_vertex(0, 1) = std::min(uniform_grid[u][v].edge_vertex(0, 1), projected_pts(i, 1));
                    uniform_grid[u][v].edge_vertex(0, 0) = std::max(uniform_grid[u][v].edge_vertex(0, 0), projected_pts(i, 1));
                }
                if ((u_d - u) > 0.9)   // yMin_0 yMax_0
                {
                    uniform_grid[u][v].edge_vertex(2, 0) = std::min(uniform_grid[u][v].edge_vertex(2, 0), projected_pts(i, 1));
                    uniform_grid[u][v].edge_vertex(2, 1) = std::max(uniform_grid[u][v].edge_vertex(2, 1), projected_pts(i, 1));
                }
                if ((v_d - v) < 0.1)   // zMin_0 zMax_0
                {
                    uniform_grid[u][v].edge_vertex(3, 1) = std::min(uniform_grid[u][v].edge_vertex(3, 1), projected_pts(i, 2));
                    uniform_grid[u][v].edge_vertex(3, 0) = std::max(uniform_grid[u][v].edge_vertex(3, 0), projected_pts(i, 2));
                }
                if ((v_d - v) > 0.9)   // zMin_n zMax_n
                {
                    uniform_grid[u][v].edge_vertex(1, 0) = std::min(uniform_grid[u][v].edge_vertex(1, 0), projected_pts(i, 2));
                    uniform_grid[u][v].edge_vertex(1, 1) = std::max(uniform_grid[u][v].edge_vertex(1, 1), projected_pts(i, 2));
                }
            }

//            MC_mesh();
            MC_mesh_fast();

//            vertex_and_face_list_extract();
//            vertex_and_face_list_update();
            vertex_and_face_list_update_fast();

            if_first_decimate = false;
            last_grid_decimate_num = ptcl_all->size();
        }
        // 平面更新后再次处理
        else
        {
            // 新添加的平面点投影到二维平面，按第一次投影定义的坐标变换
            int new_pt_num = ptcl_all->size() - last_grid_decimate_num;
            Eigen::MatrixXd eigen_cloud(new_pt_num, 3);
            for (int i = last_grid_decimate_num; i < ptcl_all->size(); ++i)
            {
                eigen_cloud(i-last_grid_decimate_num, 0) = ptcl_all->points[i].x;
                eigen_cloud(i-last_grid_decimate_num, 1) = ptcl_all->points[i].y;
                eigen_cloud(i-last_grid_decimate_num, 2) = ptcl_all->points[i].z;
            }
            Eigen::MatrixXd projected_pts(new_pt_num, 3);
            projected_pts = (eigen_cloud.rowwise() - quadtree_center.transpose()) * quadtree_axes;
            double new_yMax = projected_pts.col(1).maxCoeff();
            double new_yMin = projected_pts.col(1).minCoeff();
            double new_zMax = projected_pts.col(2).maxCoeff();
            double new_zMin = projected_pts.col(2).minCoeff();

            yMax = std::max(yMax, new_yMax);
            yMin = std::min(yMin, new_yMin);
            zMax = std::max(zMax, new_zMax);
            zMin = std::min(zMin, new_zMin);

            int new_topH = int((zMax - z0) / dis_resolution) + 1;
            int new_rightW = int((yMax - y0) / dis_resolution) + 1;
            int new_bottomH = int((z0 - zMin) / dis_resolution) + 1;
            int new_leftW = int((y0 - yMin) / dis_resolution) + 1;

            int H = new_topH + new_bottomH;
            int W = new_rightW + new_leftW;

            int u_start = new_topH - topH;
            int u_end = H - (new_bottomH - bottomH) - 1;
            int v_start = new_leftW - leftW;
            int v_end = W - (new_rightW - rightW) - 1;

            y_start = y0 - dis_resolution * new_leftW;
            z_start = z0 + dis_resolution * new_topH;

            for (int i = 0; i < H; i++)
            {
                if (i < u_start)
                {
                    std::deque<grid> temp_row;
                    for (int j = 0; j < W; j++)
                    {
                        temp_row.push_back(build_new_grid(y_start, z_start, u_start - i - 1, j));
                    }
                    uniform_grid.push_front(temp_row);
                }
                else if (i > u_end)
                {
                    std::deque<grid> temp_row;
                    for (int j = 0; j < W; j++)
                    {
                        temp_row.push_back(build_new_grid(y_start, z_start, i, j));
                    }
                    uniform_grid.push_back(temp_row);
                }
                else
                {
                    for (int j = 0; j < W; j++)
                    {
                        if (j < v_start)
                        {
                            uniform_grid[i].push_front(build_new_grid(y_start, z_start, i, v_start - j - 1));
                        }
                        else if (j > v_end)
                        {
                            uniform_grid[i].push_back(build_new_grid(y_start, z_start, i, j));
                        }
                        else
                        {
                            // 内部不处理
                        }
                    }
                }
            }

            // 记录网格顶点在点列表中的索引 -- 更新大小
            int u_start_2 = 2 * (new_topH - topH);
            int u_end_2 = 2 * H - 2 * (new_bottomH - bottomH);
            int v_start_2 = 2 * (new_leftW - leftW);
            int v_end_2 = 2 * W - 2 * (new_rightW - rightW);
            for (int i = 0; i < 2*H+1; i++)
            {
                if (i < u_start_2)
                {
                    std::deque<int> temp_row;
                    for (int j = 0; j < 2*W+1; j++)
                    {
                        temp_row.push_back(-1);
                    }
                    tri_vertex_idx_grid.push_front(temp_row);
                }
                else if (i > u_end_2)
                {
                    std::deque<int> temp_row;
                    for (int j = 0; j < 2*W+1; j++)
                    {
                        temp_row.push_back(-1);
                    }
                    tri_vertex_idx_grid.push_back(temp_row);
                }
                else
                {
                    for (int j = 0; j < 2*W+1; j++)
                    {
                        if (j < v_start_2)
                        {
                            tri_vertex_idx_grid[i].push_front(-1);
                        }
                        else if (j > v_end_2)
                        {
                            tri_vertex_idx_grid[i].push_back(-1);
                        }
                        else
                        {
                            // 内部不处理
                        }
                    }
                }
            }

            for (int i = 0; i < projected_pts.rows(); i++)
            {
                double u_d = (z_start - projected_pts(i, 2)) / dis_resolution;
                double v_d = (projected_pts(i, 1) - y_start) / dis_resolution;

                int u = int(u_d); // [)
                int v = int(v_d);

                // 更新x方向均值、点数
                uniform_grid[u][v].x_avg =
                        (uniform_grid[u][v].pts_num * uniform_grid[u][v].x_avg + projected_pts(i, 0))
                        / (uniform_grid[u][v].pts_num + 1);
                uniform_grid[u][v].pts_num++;

                // 更新网格类型 外部->内部/边缘
                if (uniform_grid[u][v].node_type == 0 && uniform_grid[u][v].pts_num > grid_ptsNum_threshold)
                {
                    // 0->1
                    uniform_grid[u][v].node_type = 1;
                    if (!uniform_grid[u][v].if_grid_need_update)
                    {
                        uniform_grid[u][v].if_grid_need_update = true;
                        std::array<int, 2> temp_grid_uv;
                        temp_grid_uv[0] = u;
                        temp_grid_uv[1] = v;
                        grid_need_update_list.push_back(temp_grid_uv);
                    }
                    // 2->1?
                    for (int j = 0; j < 9; j++)
                    {
                        if (j == 4) {continue;}

                        int m = j / 3 - 1;
                        int n = j % 3 - 1;

                        if (u + m < 0 || u + m >= H || v + n < 0 || v + n >= W)
                        {
                            continue;
                        }
                        if (uniform_grid[u+m][v+n].node_type == 2)
                        {
                            if (!uniform_grid[u+m][v+n].if_grid_need_update)
                            {
                                uniform_grid[u+m][v+n].if_grid_need_update = true;
                                std::array<int, 2> temp_grid_uv;
                                temp_grid_uv[0] = u+m;
                                temp_grid_uv[1] = v+n;
                                grid_need_update_list.push_back(temp_grid_uv);
                            }
                        }
                    }
                }

                // 更新网格边界插值点的位置 (0, 1)
                if ((u_d - u) < 0.1)   // yMin_n
                {
                    uniform_grid[u][v].edge_vertex(0, 1) = std::min(uniform_grid[u][v].edge_vertex(0, 1), projected_pts(i, 1));
                    uniform_grid[u][v].edge_vertex(0, 0) = std::max(uniform_grid[u][v].edge_vertex(0, 0), projected_pts(i, 1));
                }
                if ((u_d - u) > 0.9)   // yMin_0
                {
                    uniform_grid[u][v].edge_vertex(2, 0) = std::min(uniform_grid[u][v].edge_vertex(2, 0), projected_pts(i, 1));
                    uniform_grid[u][v].edge_vertex(2, 1) = std::max(uniform_grid[u][v].edge_vertex(2, 1), projected_pts(i, 1));
                }
                if ((v_d - v) < 0.1)   // zMin_0
                {
                    uniform_grid[u][v].edge_vertex(3, 1) = std::min(uniform_grid[u][v].edge_vertex(3, 1), projected_pts(i, 2));
                    uniform_grid[u][v].edge_vertex(3, 0) = std::max(uniform_grid[u][v].edge_vertex(3, 0), projected_pts(i, 2));
                }
                if ((v_d - v) > 0.9)   // zMin_n
                {
                    uniform_grid[u][v].edge_vertex(1, 0) = std::min(uniform_grid[u][v].edge_vertex(1, 0), projected_pts(i, 2));
                    uniform_grid[u][v].edge_vertex(1, 1) = std::max(uniform_grid[u][v].edge_vertex(1, 1), projected_pts(i, 2));
                }
            }

            topH = new_topH;
            bottomH = new_bottomH;
            leftW = new_leftW;
            rightW = new_rightW;

//            MC_mesh();
            MC_mesh_fast();

//            vertex_and_face_list_extract();
//            vertex_and_face_list_update();
            vertex_and_face_list_update_fast();
            last_grid_decimate_num = ptcl_all->size();
        }
    }
}

// 借鉴Marching Cube算法思想提取mesh
// 首先更新网格的状态，然后根据周围8个网格的状态确定：三角形顶点列表
void MeshFragment::MC_mesh_fast()
{
    int H = topH + bottomH;
    int W = rightW + leftW;
    faces_to_delete_frame.clear();
    for (auto it_uv : grid_need_update_list)
    {
        int u = it_uv[0];
        int v = it_uv[1];
        uniform_grid[u][v].if_grid_need_update = false;

        // 如果是边缘网格则先删除已有的三角形
        if (uniform_grid[u][v].node_type == 2)
        {
            uniform_grid[u][v].if_mesh_need_delete = true;    // TODO
            for (auto it : uniform_grid[u][v].tri_idx_list)
            {
                faces_to_delete.push_back(it);
                faces_to_delete_frame.push_back(it);
            }
            uniform_grid[u][v].tri_idx_list.clear();
        }

        // 添加三角形
        uniform_grid[u][v].if_mesh_need_add = true;
        // 判断内部还是边缘
        bool if_inside = true;
        Eigen::Vector4i grid_vertex_state = {1, 1, 1, 1};   // 左上角开始 顺时针 TODO 改为根据周围四个node类型判断
        for (int j = 0; j < 9; j++)
        {
            if (j == 4) {continue;}

            int m = j / 3 - 1;
            int n = j % 3 - 1;

            if (u + m < 0 || u + m >= H || v + n < 0 || v + n >= W)
            {
                if_inside = false;
                uniform_grid[u][v].node_type = 2;
                switch (j) {
                    case 0:
                        grid_vertex_state(0) = 0;
                        break;
                    case 1:
                        grid_vertex_state(0) = 0;
                        grid_vertex_state(1) = 0;
                        break;
                    case 2:
                        grid_vertex_state(1) = 0;
                        break;
                    case 3:
                        grid_vertex_state(0) = 0;
                        grid_vertex_state(3) = 0;
                        break;
                    case 5:
                        grid_vertex_state(1) = 0;
                        grid_vertex_state(2) = 0;
                        break;
                    case 6:
                        grid_vertex_state(3) = 0;
                        break;
                    case 7:
                        grid_vertex_state(2) = 0;
                        grid_vertex_state(3) = 0;
                        break;
                    case 8:
                        grid_vertex_state(2) = 0;
                        break;
                    default:
                        break;
                }
                continue;
            }

            if (uniform_grid[u+m][v+n].node_type == 0)
            {
                if_inside = false;
                uniform_grid[u][v].node_type = 2;
                switch (j) {
                    case 0:
                        grid_vertex_state(0) = 0;
                        break;
                    case 1:
                        grid_vertex_state(0) = 0;
                        grid_vertex_state(1) = 0;
                        break;
                    case 2:
                        grid_vertex_state(1) = 0;
                        break;
                    case 3:
                        grid_vertex_state(0) = 0;
                        grid_vertex_state(3) = 0;
                        break;
                    case 5:
                        grid_vertex_state(1) = 0;
                        grid_vertex_state(2) = 0;
                        break;
                    case 6:
                        grid_vertex_state(3) = 0;
                        break;
                    case 7:
                        grid_vertex_state(2) = 0;
                        grid_vertex_state(3) = 0;
                        break;
                    case 8:
                        grid_vertex_state(2) = 0;
                        break;
                    default:
                        break;
                }
            }
        }

        // 根据内部、边缘设计不同的mesh策略
        if (if_inside)
        {
            uniform_grid[u][v].node_type = 1;
            uniform_grid[u][v].vertex_list.clear();
            for (int j = 0; j < 4; j++)
            {
                Eigen::Vector2d temp_vertex;
                temp_vertex << uniform_grid[u][v].grid_vertex_list(j, 0),
                        uniform_grid[u][v].grid_vertex_list(j, 1);
                uniform_grid[u][v].vertex_list.push_back(temp_vertex);
            }
        }
        else
        {
            uniform_grid[u][v].node_type = 2;
            uniform_grid[u][v].vertex_list.clear();
            for (int j = 0; j < 4; j++)
            {
                if (grid_vertex_state(j) == 1)
                {
                    // 1-
                    Eigen::Vector2d temp_vertex;
                    temp_vertex(0) = uniform_grid[u][v].grid_vertex_list(j, 0);
                    temp_vertex(1) = uniform_grid[u][v].grid_vertex_list(j, 1);
                    uniform_grid[u][v].vertex_list.push_back(temp_vertex);
                    // 1-0
                    if (grid_vertex_state( (j+1)%4 ) == 0)
                    {
                        Eigen::Vector2d temp_vertex_1;
                        if (j == 0 || j == 2)
                        {
                            // 若还没有相应的数据则取中点
                            if (uniform_grid[u][v].edge_vertex(j, 0) == uniform_grid[u][v].grid_vertex_list(j, 0))
                            {
                                temp_vertex_1(0) = uniform_grid[u][v].grid_center(0);
                            }
                            else
                            {
                                temp_vertex_1(0) = uniform_grid[u][v].edge_vertex(j, 0);
                            }

                            temp_vertex_1(1) = uniform_grid[u][v].grid_vertex_list(j, 1);
                        }
                        else
                        {
                            temp_vertex_1(0) = uniform_grid[u][v].grid_vertex_list(j, 0);
                            if (uniform_grid[u][v].edge_vertex(j, 0) == uniform_grid[u][v].grid_vertex_list(j, 1))
                            {
                                temp_vertex_1(1) = uniform_grid[u][v].grid_center(1);
                            }
                            else
                            {
                                temp_vertex_1(1) = uniform_grid[u][v].edge_vertex(j, 0);
                            }
                        }
                        uniform_grid[u][v].vertex_list.push_back(temp_vertex_1);
                    }
                }
                else
                {
                    // 0-1
                    if (grid_vertex_state( (j+1)%4 ) == 1)
                    {
                        Eigen::Vector2d temp_vertex;

                        if (j == 0 || j == 2)
                        {
                            if (uniform_grid[u][v].edge_vertex(j, 1) == uniform_grid[u][v].grid_vertex_list((j+1)%4, 0))
                            {
                                temp_vertex(0) = uniform_grid[u][v].grid_center(0);
                            }
                            else
                            {
                                temp_vertex(0) = uniform_grid[u][v].edge_vertex(j, 1);
                            }
                            temp_vertex(1) = uniform_grid[u][v].grid_vertex_list(j, 1);
                        }
                        else
                        {
                            temp_vertex(0) = uniform_grid[u][v].grid_vertex_list(j, 0);
                            if (uniform_grid[u][v].edge_vertex(j, 1) == uniform_grid[u][v].grid_vertex_list((j+1)%4, 1))
                            {
                                temp_vertex(1) = uniform_grid[u][v].grid_center(1);
                            }
                            else
                            {
                                temp_vertex(1) = uniform_grid[u][v].edge_vertex(j, 1);
                            }
                        }
                        uniform_grid[u][v].vertex_list.push_back(temp_vertex);
                    }
                }
            }
        }
    }
}

// 每次增量提取点列表、面列表
// 首先提取点列表，要求不同网格公用顶点坐标一致，然后根据每个网格的三角形顶点列表中点的数量(3 4 5 6)构建三角剖分
void MeshFragment::vertex_and_face_list_update_fast()
{

    int H = topH + bottomH;
    int W = rightW + leftW;
    double pt_intensity = 100.0;
    pts_inner_add.clear();
    faces_with_edge_vertex.clear();

    for (auto it_uv : grid_need_update_list)
    {
        int u = it_uv[0];
        int v = it_uv[1];

        // 添加/更新点列表，并确定点的索引
        std::vector<int> temp_ind_list;
        for (int i = 0; i < uniform_grid[u][v].vertex_list.size(); i++)
        {
            // 确定该顶点在 tri_vertex_idx_grid 中的索引
            int idx_u, idx_v;
            bool if_grid_vertex = true;
            if (uniform_grid[u][v].vertex_list[i](1) == uniform_grid[u][v].grid_vertex_list(0,1))
            {
                idx_u = u * 2;
            }
            else if (uniform_grid[u][v].vertex_list[i](1) == uniform_grid[u][v].grid_vertex_list(2,1))
            {
                idx_u = (u + 1) * 2;
            }
            else
            {
                idx_u = u * 2 + 1;
                if_grid_vertex = false;
            }

            if (uniform_grid[u][v].vertex_list[i](0) == uniform_grid[u][v].grid_vertex_list(0,0))
            {
                idx_v = v * 2;
            }
            else if (uniform_grid[u][v].vertex_list[i](0) == uniform_grid[u][v].grid_vertex_list(1,0))
            {
                idx_v = (v + 1) * 2;
            }
            else
            {
                idx_v = v * 2 + 1;
                if_grid_vertex = false;
            }

            // 添加新点 / 更新点的坐标
            Eigen::Vector3d pt_local, pt_global;
            if (if_grid_vertex)
            {
                // 考虑边缘
                if (idx_u == 0)
                {
                    if (idx_v == 0)
                    {
                        pt_local(0) = uniform_grid[0][0].x_avg;
                    }
                    else if (idx_v == 2*W)
                    {
                        pt_local(0) = uniform_grid[0][W-1].x_avg;
                    }
                    else
                    {
                        pt_local(0) = (uniform_grid[0][idx_v/2].x_avg
                                       + uniform_grid[0][idx_v/2-1].x_avg) / 2;
                    }
                }
                else if (idx_u == 2*H)
                {
                    if (idx_v == 0)
                    {
                        pt_local(0) = uniform_grid[H-1][0].x_avg;
                    }
                    else if (idx_v == 2*W)
                    {
                        pt_local(0) = uniform_grid[H-1][W-1].x_avg;
                    }
                    else
                    {
                        pt_local(0) = (uniform_grid[H-1][idx_v/2].x_avg
                                       + uniform_grid[H-1][idx_v/2-1].x_avg) / 2;
                    }
                }
                else
                {
                    if (idx_v == 0)
                    {
                        pt_local(0) = (uniform_grid[idx_v/2][0].x_avg
                                       + uniform_grid[idx_v/2-1][0].x_avg) / 2;
                    }
                    else if (idx_v == 2*W)
                    {
                        pt_local(0) = (uniform_grid[idx_v/2][W-1].x_avg
                                       + uniform_grid[idx_v/2-1][W-1].x_avg) / 2;
                    }
                    else
                    {
                        pt_local(0) = (uniform_grid[idx_u/2][idx_v/2].x_avg
                                       + uniform_grid[idx_u/2-1][idx_v/2].x_avg
                                       + uniform_grid[idx_u/2][idx_v/2-1].x_avg
                                       + uniform_grid[idx_u/2-1][idx_v/2-1].x_avg) / 4;
                    }
                }

                // 坐标转到世界坐标系下
                pt_local(1) = uniform_grid[u][v].vertex_list[i](0);
                pt_local(2) = uniform_grid[u][v].vertex_list[i](1);
                pt_global = quadtree_axes * pt_local + quadtree_center;

                // 判断是不是新点
                if (tri_vertex_idx_grid[idx_u][idx_v] == -1)
                {
                    // 若是则添加新点，
                    pcl::PointXYZINormal temp_pt;
                    temp_pt.x = pt_global(0);
                    temp_pt.y = pt_global(1);
                    temp_pt.z = pt_global(2);
                    temp_pt.normal_x = normal_vector(0);
                    temp_pt.normal_y = normal_vector(1);
                    temp_pt.normal_z = normal_vector(2);
                    temp_pt.intensity = pt_intensity;

                    mtx_ptcl_grid.lock();
                    ptcl_grid->push_back(temp_pt);
                    mtx_ptcl_grid.unlock();
                    pts_state.push_back(1);
                    tri_vertex_idx_grid[idx_u][idx_v] = ptcl_grid->size() - 1;
                    temp_ind_list.push_back(tri_vertex_idx_grid[idx_u][idx_v]);

                    pts_inner_add.push_back(tri_vertex_idx_grid[idx_u][idx_v]);
                }
                else
                {
                    // 若不是则更新点的坐标
                    mtx_ptcl_grid.lock();
                    int temp_pt_idx = tri_vertex_idx_grid[idx_u][idx_v];
                    ptcl_grid->points[temp_pt_idx].x = pt_global(0);
                    mtx_ptcl_grid.unlock();
                    temp_ind_list.push_back(temp_pt_idx);
                }
            }
            // 不是网格顶点
            else
            {
                // TODO 取平均 会有之前点的影响 不是很严谨
                // 修改：根据idx_u、idx_v确定哪条边顺时针0-3，根据uniform_grid[u][v].grid_vertex_state确定min、max
                pt_local(0) = uniform_grid[u][v].x_avg;
                pt_local(1) = uniform_grid[u][v].vertex_list[i](0);
                pt_local(2) = uniform_grid[u][v].vertex_list[i](1);
                pt_global = quadtree_axes * pt_local + quadtree_center;

                // 判断是不是新点
                if (tri_vertex_idx_grid[idx_u][idx_v] == -1)
                {
                    pcl::PointXYZINormal temp_pt;
                    temp_pt.x = pt_global(0);
                    temp_pt.y = pt_global(1);
                    temp_pt.z = pt_global(2);
                    temp_pt.normal_x = normal_vector(0);
                    temp_pt.normal_y = normal_vector(1);
                    temp_pt.normal_z = normal_vector(2);
                    temp_pt.intensity = pt_intensity;

                    mtx_ptcl_grid.lock();
                    ptcl_grid->push_back(temp_pt);
                    mtx_ptcl_grid.unlock();
                    pts_state.push_back(2);
                    tri_vertex_idx_grid[idx_u][idx_v] = ptcl_grid->size() - 1;
                    temp_ind_list.push_back(tri_vertex_idx_grid[idx_u][idx_v]);

                    pts_inner_add.push_back(tri_vertex_idx_grid[idx_u][idx_v]);
                }
                else
                {
                    // 若不是则更新点的坐标
                    int temp_pt_idx = tri_vertex_idx_grid[idx_u][idx_v];

                    mtx_ptcl_grid.lock();
                    ptcl_grid->points[temp_pt_idx].x += pt_global(0);
                    ptcl_grid->points[temp_pt_idx].y += pt_global(1);
                    ptcl_grid->points[temp_pt_idx].z += pt_global(2);
                    ptcl_grid->points[temp_pt_idx].x /= 2;
                    ptcl_grid->points[temp_pt_idx].y /= 2;
                    ptcl_grid->points[temp_pt_idx].z /= 2;
                    mtx_ptcl_grid.unlock();

                    temp_ind_list.push_back(temp_pt_idx);
                }
            }
        }

        // 添加三角面片，逆时针
        if (uniform_grid[u][v].vertex_list.size() == 3)
        {
            std::array<int, 3> temp_tri_idx;
            temp_tri_idx[0] = temp_ind_list[0];
            temp_tri_idx[1] = temp_ind_list[2];
            temp_tri_idx[2] = temp_ind_list[1];
            faces_list.push_back(temp_tri_idx);
            uniform_grid[u][v].tri_idx_list.push_back(faces_list.size() - 1);

            if (uniform_grid[u][v].node_type == 2)
            {
                faces_with_edge_vertex.push_back(faces_list.size() - 1);
            }
        }
        else if (uniform_grid[u][v].vertex_list.size() == 4)
        {
            std::array<int, 3> temp_tri_idx_1, temp_tri_idx_2;
            temp_tri_idx_1[0] = temp_ind_list[0];
            temp_tri_idx_1[1] = temp_ind_list[2];
            temp_tri_idx_1[2] = temp_ind_list[1];
            faces_list.push_back(temp_tri_idx_1);

            temp_tri_idx_2[0] = temp_ind_list[2];
            temp_tri_idx_2[1] = temp_ind_list[0];
            temp_tri_idx_2[2] = temp_ind_list[3];
            faces_list.push_back(temp_tri_idx_2);
            uniform_grid[u][v].tri_idx_list.push_back(faces_list.size() - 2);
            uniform_grid[u][v].tri_idx_list.push_back(faces_list.size() - 1);

            if (uniform_grid[u][v].node_type == 1)
            {
                uniform_grid[u][v].ver_idx_list.clear();
                uniform_grid[u][v].ver_idx_list.push_back(temp_ind_list[0]);
                uniform_grid[u][v].ver_idx_list.push_back(temp_ind_list[1]);
                uniform_grid[u][v].ver_idx_list.push_back(temp_ind_list[2]);
                uniform_grid[u][v].ver_idx_list.push_back(temp_ind_list[3]);
            }

            if (uniform_grid[u][v].node_type == 2)
            {
                faces_with_edge_vertex.push_back(faces_list.size() - 2);
                faces_with_edge_vertex.push_back(faces_list.size() - 1);
            }
        }
        else if (uniform_grid[u][v].vertex_list.size() == 5)
        {
            std::array<int, 3> temp_tri_idx_1, temp_tri_idx_2, temp_tri_idx_3;
            temp_tri_idx_1[0] = temp_ind_list[0];
            temp_tri_idx_1[1] = temp_ind_list[2];
            temp_tri_idx_1[2] = temp_ind_list[1];
            faces_list.push_back(temp_tri_idx_1);

            temp_tri_idx_2[0] = temp_ind_list[2];
            temp_tri_idx_2[1] = temp_ind_list[4];
            temp_tri_idx_2[2] = temp_ind_list[3];
            faces_list.push_back(temp_tri_idx_2);

            temp_tri_idx_3[0] = temp_ind_list[0];
            temp_tri_idx_3[1] = temp_ind_list[4];
            temp_tri_idx_3[2] = temp_ind_list[2];
            faces_list.push_back(temp_tri_idx_3);
            uniform_grid[u][v].tri_idx_list.push_back(faces_list.size() - 3);
            uniform_grid[u][v].tri_idx_list.push_back(faces_list.size() - 2);
            uniform_grid[u][v].tri_idx_list.push_back(faces_list.size() - 1);

            if (uniform_grid[u][v].node_type == 2)
            {
                faces_with_edge_vertex.push_back(faces_list.size() - 3);
                faces_with_edge_vertex.push_back(faces_list.size() - 2);
                faces_with_edge_vertex.push_back(faces_list.size() - 1);
            }
        }
        else if (uniform_grid[u][v].vertex_list.size() == 6)
        {
            std::array<int, 3> temp_tri_idx_1, temp_tri_idx_2, temp_tri_idx_3, temp_tri_idx_4;
            temp_tri_idx_1[0] = temp_ind_list[0];
            temp_tri_idx_1[1] = temp_ind_list[2];
            temp_tri_idx_1[2] = temp_ind_list[1];
            faces_list.push_back(temp_tri_idx_1);

            temp_tri_idx_2[0] = temp_ind_list[2];
            temp_tri_idx_2[1] = temp_ind_list[4];
            temp_tri_idx_2[2] = temp_ind_list[3];
            faces_list.push_back(temp_tri_idx_2);

            temp_tri_idx_3[0] = temp_ind_list[4];
            temp_tri_idx_3[1] = temp_ind_list[0];
            temp_tri_idx_3[2] = temp_ind_list[5];
            faces_list.push_back(temp_tri_idx_3);

            temp_tri_idx_4[0] = temp_ind_list[0];
            temp_tri_idx_4[1] = temp_ind_list[4];   // TODO 修改 3-》2
            temp_tri_idx_4[2] = temp_ind_list[2];
            faces_list.push_back(temp_tri_idx_4);
            uniform_grid[u][v].tri_idx_list.push_back(faces_list.size() - 4);
            uniform_grid[u][v].tri_idx_list.push_back(faces_list.size() - 3);
            uniform_grid[u][v].tri_idx_list.push_back(faces_list.size() - 2);
            uniform_grid[u][v].tri_idx_list.push_back(faces_list.size() - 1);

            if (uniform_grid[u][v].node_type == 2)
            {
                faces_with_edge_vertex.push_back(faces_list.size() - 4);
                faces_with_edge_vertex.push_back(faces_list.size() - 3);
                faces_with_edge_vertex.push_back(faces_list.size() - 2);
                faces_with_edge_vertex.push_back(faces_list.size() - 1);
            }
        }
        else
        {
            // 不处理
        }
    }

    // 提取传递给非平面部分的数据
    std::set<int> pts_edge_delete_wo_rep, pts_edge_add_wo_rep;
    pts_edge_delete.clear();
    pts_edge_add.clear();
    pts_edge_update.clear();
    faces_with_2_edge_vertex.clear();
    tri_delete_string_idx.clear();

    for (auto it : faces_to_delete_frame)
    {
        pts_edge_delete_wo_rep.insert(faces_list[it][0]);
        pts_edge_delete_wo_rep.insert(faces_list[it][1]);
        pts_edge_delete_wo_rep.insert(faces_list[it][2]);

        // 为要删除的三角形 根据顶点索引 构建索引
        std::vector<int> temp_tri_v_idx;
        for (int l = 0; l < 3; l++)
        {
            if (plane_vidx_to_noplnae_vidx.find(faces_list[it][l]) == plane_vidx_to_noplnae_vidx.end())
            {
                break;
            }
            temp_tri_v_idx.push_back(plane_vidx_to_noplnae_vidx[faces_list[it][l]]);
        }
        if (temp_tri_v_idx.size() != 3) { continue; }
        // 从小到大排序
        std::sort(temp_tri_v_idx.begin(), temp_tri_v_idx.end());
        std::string tri_key = std::to_string(temp_tri_v_idx[0]) + ' '
                              + std::to_string(temp_tri_v_idx[1]) + ' ' + std::to_string(temp_tri_v_idx[2]);
        tri_delete_string_idx[tri_key] = 1;
    }

    for (auto it : faces_with_edge_vertex)
    {
        faces_with_2_edge_vertex.push_back(it);
        pts_edge_add_wo_rep.insert(faces_list[it][0]);
        pts_edge_add_wo_rep.insert(faces_list[it][1]);
        pts_edge_add_wo_rep.insert(faces_list[it][2]);
    }

    // 寻找相同点-》需更新位置的点
    for (auto it : pts_edge_delete_wo_rep)
    {
        if (pts_state[it] == 1)
        {
            continue;
        }

        if (pts_edge_add_wo_rep.find(it) != pts_edge_add_wo_rep.end())
        {
            pts_edge_update.push_back(it);
        }
        else
        {
            pts_edge_delete.push_back(it);
        }
    }

    for (auto it : pts_edge_add_wo_rep)
    {
        auto temp_it = std::find(pts_edge_update.begin(), pts_edge_update.end(), it);
        if (temp_it == pts_edge_update.end())
        {
            pts_edge_add.push_back(it);
        }
    }

}

// 四叉树简化
void MeshFragment::quadtree_decimate()
{
    if (ptcl_grid->points.size() > 0)
    {
        QuadTree* QT = new QuadTree();
        QT->uniform_grid = &uniform_grid;
        QT->build();
        QT->mergeQuadTree(QT->Root);
        QT->faces_extract(QT->Root);

        for (auto it : QT->faces_idx_to_delete)
        {
            faces_to_delete.push_back(it);
        }

        for (auto it : QT->faces_to_add)
        {
            faces_list.push_back(it);
        }

        std::vector<std::vector<int>> grid_vertex_count;
        for (int i = 0; i < uniform_grid.size()+1; i++)
        {
            std::vector<int> temp_row;
            for (int j = 0; j < uniform_grid[0].size()+1; j++)
            {
                temp_row.push_back(0);
            }
            grid_vertex_count.push_back(temp_row);
        }
        QT->grid_vertex_count = &grid_vertex_count;
        QT->tri_vertex_idx_grid = &tri_vertex_idx_grid;
        QT->grid_vertex_extract(QT->Root);

        QT->clearTree(QT->Root);

        // 调整位于边上的顶点的坐标：两端顶点坐标的线性插值
        // 整体思路（1、确定简化后的网格顶点；2、遍历简化后的网格顶点，根据出现次数确定需要调整的顶点；3、找相应的大三角形的顶点插值）
        // 另外需要注意的是，有些顶点缺少相应的固定顶点，因此需要调整两次
        std::vector<int> another_adjust; // 需要再调整一轮的点信息
        int temp_H = uniform_grid.size();
        int temp_W = uniform_grid[0].size();
        for (int i = 0; i < temp_H; i++)
        {
            for (int j = 0; j < temp_W; j++)
            {
                // 不需要调整
                if (grid_vertex_count[i][j] != 2) { continue; }

                // 判断搜索方向，分为四种情况
                int max_level = uniform_grid[i][j].cell_level;
                if (max_level < uniform_grid[i-1][j].cell_level) { max_level = uniform_grid[i-1][j].cell_level; }
                if (max_level < uniform_grid[i-1][j-1].cell_level) { max_level = uniform_grid[i-1][j-1].cell_level; }
                if (max_level < uniform_grid[i][j-1].cell_level) { max_level = uniform_grid[i][j-1].cell_level; }

                if (uniform_grid[i][j].cell_level == max_level && uniform_grid[i][j-1].cell_level == max_level)
                {
                    bool if_need_another = false;
                    int l_1 = 1, l_2 = 1;
                    while (j - l_1 >= 0)
                    {
                        if (uniform_grid[i][j-l_1-1].cell_level != max_level) {
                            if_need_another = true;  // 需要二次调整
                            break;
                        }
                        if (grid_vertex_count[i][j-l_1] == 4) { break; }
                        l_1++;
                    }
                    while (j + l_2 < temp_W)
                    {
                        if (uniform_grid[i][j+l_2].cell_level != max_level) {
                            if_need_another = true;
                            break;
                        }
                        if (grid_vertex_count[i][j+l_2] == 4) { break; }
                        l_2++;
                    }
                    if ((l_1 + l_2) == std::pow(2, max_level-1))
                    {
                        int temp_v0_idx = tri_vertex_idx_grid[i*2][j*2];
                        int temp_v1_idx = tri_vertex_idx_grid[i*2][j*2-2*l_1];
                        int temp_v2_idx = tri_vertex_idx_grid[i*2][j*2+2*l_2];

                        if (if_need_another) // 保存需要二次调整的顶点的信息
                        {
                            another_adjust.push_back(temp_v0_idx);
                            another_adjust.push_back(temp_v1_idx);
                            another_adjust.push_back(temp_v2_idx);
                            another_adjust.push_back(l_1);
                            another_adjust.push_back(l_2);
                        }

                        ptcl_grid->points[temp_v0_idx].x = (ptcl_grid->points[temp_v1_idx].x * l_2 + ptcl_grid->points[temp_v2_idx].x * l_1) / (l_1 + l_2);
                        ptcl_grid->points[temp_v0_idx].y = (ptcl_grid->points[temp_v1_idx].y * l_2 + ptcl_grid->points[temp_v2_idx].y * l_1) / (l_1 + l_2);
                        ptcl_grid->points[temp_v0_idx].z = (ptcl_grid->points[temp_v1_idx].z * l_2 + ptcl_grid->points[temp_v2_idx].z * l_1) / (l_1 + l_2);
                    }
                }
                if (uniform_grid[i][j].cell_level == max_level && uniform_grid[i-1][j].cell_level == max_level)
                {
                    bool if_need_another = false;
                    int l_1 = 1, l_2 = 1;
                    while (i - l_1 >= 0)
                    {
                        if (uniform_grid[i-l_1-1][j].cell_level != max_level) {
                            if_need_another = true;
                            break;
                        }
                        if (grid_vertex_count[i-l_1][j] == 4) { break; }
                        l_1++;
                    }
                    while (i + l_2 < temp_H)
                    {
                        if (uniform_grid[i+l_2][j].cell_level != max_level) {
                            if_need_another = true;
                            break;
                        }
                        if (grid_vertex_count[i+l_2][j] == 4) { break; }
                        l_2++;
                    }

                    if ((l_1 + l_2) == std::pow(2, max_level-1))
                    {
                        int temp_v0_idx = tri_vertex_idx_grid[i*2][j*2];
                        int temp_v1_idx = tri_vertex_idx_grid[i*2-2*l_1][j*2];
                        int temp_v2_idx = tri_vertex_idx_grid[i*2+2*l_2][j*2];

                        if (if_need_another)
                        {
                            another_adjust.push_back(temp_v0_idx);
                            another_adjust.push_back(temp_v1_idx);
                            another_adjust.push_back(temp_v2_idx);
                            another_adjust.push_back(l_1);
                            another_adjust.push_back(l_2);
                        }

                        ptcl_grid->points[temp_v0_idx].x = (ptcl_grid->points[temp_v1_idx].x * l_2 + ptcl_grid->points[temp_v2_idx].x * l_1) / (l_1 + l_2);
                        ptcl_grid->points[temp_v0_idx].y = (ptcl_grid->points[temp_v1_idx].y * l_2 + ptcl_grid->points[temp_v2_idx].y * l_1) / (l_1 + l_2);
                        ptcl_grid->points[temp_v0_idx].z = (ptcl_grid->points[temp_v1_idx].z * l_2 + ptcl_grid->points[temp_v2_idx].z * l_1) / (l_1 + l_2);
                    }
                }
                if (uniform_grid[i-1][j-1].cell_level == max_level && uniform_grid[i][j-1].cell_level == max_level)
                {
                    bool if_need_another = false;
                    int l_1 = 1, l_2 = 1;
                    while (i - l_1 >= 0)
                    {
                        if (uniform_grid[i-l_1-1][j-1].cell_level != max_level) {
                            if_need_another = true;
                            break;
                        }
                        if (grid_vertex_count[i-l_1][j] == 4) { break; }
                        l_1++;
                    }
                    while (i + l_2 < temp_H)
                    {
                        if (uniform_grid[i+l_2][j-1].cell_level != max_level) {
                            if_need_another = true;
                            break;
                        }
                        if (grid_vertex_count[i+l_2][j] == 4) { break; }
                        l_2++;
                    }

                    if ((l_1 + l_2) == std::pow(2, max_level-1))
                    {
                        int temp_v0_idx = tri_vertex_idx_grid[i*2][j*2];
                        int temp_v1_idx = tri_vertex_idx_grid[i*2-2*l_1][j*2];
                        int temp_v2_idx = tri_vertex_idx_grid[i*2+2*l_2][j*2];

                        if (if_need_another)
                        {
                            another_adjust.push_back(temp_v0_idx);
                            another_adjust.push_back(temp_v1_idx);
                            another_adjust.push_back(temp_v2_idx);
                            another_adjust.push_back(l_1);
                            another_adjust.push_back(l_2);
                        }

                        ptcl_grid->points[temp_v0_idx].x = (ptcl_grid->points[temp_v1_idx].x * l_2 + ptcl_grid->points[temp_v2_idx].x * l_1) / (l_1 + l_2);
                        ptcl_grid->points[temp_v0_idx].y = (ptcl_grid->points[temp_v1_idx].y * l_2 + ptcl_grid->points[temp_v2_idx].y * l_1) / (l_1 + l_2);
                        ptcl_grid->points[temp_v0_idx].z = (ptcl_grid->points[temp_v1_idx].z * l_2 + ptcl_grid->points[temp_v2_idx].z * l_1) / (l_1 + l_2);
                    }
                }
                if (uniform_grid[i-1][j-1].cell_level == max_level && uniform_grid[i-1][j].cell_level == max_level)
                {
                    bool if_need_another = false;
                    int l_1 = 1, l_2 = 1;
                    while (j - l_1 >= 0)
                    {
                        if (uniform_grid[i-1][j-l_1-1].cell_level != max_level) {
                            if_need_another = true;
                            break;
                        }
                        if (grid_vertex_count[i][j-l_1] == 4) { break; }
                        l_1++;
                    }
                    while (j + l_2 < temp_W)
                    {
                        if (uniform_grid[i-1][j+l_2].cell_level != max_level) {
                            if_need_another = true;
                            break;
                        }
                        if (grid_vertex_count[i][j+l_2] == 4) { break; }
                        l_2++;
                    }

                    if ((l_1 + l_2) == std::pow(2, max_level-1))
                    {
                        int temp_v0_idx = tri_vertex_idx_grid[i*2][j*2];
                        int temp_v1_idx = tri_vertex_idx_grid[i*2][j*2-2*l_1];
                        int temp_v2_idx = tri_vertex_idx_grid[i*2][j*2+2*l_2];

                        if (if_need_another)
                        {
                            another_adjust.push_back(temp_v0_idx);
                            another_adjust.push_back(temp_v1_idx);
                            another_adjust.push_back(temp_v2_idx);
                            another_adjust.push_back(l_1);
                            another_adjust.push_back(l_2);
                        }

                        ptcl_grid->points[temp_v0_idx].x = (ptcl_grid->points[temp_v1_idx].x * l_2 + ptcl_grid->points[temp_v2_idx].x * l_1) / (l_1 + l_2);
                        ptcl_grid->points[temp_v0_idx].y = (ptcl_grid->points[temp_v1_idx].y * l_2 + ptcl_grid->points[temp_v2_idx].y * l_1) / (l_1 + l_2);
                        ptcl_grid->points[temp_v0_idx].z = (ptcl_grid->points[temp_v1_idx].z * l_2 + ptcl_grid->points[temp_v2_idx].z * l_1) / (l_1 + l_2);
                    }
                }
            }
        }

        // 二次调整
        for (int i = 0; i < another_adjust.size(); i+=5)
        {
            int temp_v0_idx = another_adjust[i];
            int temp_v1_idx = another_adjust[i+1];
            int temp_v2_idx = another_adjust[i+2];
            int l_1 = another_adjust[i+3];
            int l_2 = another_adjust[i+4];

            ptcl_grid->points[temp_v0_idx].x = (ptcl_grid->points[temp_v1_idx].x * l_2 + ptcl_grid->points[temp_v2_idx].x * l_1) / (l_1 + l_2);
            ptcl_grid->points[temp_v0_idx].y = (ptcl_grid->points[temp_v1_idx].y * l_2 + ptcl_grid->points[temp_v2_idx].y * l_1) / (l_1 + l_2);
            ptcl_grid->points[temp_v0_idx].z = (ptcl_grid->points[temp_v1_idx].z * l_2 + ptcl_grid->points[temp_v2_idx].z * l_1) / (l_1 + l_2);
        }
    }
}

// 融合两个mesh fragment的顶点列表、面片列表
void MeshFragment::vertex_and_face_list_merge(const std::shared_ptr<MeshFragment>& new_plane)
{
    int last_pts_num = ptcl_grid->size();
    int last_tri_num = faces_list.size();
    *ptcl_grid += *(new_plane->ptcl_grid);

    for (int i = 0; i < new_plane->faces_list.size(); i++)
    {
        std::array<int, 3> new_tri_idx;
        new_tri_idx[0] = new_plane->faces_list[i][0] + last_pts_num;
        new_tri_idx[1] = new_plane->faces_list[i][1] + last_pts_num;
        new_tri_idx[2] = new_plane->faces_list[i][2] + last_pts_num;
        faces_list.push_back(new_tri_idx);
    }

    for (int i = 0; i < new_plane->faces_to_delete.size(); i++)
    {
        int delete_tri_idx = new_plane->faces_to_delete[i];
        faces_to_delete.push_back(delete_tri_idx + last_tri_num);
    }
}

// 将mesh保存为.ply格式的文件，用pcl库
bool MeshFragment::save_to_ply(const std::string &fileName)
{
    pcl::PolygonMesh mesh;
    mesh.header = pcl::PCLHeader();
    mesh.cloud = pcl::PCLPointCloud2();
    mesh.polygons = std::vector<pcl::Vertices>();

    pcl::toPCLPointCloud2(*ptcl_grid, mesh.cloud);

    std::sort(faces_to_delete.begin(), faces_to_delete.end());   // 从小到大排序
    int delete_i = 0;
    for (int i = 0; i < faces_list.size(); i++)
    {
        if (faces_to_delete.size() > 0)
        {
            if (i == faces_to_delete[delete_i])
            {
                while (delete_i != faces_to_delete.size() - 1 && faces_to_delete[delete_i] <= i)
                {
                    delete_i++;
                }
//                assert(faces_to_delete[delete_i] > i);
                continue;
            }
        }

        pcl::Vertices vertices;
        vertices.vertices.push_back(faces_list[i][0]);   // meshlab点序逆时针
        vertices.vertices.push_back(faces_list[i][1]);
        vertices.vertices.push_back(faces_list[i][2]);
        mesh.polygons.push_back(vertices);
    }

//    pcl::io::savePLYFile(fileName, mesh);
    pcl::io::savePLYFileBinary(fileName, mesh);
    std::cout << "Planar mesh saved to " << fileName << std::endl;

    return true;
}

// 保存时删掉重复的vertex
bool MeshFragment::save_to_ply_without_redundancy(const std::string &fileName)
{
    pcl::PolygonMesh mesh;
    mesh.header = pcl::PCLHeader();
    mesh.cloud = pcl::PCLPointCloud2();
    mesh.polygons = std::vector<pcl::Vertices>();

    std::vector<int> face_exist;
    std::sort(faces_to_delete.begin(), faces_to_delete.end());   // 从小到大排序
    int delete_i = 0;
    for (int i = 0; i < faces_list.size(); i++)
    {
        if (faces_to_delete.size() > 0)
        {
            if (i == faces_to_delete[delete_i])
            {
                while (delete_i != faces_to_delete.size() - 1 && faces_to_delete[delete_i] <= i)
                {
                    delete_i++;
                }
//                assert(faces_to_delete[delete_i] > i);
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
        vertex_exist_single.insert(faces_list[it][0]);
        vertex_exist_single.insert(faces_list[it][1]);
        vertex_exist_single.insert(faces_list[it][2]);
    }

    int i_ = 0;
    for (auto it : vertex_exist_single)
    {
        pcl::PointXYZ temp_pt;
        temp_pt.x = ptcl_grid->points[it].x;
        temp_pt.y = ptcl_grid->points[it].y;
        temp_pt.z = ptcl_grid->points[it].z;

        ptcl_processed->push_back(temp_pt);
        vertex_idx_old2new[it] = i_;
        i_++;
    }
    pcl::toPCLPointCloud2(*ptcl_processed, mesh.cloud);

    for (auto it : face_exist)
    {
        int x_ = vertex_idx_old2new[faces_list[it][0]];
        int y_ = vertex_idx_old2new[faces_list[it][1]];
        int z_ = vertex_idx_old2new[faces_list[it][2]];

        pcl::Vertices vertices;
        vertices.vertices.push_back(x_);   // meshlab点序逆时针
        vertices.vertices.push_back(y_);
        vertices.vertices.push_back(z_);
        mesh.polygons.push_back(vertices);
    }

//    pcl::io::savePLYFile(fileName, mesh);
    pcl::io::savePLYFileBinary(fileName, mesh);
    std::cout << "Plane mesh saved to " << fileName << std::endl;

    return true;
}

































