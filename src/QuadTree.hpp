//
// Created by neo on 2024/9/23.
//

#ifndef SRC_QUADTREE_H
#define SRC_QUADTREE_H
#include <iostream>
#include <string>
// Eigen
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

// 均匀网格
struct grid {
    int             pts_num;                        // 点数量
    int             node_type;                      // 类型 0 1 2 外部 内部 边缘
    double          x_avg;                          // 沿平面法向的偏移距离
    bool            if_grid_need_update;            // 有新点加入时需更新网格状态、顶点状态、顶点列表
    bool            if_mesh_need_delete;            // 原来的顶点、三角形需要删除
    bool            if_mesh_need_add;               // 需要加入新的点、三角形

    std::vector<int>                tri_idx_list;           // 属于该网格的三角形索引
    std::vector<int>                ver_idx_list;           // 内部点，四个顶点索引
    Eigen::Vector4i                 grid_vertex_state;      // 网格顶点状态 1-内部 0-外部
    Eigen::Vector2d                 grid_center;            // yAvg zAvg
    Eigen::Matrix<double, 4, 2>     grid_vertex_list;       // y z 左上角开始 顺时针
    Eigen::Matrix<double, 4, 2>     edge_vertex;            // yMax_n, yMin_n, zMin_n, zMax_n, yMin_0, yMax_0, zMax_0, zMin_0
    std::vector<Eigen::Vector2d>    vertex_list;            // local坐标系下 y z
    int                             cell_level;             // 用于四叉树简化中调整边上的顶点的x值
};

class QuadTreeNode {
public:
    int                 type;           // 节点类型
    bool                if_leaf_node;   // 是否叶子节点
    int                 node_umin;      // 节点在grid中的范围
    int                 node_umax;      //
    int                 node_vmin;      //
    int                 node_vmax;       //
    std::vector<int>    vertex_list;    // 顶点列表
    std::vector<int>    face_list;      // 三角形索引列表
    QuadTreeNode*       topLeft;        // 子节点指针
    QuadTreeNode*       topRight;
    QuadTreeNode*       bottomLeft;
    QuadTreeNode*       bottomRight;

    QuadTreeNode(int t, bool s)
    {
        type = t;
        if_leaf_node = s;
        topLeft = nullptr;
        topRight = nullptr;
        bottomLeft = nullptr;
        bottomRight = nullptr;
    };
    ~QuadTreeNode() {};

    // 判断四个子节点是否类型一致
    bool isUniform() {
        if (!topLeft || !topRight || !bottomLeft || !bottomRight) {
            return false; // 如果不是四个子节点都有值，则不统一
        }
        return topLeft->type == 1 &&
               topRight->type == 1 &&
               bottomLeft->type == 1 &&
               bottomRight->type == 1;
    }

    // 合并四个子节点
    void mergeChildren() {

        delete topLeft;
        delete topRight;
        delete bottomLeft;
        delete bottomRight;
        topLeft = nullptr;
        topRight = nullptr;
        bottomLeft = nullptr;
        bottomRight = nullptr;
    }
};


class QuadTree {
public:
    QuadTreeNode*                       Root;                   // 根节点
    std::vector<int>                    faces_idx_to_delete;    // 要删除的三角形
    std::vector<std::array<int, 3>>     faces_to_add;           // 要添加的三角形
    std::deque<std::deque<grid>>*       uniform_grid;
    std::vector<std::vector<int>>*      grid_vertex_count;
    std::deque<std::deque<int>>*        tri_vertex_idx_grid;
    int H, W, H_start, W_start, extend_side;

    QuadTree()
    {
        Root = nullptr;
    }

    // 递归建立四叉树子节点
    void build_ChildNode(QuadTreeNode* node, int u_min, int v_min, int u_max, int v_max)
    {
        if ( (u_max - u_min) == 0 && (v_max - v_min) == 0 )
        {
            node->if_leaf_node = true;
            if (u_min >= H_start && u_min < (H_start + H) && v_min >= W_start && v_min < (W_start + W))
            {
                if ((*uniform_grid)[u_min-H_start][v_min-W_start].node_type == 1)
                {
                    node->type = 1;
                    // 添加三角形列表 - 为了删除
                    for (auto it : (*uniform_grid)[u_min-H_start][v_min-W_start].tri_idx_list)
                    {
                        node->face_list.push_back(it);
                    }
                    // 添加顶点列表，顺时针，左上开始 - 为了提取合并后的顶点
                    for (auto it : (*uniform_grid)[u_min-H_start][v_min-W_start].ver_idx_list)
                    {
                        node->vertex_list.push_back(it);
                    }
                }
                if ((*uniform_grid)[u_min-H_start][v_min-W_start].node_type == 2)
                {
                    node->type = 2;
                }
                node->node_umin = u_min-H_start;
                node->node_umax = u_min-H_start;
                node->node_vmin = v_min-W_start;
                node->node_vmax = v_min-W_start;
            }

            return;
        }

        node->topLeft = new QuadTreeNode(0, false);
        node->topRight = new QuadTreeNode(0, false);
        node->bottomLeft = new QuadTreeNode(0, false);
        node->bottomRight = new QuadTreeNode(0, false);

        int u_mid = int((u_min + u_max) / 2);
        int v_mid = int((v_min + v_max) / 2);
        build_ChildNode(node->topLeft, u_min, v_min, u_mid, v_mid);
        build_ChildNode(node->topRight, u_min, v_mid+1, u_mid, v_max);
        build_ChildNode(node->bottomLeft, u_mid+1, v_min, u_max, v_mid);
        build_ChildNode(node->bottomRight, u_mid+1, v_mid+1, u_max, v_max);
    }


    // 建立四叉树
    void build()
    {
        H = (*uniform_grid).size();
        W = (*uniform_grid)[0].size();
        int large_side = std::max(H, W);
        int log2n = int(std::log2(large_side)) + 1;
        extend_side = int(std::pow(2.0, log2n));
        H_start = int((extend_side - H) / 2);
        W_start = int((extend_side - W) / 2);

        Root = new QuadTreeNode(0, false);
        build_ChildNode(Root, 0, 0, extend_side-1, extend_side-1);
    }

    // 递归合并四个子节点的函数
    void mergeQuadTree(QuadTreeNode* node) {
        if (node == nullptr) return;

        if (node->topLeft) mergeQuadTree(node->topLeft);
        if (node->topRight) mergeQuadTree(node->topRight);
        if (node->bottomLeft) mergeQuadTree(node->bottomLeft);
        if (node->bottomRight) mergeQuadTree(node->bottomRight);

        if (node->isUniform()) {
            node->type = 1;
            node->node_umin = node->topLeft->node_umin;
            node->node_umax = node->bottomRight->node_umax;
            node->node_vmin = node->topLeft->node_vmin;
            node->node_vmax = node->bottomRight->node_vmax;
            for (int i = node->node_umin; i <= node->node_umax; i++)
            {
                for (int j = node->node_vmin; j <= node->node_vmax; j++)
                {
                    (*uniform_grid)[i][j].cell_level++;
                }
            }

            // 删除的三角形
            for (auto it : node->topLeft->face_list) {
                faces_idx_to_delete.push_back(it);
            }
            for (auto it : node->topRight->face_list) {
                faces_idx_to_delete.push_back(it);
            }
            for (auto it : node->bottomLeft->face_list) {
                faces_idx_to_delete.push_back(it);
            }
            for (auto it : node->bottomRight->face_list) {
                faces_idx_to_delete.push_back(it);
            }

            // 新的三角形顶点列表
            node->vertex_list.push_back(node->topLeft->vertex_list[0]);
            node->vertex_list.push_back(node->topRight->vertex_list[1]);
            node->vertex_list.push_back(node->bottomRight->vertex_list[2]);
            node->vertex_list.push_back(node->bottomLeft->vertex_list[3]);

            delete node->topLeft;
            delete node->topRight;
            delete node->bottomLeft;
            delete node->bottomRight;
            node->topLeft = nullptr;
            node->topRight = nullptr;
            node->bottomLeft = nullptr;
            node->bottomRight = nullptr;
        }
    }

    // 提取合并节点的三角形
    void faces_extract(QuadTreeNode* node)
    {
        if (node == nullptr) return;

        if (node->topLeft) faces_extract(node->topLeft);
        if (node->topRight) faces_extract(node->topRight);
        if (node->bottomLeft) faces_extract(node->bottomLeft);
        if (node->bottomRight) faces_extract(node->bottomRight);

        if (!node->if_leaf_node && node->vertex_list.size() != 0)
        {
            std::array<int, 3> temp_tri_0, temp_tri_1;
            temp_tri_0[0] = node->vertex_list[0];
            temp_tri_0[1] = node->vertex_list[2];
            temp_tri_0[2] = node->vertex_list[1];

            temp_tri_1[0] = node->vertex_list[2];
            temp_tri_1[1] = node->vertex_list[0];
            temp_tri_1[2] = node->vertex_list[3];

            faces_to_add.push_back(temp_tri_0);
            faces_to_add.push_back(temp_tri_1);
        }
    }

    // 提取简化后保留的grid顶点，以矩阵形式保存
    void grid_vertex_extract(QuadTreeNode* node)
    {
        if (node == nullptr) return;

        if (node->topLeft) grid_vertex_extract(node->topLeft);
        if (node->topRight) grid_vertex_extract(node->topRight);
        if (node->bottomLeft) grid_vertex_extract(node->bottomLeft);
        if (node->bottomRight) grid_vertex_extract(node->bottomRight);

        if (!node->if_leaf_node && node->vertex_list.size() != 0)
        {
            (*grid_vertex_count)[node->node_umin][node->node_vmin]++;
            (*grid_vertex_count)[node->node_umin][node->node_vmax+1]++;
            (*grid_vertex_count)[node->node_umax+1][node->node_vmin]++;
            (*grid_vertex_count)[node->node_umax+1][node->node_vmax+1]++;
        }
        if (node->if_leaf_node && node->type == 1)
        {
            (*grid_vertex_count)[node->node_umin][node->node_vmin]++;
            (*grid_vertex_count)[node->node_umin][node->node_vmax+1]++;
            (*grid_vertex_count)[node->node_umax+1][node->node_vmin]++;
            (*grid_vertex_count)[node->node_umax+1][node->node_vmax+1]++;
        }
        if (node->if_leaf_node && node->type == 2)
        {
            if ((*tri_vertex_idx_grid)[2*node->node_umin][2*node->node_vmin] != -1)
            {
                (*grid_vertex_count)[node->node_umin][node->node_vmin]++;
            }
            if ((*tri_vertex_idx_grid)[2*node->node_umin][2*node->node_vmin+2] != -1)
            {
                (*grid_vertex_count)[node->node_umin][node->node_vmin+1]++;
            }
            if ((*tri_vertex_idx_grid)[2*node->node_umin+2][2*node->node_vmin] != -1)
            {
                (*grid_vertex_count)[node->node_umin+1][node->node_vmin]++;
            }
            if ((*tri_vertex_idx_grid)[2*node->node_umin+2][2*node->node_vmin+2] != -1)
            {
                (*grid_vertex_count)[node->node_umin+1][node->node_vmin+1]++;
            }
        }
    }

    // 清理四叉树内存
    void clearTree(QuadTreeNode* node) {
        if (node == nullptr) return;
        clearTree(node->topLeft);
        clearTree(node->topRight);
        clearTree(node->bottomLeft);
        clearTree(node->bottomRight);
        delete node;
    }
};

#endif //SRC_QUADTREE_H
