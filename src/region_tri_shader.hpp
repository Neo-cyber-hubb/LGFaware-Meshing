//
// Created by neo on 2024/7/10.
//

#ifndef SRC_REGION_TRI_SHADER_HPP
#define SRC_REGION_TRI_SHADER_HPP
#include "MeshFragment.h"
#include "NonPlaneMesh.h"

extern vec_3f        g_axis_min_max[ 2 ];
extern std::unordered_map<int, std::shared_ptr<MeshFragment>> plane_idx2ptr;
extern std::vector<std::shared_ptr<MeshFragment>> plane_map;         // 各不相同的plane
extern std::shared_ptr<NonPlaneMesh> non_plane_map;

struct Region_triangles_shader
{
    std::vector< vec_3f >               m_triangle_pt_vec;     // 三角形顶点列表，有重复
    Common_tools::Triangle_facet_shader m_triangle_facet_shader;
    int                                 m_need_init_shader = true;
    int                                 m_need_refresh_shader = true;
    int                                 m_if_set_color = false;
    std::shared_ptr< std::mutex >       m_mutex_ptr = nullptr;

    Region_triangles_shader() { m_mutex_ptr = std::make_shared< std::mutex >(); }

    void init_openGL_shader() { m_triangle_facet_shader.init( SHADER_DIR ); }

    Common_tools::Triangle_facet_shader *get_shader_ptr() { return &m_triangle_facet_shader; }

    void init_pointcloud()
    {
        std::unique_lock< std::mutex > lock( *m_mutex_ptr );
        if ( m_if_set_color )
        {
            m_triangle_facet_shader.set_pointcloud( m_triangle_pt_vec, g_axis_min_max, 2 );
        }
        else
        {
            m_triangle_facet_shader.set_pointcloud( m_triangle_pt_vec );
        }
    }

    // 取三角形顶点，应该存在重复吧
    void unparse_triangle_set_to_vector( const Triangle_set &tri_angle_set )
    {
//        // TODO: synchronized data buffer here:
        std::unique_lock< std::mutex > lock( *m_mutex_ptr );
        m_triangle_pt_vec.resize( tri_angle_set.size() * 3 );
//        // cout << "Number of pt_size = " << m_triangle_pt_list.size() << endl;
        int count = 0;
        for ( Triangle_set::iterator it = tri_angle_set.begin(); it != tri_angle_set.end(); it++ )
        {
            vec_3f pt_a, pt_b, pt_c;
            pt_a << 0, 0, 0;
            pt_b << 0, 0, 0;
            pt_c << 0, 0, 0;

            if (( *it )->rec_source == 0)
            {
                non_plane_map->mtx_pts_list.lock();
                pt_a << non_plane_map->pts_list->points[( *it )->m_tri_pts_id[ 0 ]].x,
                        non_plane_map->pts_list->points[( *it )->m_tri_pts_id[ 0 ]].y,
                        non_plane_map->pts_list->points[( *it )->m_tri_pts_id[ 0 ]].z;
                pt_b << non_plane_map->pts_list->points[( *it )->m_tri_pts_id[ 1 ]].x,
                        non_plane_map->pts_list->points[( *it )->m_tri_pts_id[ 1 ]].y,
                        non_plane_map->pts_list->points[( *it )->m_tri_pts_id[ 1 ]].z;
                pt_c << non_plane_map->pts_list->points[( *it )->m_tri_pts_id[ 2 ]].x,
                        non_plane_map->pts_list->points[( *it )->m_tri_pts_id[ 2 ]].y,
                        non_plane_map->pts_list->points[( *it )->m_tri_pts_id[ 2 ]].z;
                non_plane_map->mtx_pts_list.unlock();
            }
            else
            {
                auto temp_it = plane_idx2ptr.find(( *it )->rec_source);
                if (temp_it == plane_idx2ptr.end()) { continue; }

                auto temp_mf = plane_idx2ptr[( *it )->rec_source];
                temp_mf->mtx_ptcl_grid.lock();
                pt_a << temp_mf->ptcl_grid->points[( *it )->m_tri_pts_id[ 0 ]].x,
                        temp_mf->ptcl_grid->points[( *it )->m_tri_pts_id[ 0 ]].y,
                        temp_mf->ptcl_grid->points[( *it )->m_tri_pts_id[ 0 ]].z;
                pt_b << temp_mf->ptcl_grid->points[( *it )->m_tri_pts_id[ 1 ]].x,
                        temp_mf->ptcl_grid->points[( *it )->m_tri_pts_id[ 1 ]].y,
                        temp_mf->ptcl_grid->points[( *it )->m_tri_pts_id[ 1 ]].z;
                pt_c << temp_mf->ptcl_grid->points[( *it )->m_tri_pts_id[ 2 ]].x,
                        temp_mf->ptcl_grid->points[( *it )->m_tri_pts_id[ 2 ]].y,
                        temp_mf->ptcl_grid->points[( *it )->m_tri_pts_id[ 2 ]].z;
                temp_mf->mtx_ptcl_grid.unlock();
            }

            m_triangle_pt_vec[ count ] = pt_a.cast< float >();
            m_triangle_pt_vec[ count + 1 ] = pt_b.cast< float >();
            m_triangle_pt_vec[ count + 2 ] = pt_c.cast< float >();
            count = count + 3;
        }
        m_triangle_pt_vec.resize( count );
    }

    void get_axis_min_max( vec_3f *axis_min_max = nullptr )
    {
        if ( axis_min_max != nullptr )
        {
            for ( int i = 0; i < m_triangle_pt_vec.size(); i++ )
            {
                if ( axis_min_max[ 0 ]( 0 ) > m_triangle_pt_vec[ i ]( 0 ) )
                {
                    axis_min_max[ 0 ]( 0 ) = m_triangle_pt_vec[ i ]( 0 );
                }
                if ( axis_min_max[ 0 ]( 1 ) > m_triangle_pt_vec[ i ]( 1 ) )
                {
                    axis_min_max[ 0 ]( 1 ) = m_triangle_pt_vec[ i ]( 1 );
                }
                if ( axis_min_max[ 0 ]( 2 ) > m_triangle_pt_vec[ i ]( 2 ) )
                {
                    axis_min_max[ 0 ]( 2 ) = m_triangle_pt_vec[ i ]( 2 );
                }
                if ( axis_min_max[ 1 ]( 0 ) < m_triangle_pt_vec[ i ]( 0 ) )
                {
                    axis_min_max[ 1 ]( 0 ) = m_triangle_pt_vec[ i ]( 0 );
                }
                if ( axis_min_max[ 1 ]( 1 ) < m_triangle_pt_vec[ i ]( 1 ) )
                {
                    axis_min_max[ 1 ]( 1 ) = m_triangle_pt_vec[ i ]( 1 );
                }
                if ( axis_min_max[ 1 ]( 2 ) < m_triangle_pt_vec[ i ]( 2 ) )
                {
                    axis_min_max[ 1 ]( 2 ) = m_triangle_pt_vec[ i ]( 2 );
                }
            }
        }
    }

    void synchronized_from_region( Sync_triangle_set *sync_triangle_set, vec_3f *axis_min_max = nullptr )
    {
        if ( sync_triangle_set == nullptr )
        {
            cout << "sync_triangle_set == nullptr" << endl;
            return;
        }

        if ( sync_triangle_set->m_if_required_synchronized )
        {
            Triangle_set triangle_set;
            sync_triangle_set->get_triangle_set( triangle_set, true );
            unparse_triangle_set_to_vector( triangle_set );   // 获取三角形顶点
            get_axis_min_max( axis_min_max );   // 根据点的坐标更新矩形框范围
            std::this_thread::sleep_for( std::chrono::microseconds( 100 ) );   // ？ 100微秒
            m_need_refresh_shader = true;  // 有新三角形更新shader
        }
    }

    void draw( const Cam_view &gl_cam )
    {
        if ( m_need_init_shader )  // 初始化shader
        {
            init_openGL_shader();
            m_need_init_shader = false;
        }
        if ( m_triangle_pt_vec.size() < 3 )
        {
            return;
        }
        if ( m_need_refresh_shader )  // 更新shader（点列表）
        {
            init_pointcloud();
            m_need_refresh_shader = false;
        }
        m_triangle_facet_shader.draw( gl_cam.m_glm_projection_mat, Common_tools::eigen2glm( gl_cam.m_camera_pose_mat44_inverse ) );
    }
};

#endif //SRC_REGION_TRI_SHADER_HPP
