//
// Created by neo on 2024/7/10.
//

#include "GL_gui.h"

//extern vec_3f        g_axis_min_max[ 2 ];

void GL_gui::print_help_window( bool *if_display_help_win )
{
    ImGui::Begin( "--- Help ---", if_display_help_win );
    ImGui::Text( "[H]     | Display/Close main windows" );
    ImGui::Text( "[C]     | Show/Close camera pose windows" );
    ImGui::Text( "[T]     | Follow the camera" );
    ImGui::Text( "[D]     | Show/close the depth image" );
    ImGui::Text( "[L]     | Show/close the LiDAR points" );
    ImGui::Text( "[M]     | Show/close the mesh" );
    ImGui::Text( "[S]     | Save camera view" );
    ImGui::Text( "[Z]     | Load camera view" );
    ImGui::Text( "[+/-]   | Increase/Decrease the line width" );
    ImGui::Text( "[F1]   | Display help window" );
    ImGui::Text( "[Space] | To pause the program" );
    ImGui::Text( "[Esc]   | Exit the program" );
    ImGui::End();
}

void GL_gui::get_last_avr_pose( int current_frame_idx, Eigen::Quaterniond &q_avr, vec_3 &t_vec )
{
    const int win_ssd = 1;
    mat_3_3   lidar_frame_to_camera_frame;
    // Clang-format off
    lidar_frame_to_camera_frame << 0, 0, -1, -1, 0, 0, 0, 1, 0;
    // Clang-format on
    q_avr = Eigen::Quaterniond::Identity();
    t_vec = vec_3::Zero();
    if ( current_frame_idx < 1 )
    {
        return;
    }
    int                frame_count = 0;
    int                frame_s = std::max( 0, current_frame_idx - win_ssd );
    vec_3              log_angle_acc = vec_3( 0, 0, 0 );
    Eigen::Quaterniond q_first;
    for ( int frame_idx = frame_s; frame_idx < current_frame_idx; frame_idx++ )
    {
        if ( g_eigen_vec_vec[ frame_idx ].second.size() != 0 )
        {
            Eigen::Quaterniond pose_q( g_eigen_vec_vec[ frame_idx ].second.head< 4 >() );
            pose_q.normalize();
            if ( frame_count == 0 )
            {
                q_first = pose_q;
            }
            q_avr = q_avr * pose_q;
            log_angle_acc += Sophus::SO3d( q_first.inverse() * pose_q ).log();
            t_vec = t_vec + g_eigen_vec_vec[ frame_idx ].second.block( 4, 0, 3, 1 );
            frame_count++;
        }
    }
    t_vec = t_vec / frame_count;
    q_avr = q_first * Sophus::SO3d::exp( log_angle_acc / frame_count ).unit_quaternion();
    q_avr.normalize();
    q_avr = q_avr * Eigen::Quaterniond( lidar_frame_to_camera_frame );
}

void GL_gui::display_current_LiDAR_pts( int current_frame_idx, double pts_size, vec_4f color )
{
//    if ( current_frame_idx < 1 )
//    {
//        return;
//    }
//    g_LiDAR_point_shader.set_point_attr( pts_size );
//    g_LiDAR_point_shader.set_pointcloud( g_eigen_vec_vec[ current_frame_idx ].first, vec_3( 1.0, 1.0, 1.0 ) );
//    g_LiDAR_point_shader.draw( g_gl_camera.m_gl_cam.m_glm_projection_mat,
//                               Common_tools::eigen2glm( g_gl_camera.m_gl_cam.m_camera_pose_mat44_inverse ) );
}

void GL_gui::display_reinforced_LiDAR_pts( std::vector< vec_3f > &pt_vec, double pts_size, vec_3f color )
{
//    g_LiDAR_point_shader.set_point_attr( pts_size );
//    g_LiDAR_point_shader.set_pointcloud( pt_vec, color.cast< double >() );
//    g_LiDAR_point_shader.draw( g_gl_camera.m_gl_cam.m_glm_projection_mat,
//                               Common_tools::eigen2glm( g_gl_camera.m_gl_cam.m_camera_pose_mat44_inverse ) );
}

void GL_gui::init_openGL_shader()
{
    g_LiDAR_point_shader.init( SHADER_DIR );
    g_path_shader.init( SHADER_DIR );
    // Init axis buffer
    g_axis_shader.init( SHADER_DIR, 1 );
    // Init ground shader
    g_ground_plane_shader.init( SHADER_DIR, 10, 10 );
    g_camera_pose_shader.init( SHADER_DIR );
    g_triangle_facet_shader.init( SHADER_DIR );
}

// gui显示三角形
void GL_gui::draw_triangle( const Cam_view &gl_cam )
{
    int region_size = g_region_triangles_shader_vec.size();  // 遍历所有的shader指针
    for ( int region_idx = 0; region_idx < region_size; region_idx++ )
    {
        g_region_triangles_shader_vec[ region_idx ]->m_triangle_facet_shader.m_if_draw_face = g_display_face; // 显示面 还是 线框
        g_region_triangles_shader_vec[ region_idx ]->m_if_set_color = g_mesh_if_color;
        g_region_triangles_shader_vec[ region_idx ]->draw( gl_cam );
    }
}

void GL_gui::display_camera_traj( float display_size )
{
    if ( pt_camera_traj.size() == 0 )
    {
        return;
    }
    vec_3 pt_color = vec_3( 0.0, 1.0, 0.0 );
    g_path_shader.set_pointcloud( pt_camera_traj, pt_color );
    g_path_shader.set_point_attr( display_size + 2, 0, 1.0 );
    g_path_shader.m_draw_points_number = pt_camera_traj.size();
    g_path_shader.draw( g_gl_camera.m_gl_cam.m_glm_projection_mat, Common_tools::eigen2glm( g_gl_camera.m_gl_cam.m_camera_pose_mat44_inverse ),
                        GL_LINE_STRIP );
}

void GL_gui::draw_camera_pose( int current_frame_idx, float pt_disp_size, float display_cam_size )
{

    Eigen::Quaterniond pose_q( g_eigen_vec_vec[ current_frame_idx ].second.head< 4 >() );
    vec_3              pose_t = g_eigen_vec_vec[ current_frame_idx ].second.block( 4, 0, 3, 1 );
    mat_3_3            lidar_frame_to_camera_frame;
    lidar_frame_to_camera_frame << 0, 0, 1, -1, 0, 0, 0, -1, 0;
    pose_q = Eigen::Quaterniond( pose_q.toRotationMatrix() * lidar_frame_to_camera_frame );

    pose_t = pose_q.inverse() * ( pose_t * -1.0 );
    pose_q = pose_q.inverse();

    g_camera_pose_shader.set_camera_pose_and_scale( pose_q, pose_t, display_cam_size );
    g_camera_pose_shader.set_point_attr( 5, 0, 1.0 );
    g_camera_pose_shader.draw( g_gl_camera.m_gl_cam.m_glm_projection_mat, Common_tools::eigen2glm( g_gl_camera.m_gl_cam.m_camera_pose_mat44_inverse ),
                               -1 );
}

void GL_gui::draw_camera_trajectory( int current_frame_idx, float pt_disp_size )
{
    pt_camera_traj.clear();
    for ( int i = 0; i < current_frame_idx; i++ )
    {
        if ( g_eigen_vec_vec[ i ].second.size() >= 7 )
        {
            pt_camera_traj.push_back( g_eigen_vec_vec[ i ].second.block( 4, 0, 3, 1 ) );
        }
    }
    display_camera_traj( pt_disp_size );
}

void GL_gui::display()
{
    // Setup window
    // PCL库将设置其日志系统，以确保所有信息都会被输出到控制台，包括调试信息、警告和错误信息
//    pcl::console::setVerbosityLevel( pcl::console::L_ALWAYS );
//    Common_tools::printf_software_version();
//    printf_program( "ImMesh: An Immediate LiDAR Localization and Meshing Framework" );

//    ros::init( argc, argv, "laserMapping" );
//    voxel_mapping.init_ros_node();

    GLFWwindow *window = g_gl_camera.init_openGL_and_ImGUI( "Mesh visualization", 1, 20 );
    // gladLoadGLLoader 函数用于从适当的源（例如，通过glfwGetProcAddress获取的函数指针）加载OpenGL函数指针。一旦函数指针被加载，就可以像使用普通函数一样使用OpenGL的API。
    if ( !gladLoadGLLoader( ( GLADloadproc ) glfwGetProcAddress ) )
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
    }
    Common_tools::Point_cloud_shader  g_pt_shader;
    init_openGL_shader();
    if ( window == nullptr )
    {
        cout << "Window == nullptr" << endl;
    }

    // 从读取的点云文件中重建mesh，跳过
//    g_enable_mesh_rec = voxel_mapping.m_if_enable_mesh_rec;
//    cout << "Offline point cloud name: " << ANSI_COLOR_GREEN_BOLD << voxel_mapping.m_pointcloud_file_name << ANSI_COLOR_RESET << endl;
//    if ( Common_tools::if_file_exist( voxel_mapping.m_pointcloud_file_name ) )
//    {
//        pcl::PointCloud< pcl::PointXYZI > offline_pts;
//        cout << "Loading data..." ;
//        fflush( stdout );
//        pcl::io::loadPCDFile( voxel_mapping.m_pointcloud_file_name, offline_pts );
//        cout << " total of pts = " << offline_pts.points.size() << endl;
//        cout << "g_map_rgb_pts_mesh.m_minimum_pts_size = " << g_map_rgb_pts_mesh.m_minimum_pts_size << endl;
//        reconstruct_mesh_from_pointcloud( offline_pts.makeShared() );
//    }
//    else if(voxel_mapping.m_pointcloud_file_name.length() > 5)
//    {
//        cout << ANSI_COLOR_RED_BOLD << "Offline point cloud file: " << voxel_mapping.m_pointcloud_file_name <<" NOT exist!!!, Please check!!!" << ANSI_COLOR_RESET << endl;
//        while(1);
//    }


//    cout << "====Loading parameter=====" << endl;
//
//    threshold_scale = voxel_mapping.m_meshing_distance_scale;
//    minimum_pts = voxel_mapping.m_meshing_points_minimum_scale * voxel_mapping.m_meshing_distance_scale;
//    g_meshing_voxel_size = voxel_mapping.m_meshing_voxel_resolution * voxel_mapping.m_meshing_distance_scale;
//    appending_pts_frame = voxel_mapping.m_meshing_number_of_pts_append_to_map;
//    region_size = voxel_mapping.m_meshing_region_size * voxel_mapping.m_meshing_distance_scale;
//    g_display_mesh = voxel_mapping.m_if_draw_mesh;
//    scope_color( ANSI_COLOR_YELLOW_BOLD );
//    cout << "=========Meshing config ========= " << endl;
//    cout << "Threshold scale = " << threshold_scale << endl;    // 1
//    cout << "Minimum pts distance = " << minimum_pts << endl;   // 0.1
//    cout << "Voxel size = " << g_meshing_voxel_size << endl;    // 0.4
//    cout << "Region size = " << region_size << endl;            // 10

    g_current_frame = -3e8;
//    g_triangles_manager.m_pointcloud_map = &g_map_rgb_pts_mesh;
//    g_map_rgb_pts_mesh.sest_minimum_dis( minimum_pts );
//    g_map_rgb_pts_mesh.set_voxel_resolution( g_meshing_voxel_size );
    g_triangles_manager.m_region_size = 10;
//    g_map_rgb_pts_mesh.m_recent_visited_voxel_activated_time = 0;
//    cout << "==== Loading parameter end =====" << endl;

    // 创建一个新的线程，并将其与线程函数关联起来
    // 第一个参数 &Voxel_mapping::service_LiDAR_update 是一个函数的地址，这个函数是 Voxel_mapping 类的成员函数 service_LiDAR_update。& 符号表示取地址操作符。
    // 第二个参数 &voxel_mapping 是一个指向 Voxel_mapping 类实例的指针。这个实例将作为 service_LiDAR_update 成员函数的 this 指针，即调用上下文中的对象。

//    std::thread thr_mapping = std::thread( &Voxel_mapping::service_LiDAR_update, &voxel_mapping );  // 处理lidar数据的线程
//    std::thread thr = std::thread( service_refresh_and_synchronize_triangle, 100 );                 // 移动数据，用于可视化线程
//    std::this_thread::sleep_for( std::chrono::milliseconds( 100 ) );
    Common_tools::Timer disp_tim;

    // 开启一个新线程，异步拷贝三角形
//    std::thread thr = std::thread( &GL_gui::service_refresh_and_synchronize_triangle, this );
    g_flag_pause = false;

    std::string gl_camera_file_name = Common_tools::get_home_folder().append( "/ImMeshing.gl_camera" );  // 没有啊
    g_gl_camera.load_camera( gl_camera_file_name );
    g_gl_camera.m_gl_cam.m_camera_z_far = 1500;
    g_gl_camera.m_gl_cam.m_camera_z_near = 0.1;

    // Rasterization configuration
    Cam_view m_depth_view_camera;
    m_depth_view_camera.m_display_w = 640;
    m_depth_view_camera.m_display_h = 480;
    m_depth_view_camera.m_camera_focus = 400;
    m_depth_view_camera.m_maximum_disp_depth = 150.0;
    m_depth_view_camera.m_draw_depth_pts_size = 2;
    m_depth_view_camera.m_draw_LiDAR_pts_size = m_depth_view_camera.m_draw_depth_pts_size;
    m_depth_view_camera.m_if_draw_depth_pts = true;
    vec_3 ext_rot_angle = vec_3( 0, 0, 0 );

    while ( !glfwWindowShouldClose( window ) )  // glfwWindowShouldClose 是 GLFW 库中的一个函数，用于检查窗口是否应该被关闭
    {
        // 渲染循环的开始，它设置了OpenGL的状态，初始化了ImGui，并准备了渲染环境。
        g_gl_camera.draw_frame_start();

        Eigen::Quaterniond q_last_avr;
        vec_3              t_last_avr;
        get_last_avr_pose( g_current_frame, q_last_avr, t_last_avr );
        if ( g_if_draw_depth )
        {
            Common_tools::Timer tim;
            tim.tic();
            m_depth_view_camera.m_camera_z_far = g_gl_camera.m_gl_cam.m_camera_z_far;
            m_depth_view_camera.m_camera_z_near = g_gl_camera.m_gl_cam.m_camera_z_near;

            if ( g_if_depth_bind_cam )
                m_depth_view_camera.set_camera_pose( q_last_avr.toRotationMatrix(), t_last_avr );
            else
                m_depth_view_camera.set_camera_pose( g_gl_camera.m_gl_cam.m_camera_rot, g_gl_camera.m_gl_cam.m_camera_pos );

            m_depth_view_camera.set_gl_matrix();

            draw_triangle( m_depth_view_camera );

            m_depth_view_camera.read_depth();
            glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );  // 清除颜色、内存缓冲区
            m_depth_view_camera.draw_depth_image();
            g_gl_camera.set_gl_camera_pose_matrix();
            if(m_depth_view_camera.m_if_draw_depth_pts)
            {
                g_draw_LiDAR_point = true;
            }
        }

        // 设置help界面
        // ImGui是一个用于构建即时模式图形界面的C++库
        if ( g_display_main_window )
        {
            // ImGui::Begin(...) 和 ImGui::End()：这两个函数调用标记了一个ImGui窗口的开始和结束。"ImMesh's Main_windows" 是窗口的标题，
            // g_display_main_window 是一个布尔值指针，用于控制窗口的显示和隐藏。
            ImGui::Begin( "Main_windows", &g_display_main_window );               // Create a window called "Hello, world!" and append into it.
            // ImGui::TextColored(...)：这个函数用于设置文本颜色，并显示文本。ImVec4 是一个包含RGBA值的向量，表示颜色。
//            ImGui::TextColored(ImVec4(1.0f, 1.0f, 1.0f, 0.75f), "ImMesh: An Immediate LiDAR Localization and Meshing Framework");
//            ImGui::TextColored(ImVec4(1.0f, 1.0f, 1.0f, 0.75f), "Github:  ");
//            // ImGui::SameLine()：这个函数用于将接下来的控件放置在同一行。
//            ImGui::SameLine();
//            ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "https://github.com/hku-mars/ImMesh");
//            ImGui::TextColored(ImVec4(1.0f, 1.0f, 1.0f, 0.75f), "Author:  ");
//            ImGui::SameLine();
//            ImGui::TextColored(ImVec4(0.5f, 0.5f, 1.0f, 1.0f), "Jiarong Lin & Chongjian Yuan");
            // ImGui::TreeNode(...) 和 ImGui::TreePop()：这两个函数调用用于创建一个可折叠的树形结构，TreeNode 开始一个新节点，TreePop 结束当前节点。
            if (ImGui::TreeNode("Help"))
            {
                // ImGui::Text(...)：显示普通文本。
                ImGui::Text( "[H]    | Display/Close main windows" );
                ImGui::Text( "[C]    | Show/Close camera pose windows" );
                ImGui::Text( "[T]    | Follow the camera" );
                ImGui::Text( "[D]    | Show/close the depth image" );
                ImGui::Text( "[L]    | Show/close the LiDAR points" );
                ImGui::Text( "[M]    | Show/close the mesh" );
                ImGui::Text( "[S]    | Save camera view" );
                ImGui::Text( "[Z]    | Load camera view" );
                ImGui::Text( "[+/-]  | Increase/Decrease the line width" );
                ImGui::Text( "[F1]   | Display help window" );
                ImGui::Text( "[Space]| To pause the program" );
                ImGui::Text( "[Esc]  | Exit the program" );
                ImGui::TreePop();
                // ImGui::Separator() 函数用于创建一个水平分隔线，它在GUI中起到视觉上分隔内容的作用
                ImGui::Separator();
            }
            // ImGui::SetNextItemOpen(...)：设置下一项打开或关闭。
            ImGui::SetNextItemOpen(true, 1);
            // ImGui::SetNextTreeNodeOpen();
            if (ImGui::TreeNode("Draw Online reconstructed mesh options:"))
            {
                // ImGui::RadioButton(...)：创建单选按钮，允许在多个选项中选择一个。
                ImGui::RadioButton("Draw mesh's Facet", &g_display_face, 1);
                ImGui::RadioButton("Draw mesh's Wireframe", &g_display_face, 0);
                // ImGui::Checkbox(...)：创建复选框，允许启用或禁用某个选项。
                // ImGui::Checkbox 函数的返回值是 bool 类型，表示复选框是否被勾选。
                // 如果返回值为 true，则表示复选框现在是勾选状态；如果为 false，则表示复选框没有被勾选。
                if(ImGui::Checkbox( "Draw mesh with color", &g_mesh_if_color ))
                {
                    g_force_refresh_triangle = true;
                }
                ImGui::TreePop();
                ImGui::Separator();
            }
            if (ImGui::TreeNode("LiDAR pointcloud reinforcement"))
            {
                ImGui::Checkbox( "Enable", &g_if_draw_depth );
                if(g_if_draw_depth)
                {
                    ImGui::Checkbox( "If depth in sensor frame", &g_if_depth_bind_cam );
                }
                // ImGui::SliderInt(...) 和 ImGui::SliderFloat(...)：创建整数和浮点数滑块，允许用户通过滑动选择一个值。
                ImGui::SliderInt( "Reinforced point size", &m_depth_view_camera.m_draw_depth_pts_size, 0, 10 );
                ImGui::SliderInt( "LiDAR point size", &m_depth_view_camera.m_draw_LiDAR_pts_size, 0, 10 );
                ImGui::TreePop();
                ImGui::Separator();
            }
            ImGui::Checkbox( "Move follow camera", &g_follow_cam );
            ImGui::Checkbox( "Mapping pause", &g_flag_pause );
            ImGui::Checkbox( "Draw LiDAR point", &g_draw_LiDAR_point );
            if(g_draw_LiDAR_point)
            {
                ImGui::SliderInt( "LiDAR point size", & m_depth_view_camera.m_draw_LiDAR_pts_size, 0, 10 );
            }
            ImGui::Checkbox( "Axis and Z_plane", &g_if_draw_z_plane );

            ImGui::SliderFloat( "Path width", &g_draw_path_size, 1.0, 10.0f );
            ImGui::SliderFloat( "Camera size", &g_display_camera_size, 0.01, 10.0, "%lf", ImGuiSliderFlags_Logarithmic );

            // ImGui::Button(...)：创建一个按钮，如果用户点击按钮，则返回true
            if ( ImGui::Button( "  Save Mesh to PLY file  " ) )
            {
//                int temp_flag = g_flag_pause;
//                g_flag_pause = true;
//                Common_tools::create_dir( data_path_file );
//                save_to_ply_file( std::string( data_path_file ).append( "/rec_mesh.ply" ), g_ply_smooth_factor, g_ply_smooth_k );
//                g_flag_pause = temp_flag;
            }

            if ( ImGui::Button( "Load Camera view" ) )
            {
                g_gl_camera.load_camera( gl_camera_file_name );
            }

            if ( ImGui::Button( "Save Camera view" ) )
            {
                cout << "Save view to " << gl_camera_file_name << endl;
                g_gl_camera.save_camera( gl_camera_file_name );
            }
            ImGui::Checkbox( "Show OpenGL camera paras", &g_display_camera_pose_window );
            // ImGui::GetIO().Framerate 获取当前的渲染帧率，它基于最近几帧的时间间隔来估算
            ImGui::Text( "Refresh rate %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate );
            if ( ImGui::Button( "      Exit Program      " ) ) // Buttons return true when clicked (most widgets return true when edited/activated)
            {
                // 设置窗口的关闭标志，如果设置为true，则窗口将关闭。
                glfwSetWindowShouldClose( window, 1 );
            }

            ImGui::End();
        }

        if(g_draw_LiDAR_point)
        {
            if ( m_depth_view_camera.m_if_draw_depth_pts )
            {
                if ( m_depth_view_camera.m_draw_LiDAR_pts_size > 0 )
                {
                    display_current_LiDAR_pts( g_current_frame, m_depth_view_camera.m_draw_LiDAR_pts_size, vec_4f( 1.0, 1.0, 1.0, 0.85 ) );
                }
                display_reinforced_LiDAR_pts( m_depth_view_camera.m_depth_pts_vec, m_depth_view_camera.m_draw_depth_pts_size, vec_3f( 1.0, 0.0, 1.0 ) );
            }
        }

        if ( g_follow_cam )
        {
            if ( g_current_frame > 1 )
            {
                g_gl_camera.tracking_camera( Eigen::Quaterniond::Identity(), t_last_avr );
            }
        }

        if ( g_display_help_win )
        {
//            // display help window
            print_help_window( &g_display_help_win );
        }

        if ( g_display_camera_pose_window )
        {
            g_gl_camera.draw_camera_window( g_display_camera_pose_window );
        }

        if ( g_if_draw_z_plane )
        {
            g_axis_shader.draw( g_gl_camera.m_gl_cam.m_glm_projection_mat,
                                Common_tools::eigen2glm( g_gl_camera.m_gl_cam.m_camera_pose_mat44_inverse ) );
            g_ground_plane_shader.draw( g_gl_camera.m_gl_cam.m_glm_projection_mat,
                                        Common_tools::eigen2glm( g_gl_camera.m_gl_cam.m_camera_pose_mat44_inverse ) );
        }

        if ( g_display_mesh )
        {
            draw_triangle( g_gl_camera.m_gl_cam );
        }

        // 显示相机位姿和轨迹
        if ( g_current_frame >= 0 )
        {
            draw_camera_pose( g_current_frame, g_draw_path_size, g_display_camera_size );
            draw_camera_trajectory( g_current_frame + 1, g_draw_path_size);
        }

        // For Key-board control
        if ( g_gl_camera.if_press_key( "H" ) )
        {
            g_display_main_window = !g_display_main_window;
        }
        if ( g_gl_camera.if_press_key( "C" ) )
        {
            g_display_camera_pose_window = !g_display_camera_pose_window;
        }
        if ( g_gl_camera.if_press_key( "F" ) )
        {
            g_display_face = !g_display_face;
        }
        if ( g_gl_camera.if_press_key( "Space" ) )
        {
            g_flag_pause = !g_flag_pause;
        }
        if ( g_gl_camera.if_press_key( "S" ) )
        {
            g_gl_camera.save_camera( gl_camera_file_name );
        }
        if ( g_gl_camera.if_press_key( "Z" ) )
        {
            g_gl_camera.load_camera( gl_camera_file_name );
        }
        if ( g_gl_camera.if_press_key( "D" ) )
        {
            g_if_draw_depth = !g_if_draw_depth;
        }
        if ( g_gl_camera.if_press_key( "M" ) )
        {
            g_display_mesh = !g_display_mesh;
        }
        if ( g_gl_camera.if_press_key( "T" ) )
        {
            g_follow_cam = !g_follow_cam;
            if ( g_current_frame > 1 )
            {
                g_gl_camera.set_last_tracking_camera_pos( q_last_avr, t_last_avr );
            }
        }
        if ( g_gl_camera.if_press_key( "Escape" ) )
        {
            glfwSetWindowShouldClose( window, 1 );
        }

        if ( g_gl_camera.if_press_key( "F1" ) )
        {
            g_display_help_win = !g_display_help_win;
        }

        g_gl_camera.set_gl_camera_pose_matrix();
        g_gl_camera.draw_frame_finish();
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow( window );
    glfwTerminate();

    /*
    ImGui_ImplOpenGL3_Shutdown()：
    这个函数调用是ImGui的OpenGL后端特定的函数，用于关闭和清理OpenGL相关的ImGui状态和资源。
     如果你的应用程序使用了ImGui的OpenGL3后端进行渲染，那么在应用程序结束前，你需要调用这个函数来确保释放所有与OpenGL相关的ImGui资源。

    ImGui_ImplGlfw_Shutdown()：
    这个函数调用是ImGui的GLFW后端特定的函数，用于关闭和清理GLFW相关的ImGui平台资源。
     如果你的应用程序使用了GLFW作为窗口创建和管理的库，并且集成了ImGui，那么在ImGui的OpenGL后端关闭后，你需要调用这个函数来清理GLFW平台上的ImGui资源。

    ImGui::DestroyContext()：
    这个函数用于销毁ImGui应用程序的全部上下文和状态。它会清理所有的ImGui资源，包括字体、样式、窗口状态等。
     一旦调用这个函数，所有的ImGui状态将被重置，如果之后还需要使用ImGui，需要重新初始化。

    glfwDestroyWindow(window)：
    这个GLFW函数用于销毁之前创建的窗口对象。window 参数是指向之前用GLFW创建的窗口的指针。
     一旦窗口被销毁，就不能再用它来渲染或接收事件了。

    glfwTerminate()：
    这个GLFW函数用于在应用程序结束时清理GLFW库，释放它所占用的所有资源。
     这包括关闭所有剩余的窗口、上下文和其他库资源。通常在程序退出前调用，确保没有资源泄露。
    */
}

// ANCHOR - synchronize_triangle_list_for_disp
// 另一个线程异步提取更新三角形列表
void GL_gui::synchronize_triangle_list_for_disp()
{
    int region_size = g_triangles_manager.m_triangle_set_vector.size();
    bool if_force_refresh = g_force_refresh_triangle;
    for ( int region_idx = 0; region_idx < region_size; region_idx++ )
    {
        Sync_triangle_set *                        sync_triangle_set_ptr = g_triangles_manager.m_triangle_set_vector[ region_idx ];
        std::shared_ptr< Region_triangles_shader > region_triangles_shader_ptr = nullptr;
        if ( g_map_region_triangles_shader.find( sync_triangle_set_ptr ) == g_map_region_triangles_shader.end() ) // 一个体素对应一个shader
        {
            // new a shader
            region_triangles_shader_ptr = std::make_shared< Region_triangles_shader >();
            g_region_triangle_shader_mutex.lock();
            g_map_region_triangles_shader.insert( std::make_pair( sync_triangle_set_ptr, region_triangles_shader_ptr ) );
            g_region_triangles_shader_vec.push_back( region_triangles_shader_ptr );
            g_region_triangle_shader_mutex.unlock();
        }
        else
        {
            region_triangles_shader_ptr = g_map_region_triangles_shader[ sync_triangle_set_ptr ];
        }

        if ( region_triangles_shader_ptr != nullptr )
        {
            if ( g_force_refresh_triangle && if_force_refresh == false )
            {
                if_force_refresh = true;
                region_idx = -1; // refresh from the start
            }
            if(if_force_refresh)
            {
                sync_triangle_set_ptr->m_if_required_synchronized = true;
            }
            region_triangles_shader_ptr->synchronized_from_region( sync_triangle_set_ptr, g_axis_min_max );
        }
    }
    if ( g_force_refresh_triangle )
    {
        g_force_refresh_triangle = false;
    }
}

void GL_gui::service_refresh_and_synchronize_triangle( )
{
    g_axis_min_max[ 0 ] = vec_3f( 1e8, 1e8, 1e8 );
    g_axis_min_max[ 1 ] = vec_3f( -1e8, -1e8, -1e8 );
//    g_axis_min_max[ 1 ] = vec_3f( 1e8, 1e8, 1e8 );
//    g_axis_min_max[ 0 ] = vec_3f( -1e8, -1e8, -1e8 );
//    int sleep_time = 100;
//    while ( 1 )
//    {
//        std::this_thread::sleep_for( std::chrono::milliseconds( sleep_time ) );  // 休眠100ms？
//        synchronize_triangle_list_for_disp();
//    }

    synchronize_triangle_list_for_disp();
}

Eigen::Vector3d GL_gui::get_tri_center(int A_idx, int B_idx, int C_idx, int idx)
{
    Eigen::Vector3d triangle_pos(0, 0, 0);
    if (idx == 0)
    {
        triangle_pos(0) = (non_plane_map->pts_list->points[A_idx].x +
                           non_plane_map->pts_list->points[B_idx].x +
                           non_plane_map->pts_list->points[C_idx].x) / 3.0;
        triangle_pos(1) = (non_plane_map->pts_list->points[A_idx].y +
                           non_plane_map->pts_list->points[B_idx].y +
                           non_plane_map->pts_list->points[C_idx].y) / 3.0;
        triangle_pos(2) = (non_plane_map->pts_list->points[A_idx].z +
                           non_plane_map->pts_list->points[B_idx].z +
                           non_plane_map->pts_list->points[C_idx].z) / 3.0;
    }
    else
    {
        auto temp_it = plane_idx2ptr.find(idx);
        if (temp_it == plane_idx2ptr.end())
        {
            cout << "\033[31m false palne idx! \033[0m" << endl;   // 输出红色字体
        }
        else
        {
            auto temp_mf = plane_idx2ptr[idx];
            triangle_pos(0) = (temp_mf->ptcl_grid->points[A_idx].x +
                                temp_mf->ptcl_grid->points[B_idx].x +
                                temp_mf->ptcl_grid->points[C_idx].x) / 3.0;
            triangle_pos(1) = (temp_mf->ptcl_grid->points[A_idx].y +
                                temp_mf->ptcl_grid->points[B_idx].y +
                                temp_mf->ptcl_grid->points[C_idx].y) / 3.0;
            triangle_pos(2) = (temp_mf->ptcl_grid->points[A_idx].z +
                                temp_mf->ptcl_grid->points[B_idx].z +
                                temp_mf->ptcl_grid->points[C_idx].z) / 3.0;
        }
    }
    return triangle_pos;
}

int GL_gui::if_index_flip(int A_idx, int B_idx, int C_idx)
{
    int less_than_num = 0;
    if (A_idx < B_idx) { less_than_num++; }
    if (B_idx < C_idx) { less_than_num++; }
    if (C_idx < A_idx) { less_than_num++; }
    if (less_than_num == 2)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

void GL_gui::triangle_manager_update(const std::vector<int>& update_plane_idx_list)
{
    /****** 非平面三角形更新 *****/
    // 确定要删除、添加的三角形索引列表
    std::vector<int> face_update_delete, face_update_add;
    for (int i = non_plane_map->last_faces_to_delete_num; i < non_plane_map->faces_to_delete.size(); i++)
    {
        face_update_delete.push_back(non_plane_map->faces_to_delete[i]);
    }
    non_plane_map->last_faces_to_delete_num = non_plane_map->faces_to_delete.size();

    for (int i = non_plane_map->last_faces_list_all_num; i < non_plane_map->faces_list_all.size(); i++)
    {
        auto temp_it = std::find(face_update_delete.begin(), face_update_delete.end(), i);
        if (temp_it == face_update_delete.end())
        {
            face_update_add.push_back(i);
        }
    }
    non_plane_map->last_faces_list_all_num = non_plane_map->faces_list_all.size();

    // 添加三角形到triangle_manager中
    for (auto cur_tri_ind : face_update_add)
    {
        // 新建Triangle对象
        int temp_tri_pts_id[ 3 ];
        temp_tri_pts_id[0] = non_plane_map->faces_list_all[cur_tri_ind][0];
        temp_tri_pts_id[1] = non_plane_map->faces_list_all[cur_tri_ind][1];
        temp_tri_pts_id[2] = non_plane_map->faces_list_all[cur_tri_ind][2];
        Triangle_ptr triangle_ptr;
        triangle_ptr = std::make_shared< Triangle >( temp_tri_pts_id[0], temp_tri_pts_id[1], temp_tri_pts_id[2] );
        triangle_ptr->m_index_flip = if_index_flip( temp_tri_pts_id[0], temp_tri_pts_id[1], temp_tri_pts_id[2] );
        triangle_ptr->rec_source = 0;

        // 插入到体素中
        Eigen::Vector3d triangle_pos = get_tri_center(temp_tri_pts_id[0], temp_tri_pts_id[1], temp_tri_pts_id[2], 0);
        int          hash_3d_x = std::round( triangle_pos( 0 ) / g_triangles_manager.m_region_size );     // 根据三角形重心点确定所属体素
        int          hash_3d_y = std::round( triangle_pos( 1 ) / g_triangles_manager.m_region_size );
        int          hash_3d_z = std::round( triangle_pos( 2 ) / g_triangles_manager.m_region_size );
        Sync_triangle_set* sync_triangle_set_ptr = g_triangles_manager.m_triangle_set_in_region.get_data( hash_3d_x, hash_3d_y, hash_3d_z );
        int temp_set_idx;
        if ( sync_triangle_set_ptr == nullptr )
        {
            sync_triangle_set_ptr = new Sync_triangle_set();
            sync_triangle_set_ptr->insert( triangle_ptr );
            sync_triangle_set_ptr->vec_idx = g_triangles_manager.m_triangle_set_vector.size();
            temp_set_idx = sync_triangle_set_ptr->vec_idx;

            g_triangles_manager.m_triangle_set_in_region.insert( hash_3d_x, hash_3d_y, hash_3d_z, *sync_triangle_set_ptr );
            g_triangles_manager.m_triangle_set_vector.push_back( g_triangles_manager.m_triangle_set_in_region.get_data( hash_3d_x, hash_3d_y, hash_3d_z ) );
        }
        else
        {
            sync_triangle_set_ptr->insert( triangle_ptr );
            temp_set_idx = sync_triangle_set_ptr->vec_idx;
        }

        // 建立Triangle_ptr Sync_triangle_set的索引
        m_triangle_hash_np[cur_tri_ind] = triangle_ptr;
        m_tri_set_np[cur_tri_ind] = temp_set_idx;
    }

    // 从triangle_manager中删除三角形
    for (auto cur_tri_ind : face_update_delete)
    {
        auto it_1 = m_tri_set_np.find(cur_tri_ind);
        auto it_2 = m_triangle_hash_np.find(cur_tri_ind);
        if (it_1 != m_tri_set_np.end() && it_2 != m_triangle_hash_np.end())
        {
            Sync_triangle_set* sync_triangle_set_ptr = g_triangles_manager.m_triangle_set_vector[m_tri_set_np[cur_tri_ind]];
            Triangle_ptr triangle_ptr = m_triangle_hash_np[cur_tri_ind];
            if ( sync_triangle_set_ptr != nullptr && triangle_ptr != nullptr)
            {
                sync_triangle_set_ptr->erase(triangle_ptr);
            }
        }
        m_triangle_hash_np[cur_tri_ind] = nullptr;
    }

    /****** 平面三角形更新  *****/
    // 平面索引去重复
    std::set<int> update_plane_idx_list_single;
    for (auto i : update_plane_idx_list)
    {
        update_plane_idx_list_single.insert(i);
    }

    for (auto plane_global_idx : update_plane_idx_list_single)
    {
        std::shared_ptr<MeshFragment> temp_mf = plane_idx2ptr[plane_global_idx];

        // 整个平面删除
        if (temp_mf->if_delete)
        {
            // 确定要删除的三角形索引列表
            std::vector<int> face_update_delete_all;
            for (int i = 0; i < temp_mf->last_faces_list_num; i++)
            {
                auto temp_it = std::find(temp_mf->faces_to_delete.begin(), temp_mf->faces_to_delete.end(), i);
                if (temp_it == temp_mf->faces_to_delete.end())
                {
                    face_update_delete_all.push_back(i);
                }
            }

            // 从triangle_manager中删除三角形
            for (auto cur_tri_ind : face_update_delete_all)
            {
                auto it_1 = m_tri_set_p[plane_global_idx].find(cur_tri_ind);
                auto it_2 = m_triangle_hash_p[plane_global_idx].find(cur_tri_ind);
                if (it_1 != m_tri_set_p[plane_global_idx].end() && it_2 != m_triangle_hash_p[plane_global_idx].end())
                {
                    Sync_triangle_set* sync_triangle_set_ptr = g_triangles_manager.m_triangle_set_vector[m_tri_set_p[plane_global_idx][cur_tri_ind]];
                    Triangle_ptr triangle_ptr = m_triangle_hash_p[plane_global_idx][cur_tri_ind];
                    if ( sync_triangle_set_ptr != nullptr && triangle_ptr != nullptr)
                    {
                        sync_triangle_set_ptr->erase(triangle_ptr);
                    }
                }
                m_triangle_hash_p[plane_global_idx][cur_tri_ind] = nullptr;
            }

            // 删除该平面所有三角形指针的索引
            m_triangle_hash_p[plane_global_idx].clear();
            m_tri_set_p[plane_global_idx].clear();

            // 删除对平面指针的索引
            auto temp_it = plane_idx2ptr.find(plane_global_idx);
            if (temp_it != plane_idx2ptr.end())
            {
                plane_idx2ptr.erase(temp_it);
            }

            continue;
        }

        // 确定要删除、添加的三角形索引列表
        std::vector<int> face_update_delete_p, face_update_add_p;
        int cur_faces_to_delete_size = temp_mf->faces_to_delete.size();
        for (int i = temp_mf->last_faces_to_delete_num; i < cur_faces_to_delete_size; i++)
        {
            face_update_delete_p.push_back(temp_mf->faces_to_delete[i]);
        }
        temp_mf->last_faces_to_delete_num = cur_faces_to_delete_size;

        int cur_faces_list_size = temp_mf->faces_list.size();
        for (int i = temp_mf->last_faces_list_num; i < cur_faces_list_size; i++)
        {
            auto temp_it = std::find(face_update_delete_p.begin(), face_update_delete_p.end(), i);
            if (temp_it == face_update_delete_p.end())
            {
                face_update_add_p.push_back(i);
            }
        }
        temp_mf->last_faces_list_num = cur_faces_list_size;

        // 添加三角形到triangle_manager中
        for (auto cur_tri_ind : face_update_add_p)
        {
            // 新建Triangle对象
            int temp_tri_pts_id[ 3 ];
            temp_tri_pts_id[0] = temp_mf->faces_list[cur_tri_ind][0];
            temp_tri_pts_id[1] = temp_mf->faces_list[cur_tri_ind][1];
            temp_tri_pts_id[2] = temp_mf->faces_list[cur_tri_ind][2];
            Triangle_ptr triangle_ptr;
            triangle_ptr = std::make_shared< Triangle >( temp_tri_pts_id[0], temp_tri_pts_id[1], temp_tri_pts_id[2] );
            triangle_ptr->m_index_flip = if_index_flip( temp_tri_pts_id[0], temp_tri_pts_id[1], temp_tri_pts_id[2] );
            triangle_ptr->rec_source = plane_global_idx;

            // 插入到体素中
            Eigen::Vector3d triangle_pos = get_tri_center(temp_tri_pts_id[0], temp_tri_pts_id[1], temp_tri_pts_id[2], plane_global_idx);
            int          hash_3d_x = std::round( triangle_pos( 0 ) / g_triangles_manager.m_region_size );     // 根据三角形重心点确定所属体素
            int          hash_3d_y = std::round( triangle_pos( 1 ) / g_triangles_manager.m_region_size );
            int          hash_3d_z = std::round( triangle_pos( 2 ) / g_triangles_manager.m_region_size );
            Sync_triangle_set* sync_triangle_set_ptr = g_triangles_manager.m_triangle_set_in_region.get_data( hash_3d_x, hash_3d_y, hash_3d_z );
            int temp_set_idx;
            if ( sync_triangle_set_ptr == nullptr )
            {
                sync_triangle_set_ptr = new Sync_triangle_set();
                sync_triangle_set_ptr->insert( triangle_ptr );
                sync_triangle_set_ptr->vec_idx = g_triangles_manager.m_triangle_set_vector.size();
                temp_set_idx = sync_triangle_set_ptr->vec_idx;

                g_triangles_manager.m_triangle_set_in_region.insert( hash_3d_x, hash_3d_y, hash_3d_z, *sync_triangle_set_ptr );
                g_triangles_manager.m_triangle_set_vector.push_back( g_triangles_manager.m_triangle_set_in_region.get_data( hash_3d_x, hash_3d_y, hash_3d_z ) );
            }
            else
            {
                sync_triangle_set_ptr->insert( triangle_ptr );
                temp_set_idx = sync_triangle_set_ptr->vec_idx;
            }

            // 建立Triangle_ptr的索引
            m_triangle_hash_p[plane_global_idx][cur_tri_ind] = triangle_ptr;
            m_tri_set_p[plane_global_idx][cur_tri_ind] = temp_set_idx;
        }

        // 从triangle_manager中删除三角形
        for (auto cur_tri_ind : face_update_delete_p)
        {
            auto it_1 = m_tri_set_p[plane_global_idx].find(cur_tri_ind);
            auto it_2 = m_triangle_hash_p[plane_global_idx].find(cur_tri_ind);
            if (it_1 != m_tri_set_p[plane_global_idx].end() && it_2 != m_triangle_hash_p[plane_global_idx].end())
            {
                Sync_triangle_set* sync_triangle_set_ptr = g_triangles_manager.m_triangle_set_vector[m_tri_set_p[plane_global_idx][cur_tri_ind]];
                Triangle_ptr triangle_ptr = m_triangle_hash_p[plane_global_idx][cur_tri_ind];
                if ( sync_triangle_set_ptr != nullptr && triangle_ptr != nullptr)
                {
                    sync_triangle_set_ptr->erase(triangle_ptr);
                }
            }
            m_triangle_hash_p[plane_global_idx][cur_tri_ind] = nullptr;
        }
    }
}