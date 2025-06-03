//
// Created by neo on 2024/4/7.
//

#ifndef SRC_HASH_IDX_HPP
#define SRC_HASH_IDX_HPP

#include "unordered_map"
#include "iostream"

template <typename data_type = float, typename T = void *>
struct Hash_map_3d_base
{
    using hash_3d_T = std::unordered_map<data_type, std::unordered_map<data_type, std::unordered_map<data_type, T>>>;
    hash_3d_T m_map_3d_hash_map;
    void insert(const data_type &x, const data_type &y, const data_type &z, const T &target)
    {
        m_map_3d_hash_map[x][y][z] = target;
    }

    int if_exist(const data_type &x, const data_type &y, const data_type &z)
    {
        if(m_map_3d_hash_map.find(x) == m_map_3d_hash_map.end()  )
        {
            return 0;
        }
        else if(m_map_3d_hash_map[x].find(y) ==  m_map_3d_hash_map[x].end() )
        {
            return 0;
        }
        else if( m_map_3d_hash_map[x][y].find(z) == m_map_3d_hash_map[x][y].end() )
        {
            return 0;
        }
        return 1;
    }

    void clear()
    {
        m_map_3d_hash_map.clear();
    }

    int total_size()
    {
        int count =0 ;
        for(auto it : m_map_3d_hash_map)
        {
            for(auto it_it: it.second)
            {
                for( auto it_it_it: it_it.second )
                {
                    count++;
                }
            }
        }
        return count;
    }
};

//
//template <typename data_type = float, typename T = void *>
//struct Hash_map_2d
//{
//    using hash_2d_T = std::unordered_map<data_type, std::unordered_map<data_type, T> >;
//    // using hash_2d_it = typename std::unordered_map<data_type, std::unordered_map<data_type, T> >::iterator ;
//    // using hash_2d_it_it = typename std::unordered_map<data_type, T>::iterator ;
//
//    hash_2d_T m_map_2d_hash_map;
//    void insert(const data_type &x, const data_type &y,  const T &target)
//    {
//        m_map_2d_hash_map[x][y] = target;
//    }
//
//    int if_exist(const data_type &x, const data_type &y )
//    {
//        if(m_map_2d_hash_map.find(x) == m_map_2d_hash_map.end()  )
//        {
//            return 0;
//        }
//        else if(m_map_2d_hash_map[x].find(y) ==  m_map_2d_hash_map[x].end() )
//        {
//            return 0;
//        }
//
//        return 1;
//    }
//
//    void clear()
//    {
//        m_map_2d_hash_map.clear();
//    }
//
//    int total_size()
//    {
//        int count =0 ;
//        //for(hash_2d_it it =  m_map_2d_hash_map.begin(); it != m_map_2d_hash_map.end(); it++)
//        for(auto it : m_map_2d_hash_map)
//        {
//            for(auto it_it: it.second)
//            {
//                count++;
//            }
//        }
//        return count;
//    }
//};







#endif //SRC_HASH_IDX_HPP
