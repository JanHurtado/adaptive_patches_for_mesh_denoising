#ifndef METRIC_H
#define METRIC_H

#include "core/util.h"
#include "core/nanoflann.hpp"
#include "core/curvature.h"



using namespace nanoflann;


/** @defgroup groupKD-Tree KD-Tree group
 *  This group contains structures, definitions and functions for KD-Tree indexation
 *  @{
 */

/**
 * @brief The PointCloud struct: This struct is used by nanoflann library to represent point clouds indexed in a KD-Tree. See nanoflann documentation for more detail.
 */
struct PointCloud
{
    struct Point
    {
        AKNumber  x,y,z;
    };

    vector<Point>  pts;

    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline AKNumber kdtree_distance2(const AKNumber *p1, const size_t idx_p2,size_t size) const
    {
        AKNumber d0=p1[0]-pts[idx_p2].x;
        AKNumber d1=p1[1]-pts[idx_p2].y;
        AKNumber d2=p1[2]-pts[idx_p2].z;
        return d0*d0+d1*d1+d2*d2;
    }

    inline AKNumber kdtree_distance(const AKNumber *p1, const size_t idx_p2,size_t size) const
    {
        AKNumber d0=p1[0]-pts[idx_p2].x;
        AKNumber d1=p1[1]-pts[idx_p2].y;
        AKNumber d2=p1[2]-pts[idx_p2].z;
        return d0*d0+d1*d1+d2*d2;
    }

    inline AKNumber kdtree_get_pt(const size_t idx, int dim) const
    {
        if (dim==0) return pts[idx].x;
        else if (dim==1) return pts[idx].y;
        else return pts[idx].z;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const { return false; }

};

/**
 * @brief generateVertexPointCloud: Point cloud generation using vertex coordinates as points
 * @param mesh: input mesh
 * @param pointCloud: output point cloud
 */
void generate_vertex_point_cloud(AKMesh & mesh,PointCloud &pointCloud);

/**
 * @brief generateFaceCentroidPointCloud: Point cloud generation using face centroids as points
 * @param mesh: input mesh
 * @param pointCloud: output point cloud
 */
void generate_face_centroid_point_cloud(AKMesh & mesh,PointCloud &pointCloud);

/**
 * @brief AKKDTree: KD-Tree data structure definition with dimension = 3 and based on L2 distance.
 */
typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<AKNumber, PointCloud > ,
    PointCloud,
    3 /* dim */
    > AKKDTree;

/** @} */ // end of group1



AKNumber compute_mean_square_angle_error(AKMesh &desired_mesh , AKMesh & mesh);

AKNumber compute_area_error(AKNumber desiredArea, AKMesh & mesh);

AKNumber compute_volume_error(AKNumber desiredVol,AKMesh & mesh);

AKNumber compute_mean_vertex_distance_error(AKMesh & desired_mesh, AKKDTree & dmKDTreeVertices, AKKDTree &dmKDTreeFaces, AKMesh & mesh);

AKNumber compute_L2_vertex_based_error(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices, AKKDTree & dmKDTreeFaces, AKMesh & mesh);

AKNumber compute_mean_square_angle_error_2(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices, AKKDTree & dmKDTreeFaces, AKMesh & mesh);

AKNumber compute_mean_quadric_error(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices, AKMesh & mesh);

AKNumber compute_L2_normal_based_error(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices, AKKDTree & dmKDTreeFaces, AKMesh & mesh);

AKNumber compute_mean_tangential_error(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices, AKMesh & mesh);

AKNumber compute_mean_discrete_curvature_error(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices,vector<AKNumber> & dmCurvature, AKMesh & mesh, vector<AKNumber> & mCurvature);


inline AKNumber compute_mean_square_angle_error(AKMesh & desired_mesh,AKKDTree & dmKDTreeVertices, AKKDTree &dmKDTreeFaces,vector<AKNumber> & dmCurvature,AKMesh & mesh, vector<AKNumber> & mCurvature)
{
    return compute_mean_square_angle_error(desired_mesh , mesh);
}

inline AKNumber compute_area_error(AKMesh & desired_mesh,AKKDTree & dmKDTreeVertices, AKKDTree &dmKDTreeFaces,vector<AKNumber> & dmCurvature,AKMesh & mesh, vector<AKNumber> & mCurvature)
{
    return compute_area_error(compute_area(desired_mesh),mesh);
}

inline AKNumber compute_volume_error(AKMesh & desired_mesh,AKKDTree & dmKDTreeVertices, AKKDTree &dmKDTreeFaces,vector<AKNumber> & dmCurvature,AKMesh & mesh, vector<AKNumber> & mCurvature)
{
    return compute_volume_error(compute_volume(desired_mesh),mesh);
}

inline AKNumber compute_mean_vertex_distance_error(AKMesh & desired_mesh,AKKDTree & dmKDTreeVertices, AKKDTree &dmKDTreeFaces,vector<AKNumber> & dmCurvature,AKMesh & mesh, vector<AKNumber> & mCurvature)
{
    return compute_mean_vertex_distance_error(desired_mesh, dmKDTreeVertices,dmKDTreeFaces, mesh);
}

inline AKNumber compute_L2_vertex_based_error(AKMesh & desired_mesh,AKKDTree & dmKDTreeVertices, AKKDTree &dmKDTreeFaces,vector<AKNumber> & dmCurvature,AKMesh & mesh, vector<AKNumber> & mCurvature)
{
    return compute_L2_vertex_based_error(desired_mesh,dmKDTreeVertices,dmKDTreeFaces,mesh);
}

inline AKNumber compute_mean_square_angle_error_2(AKMesh & desired_mesh,AKKDTree & dmKDTreeVertices, AKKDTree &dmKDTreeFaces,vector<AKNumber> & dmCurvature,AKMesh & mesh, vector<AKNumber> & mCurvature)
{
    return compute_mean_square_angle_error_2(desired_mesh,dmKDTreeVertices,dmKDTreeFaces,mesh);
}

inline AKNumber compute_mean_quadric_error(AKMesh & desired_mesh,AKKDTree & dmKDTreeVertices, AKKDTree &dmKDTreeFaces,vector<AKNumber> & dmCurvature,AKMesh & mesh, vector<AKNumber> & mCurvature)
{
    return compute_mean_quadric_error(desired_mesh,dmKDTreeVertices,mesh);
}

inline AKNumber compute_L2_normal_based_error(AKMesh & desired_mesh,AKKDTree & dmKDTreeVertices, AKKDTree &dmKDTreeFaces,vector<AKNumber> & dmCurvature,AKMesh & mesh, vector<AKNumber> & mCurvature)
{
    return compute_L2_normal_based_error(desired_mesh,dmKDTreeVertices,dmKDTreeFaces,mesh);
}

inline AKNumber compute_mean_tangential_error(AKMesh & desired_mesh,AKKDTree & dmKDTreeVertices, AKKDTree &dmKDTreeFaces,vector<AKNumber> & dmCurvature,AKMesh & mesh, vector<AKNumber> & mCurvature)
{
    return compute_mean_tangential_error(desired_mesh,dmKDTreeVertices,mesh);
}

inline AKNumber compute_mean_discrete_curvature_error(AKMesh & desired_mesh,AKKDTree & dmKDTreeVertices, AKKDTree &dmKDTreeFaces,vector<AKNumber> & dmCurvature,AKMesh & mesh, vector<AKNumber> & mCurvature)
{
    return compute_mean_discrete_curvature_error(desired_mesh,dmKDTreeVertices,dmCurvature,mesh,mCurvature);
}

#endif // METRIC_H
