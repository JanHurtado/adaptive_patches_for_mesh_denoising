#ifndef UTIL_H
#define UTIL_H

#include "core/mesh.h"
#include "core/neighborhood.h"
#include "core/solver.h"
#include <vector>
#include <queue>
#include <ilcplex/cplex.h>


using namespace std;

#define PI 3.14159265359

/** @addtogroup mesh_processing
  * @brief Mesh processing auxiliar functions.
  *
  * @{
  */

/**
 * @brief computeArea - area computation of a mesh (sum of triangle areas)
 * @param mesh input mesh
 * @return total area
 */
AKNumber compute_area (AKMesh & mesh);

/**
 * @brief computeVolume - volume computation of a mesh
 * @param mesh input mesh
 * @return total area
 */
AKNumber compute_volume ( AKMesh & mesh );

/**
 * @brief computeAverageEdgeLength - average edge length computation of a mesh
 * @param mesh input mesh
 * @return average edge length
 */
AKNumber compute_average_edge_length(AKMesh &mesh);

/**
 * @brief computeAllFaceAreas - area computation for all faces
 * @param mesh input mesh
 * @param areas output vector containing areas for all faces
 */
void compute_all_face_areas(AKMesh &mesh, vector<AKNumber> &areas);

/**
 * @brief computeAllFaceCentroids - centroid computation for all faces
 * @param mesh input mesh
 * @param centroids output vector containing centroids for all faces
 */
void compute_all_face_centroids(AKMesh &mesh, vector<AKMesh::Point> &centroids);

/**
 * @brief computeAllFaceNormals - normal computation for all faces
 * @param mesh input mesh
 * @param normals vector containing normals for all faces
 */
void compute_all_face_normals(AKMesh &mesh, vector<AKMesh::Normal> &normals);

/**
 * @brief computeFaceVertexAngle - angle computation for a given vertex included in a face
 * @param mesh input mesh
 * @param fh corresponding face handle
 * @param vh corresponding vertex handle
 * @param angle angle
 */
void compute_face_vertex_angle(AKMesh &mesh,AKMesh::FaceHandle fh,AKMesh::VertexHandle vh, AKNumber & angle);

/**
 * @brief computeVertexArea - vertex area computation (barycentric area)
 * @param mesh input mesh
 * @param vh input vertex (vertex handle)
 * @param areas all precomputed triangle areas
 * @return vertex area
 */
AKNumber compute_vertex_area(AKMesh & mesh,AKMesh::VertexHandle vh,vector<AKNumber> & areas);

/**
 * @brief computeAllVertexAreas - vertex area computation for all vertices (barycentric areas)
 * @param mesh input mesh
 * @param areas output vector containing all vertex areas
 */
void compute_all_vertex_areas(AKMesh & mesh, vector<AKNumber> & areas);

/**
 * @brief computeAllPoints - compute all vertex coordinates
 * @param mesh input mesh
 * @param points all vertex coordinates
 */
void compute_all_points(AKMesh &mesh, vector<AKMesh::Point> & points);

/**
 * @brief GaussianWeight - gaussian function computation used for bilateral filtering
 * @param distance input distance
 * @param sigma input sigma
 * @return function value for inputs (distance and sigma)
 */
inline AKNumber compute_gaussian_weight(AKNumber distance, AKNumber sigma)
{
    return static_cast<AKNumber>(exp( -0.5 * distance * distance / (sigma * sigma)));
}

/**
 * @brief NormalDistance - computation of normal distance: \f$ |n1-n2| \f$
 * @param n1 - input normal 1
 * @param n2 - input normal 2
 * @return normal distance
 */
inline AKNumber normal_distance(const AKMesh::Normal &n1, const AKMesh::Normal &n2)
{
    return (n1 - n2).length();
}

/** @} */

void generate_data_vectors(AKMesh & _mesh, std::vector<AKNumber> & _points, std::vector<unsigned> & _faces);

void normalize_mesh(AKMesh & _mesh);

void normalize_mesh(AKMesh & _mesh, AKNumber & scale_factor);

void retrieve_mesh(AKMesh & _mesh, AKNumber scale_factor);

void writeData(AKMatrix & H, AKMatrix & A, AKNumber area_const);

void computeAllFaceNormals2(AKMesh &mesh, vector<AKMesh::Normal> &normals);

#endif // UTIL_H
