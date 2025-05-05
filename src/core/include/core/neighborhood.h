#ifndef NEIGHBORHOOD_H
#define NEIGHBORHOOD_H

#include <core/mesh.h>
#include <vector>
#include <queue>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

/** @addtogroup mesh_processing
  * @brief Mesh processing neighborhood functions.
  *
  * @{
  */

/////////////////////////////////////////////////
/// Vertex neighborhood
/////////////////////////////////////////////////

/**
 * @brief computeVertexNeighbors: compute all neighboring vertices of a given vertex regarding a depth k (k-ring)
 * @param mesh: input mesh
 * @param vh: evaluated vertex
 * @param k: depth
 * @param vertex_neighbors: vector containing neighbors (output)
 */
void compute_vertex_neighbors(AKMesh &mesh, AKMesh::VertexHandle vh, int k, vector<AKMesh::VertexHandle> &vertex_neighbors);

/**
 * @brief computeAdaptiveVertexNeighbors: compute all neighboring vertices of a given vertex regarding a radius
 * @param mesh: input mesh
 * @param vh: evaluated vertex
 * @param radius: radius
 * @param vertex_neighbors: vector containing neighbors (output)
 */
void compute_adaptive_vertex_neighbors(AKMesh &mesh, AKMesh::VertexHandle vh, AKNumber radius, vector<AKMesh::VertexHandle> &vertex_neighbors);

/**
 * @brief computeAllVertexNeighbors: compute all neighboring vertices of all vertices regarding a depth k (k-ring)
 * @param mesh: input mesh
 * @param k: depth
 * @param all_vertex_neighbors: vector containing all neighbor vectors (output)
 */
void compute_all_vertex_neighbors(AKMesh &mesh, int k,vector<vector<AKMesh::VertexHandle> > &all_vertex_neighbors);

/**
 * @brief computeAllAdaptiveVertexNeighbors: compute all neighboring vertices of all vertices regarding a radius
 * @param mesh: input mesh
 * @param radius: radius
 * @param all_vertex_neighbors: vector containing all neighbor vectors (output)
 */
void compute_all_adaptive_vertex_neighbors(AKMesh &mesh, AKNumber radius, vector<vector<AKMesh::VertexHandle> > &all_vertex_neighbors);


/////////////////////////////////////////////////
/// Face neighborhood
/////////////////////////////////////////////////

/**
 * @brief computeFaceNeighbors_EdgeBased: compute neighboring faces based on face edges
 * @param mesh: input mesh
 * @param fh: evaluated face
 * @param face_neighbors: vector containing neighbors (output)
 */
void compute_edge_based_face_neighbors(AKMesh &mesh, AKMesh::FaceHandle fh, vector<AKMesh::FaceHandle> &face_neighbors);

/**
 * @brief computeFaceNeighbors_VertexBased: compute neighboring faces based on face vertices
 * @param mesh: input mesh
 * @param fh: evaluated face
 * @param face_neighbors: vector containing neighbors (output)
 */
void compute_vertex_based_face_neighbors(AKMesh &mesh, AKMesh::FaceHandle fh, vector<AKMesh::FaceHandle> &face_neighbors);

/**
 * @brief computeFaceNeighbors_RadiusBased: compute neighboring faces regarding a radius
 * @param mesh: input mesh
 * @param fh: evaluated face
 * @param radius: radius
 * @param face_neighbors: vector containing neighbors (output)
 */
void compute_radius_based_face_neighbors(AKMesh &mesh, AKMesh::FaceHandle fh, AKNumber radius, vector<AKMesh::FaceHandle> &face_neighbors);

/**
 * @brief computeAllFaceNeighbors_EdgeBased: compute all neighboring faces of all faces (edge based)
 * @param mesh: input mesh
 * @param include_tarcompute_face: include fh as neighbor
 * @param all_face_neighbors: vector containing neighbors (output)
 */
void compute_all_edge_based_face_neighbors(AKMesh &mesh, bool include_tarcompute_face, vector<vector<AKMesh::FaceHandle> > &all_face_neighbors);

/**
 * @brief computeAllFaceNeighbors_VertexBased: compute all neighboring faces of all faces (vertex based)
 * @param mesh: input mesh
 * @param include_tarcompute_face: include fh as neighbor
 * @param all_face_neighbors: vector containing neighbors (output)
 */
void compute_all_vertex_based_face_neighbors(AKMesh &mesh, bool include_tarcompute_face, vector<vector<AKMesh::FaceHandle> > &all_face_neighbors);

/** @} */

#endif // NEIGHBORHOOD_H
