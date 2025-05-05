#ifndef DENOISING_H
#define DENOISING_H

#include "core/neighborhood.h"
#include "core/util.h"
#include "core/adaptive_kernel.h"


/** @addtogroup mesh_processing
  * @brief Mesh denoising functions
  *
  * @{
  */

/**
 * @brief The FaceNeighborType enum - determines the face neighborhood features
 */
enum FaceNeighborType{
    kVertexBased, /**< vertex based face neighborhood */
    kEdgeBased, /**< edge based face neighborhood */
    kRadiusBased /**< radius based face neighborhood */
};

/**
 * @brief updateVertexPositions - computes new vertex positions adapted to an input normal field.
 * @param mesh input mesh.
 * @param filtered_normals input normal field.
 * @param iteration_number number of iterations of the optimization algorithm (gradient descent approach).
 * @param fixed_boundary determines if a boundary vertex will be updated
 */
void update_vertex_positions(AKMesh &mesh, vector<AKMesh::Normal> &filtered_normals, int iteration_number, bool fixed_boundary);

/**
 * @brief computeSigmaC - sigma_c computation based on the average of face centroid distances (only between adjacent faces), and multiplied by a scalar
 * @param mesh input mesh
 * @param face_centroids precomputed face centroids
 * @param sigma_c_scalar scalar
 * @return value of sigma_c
 */
AKNumber compute_sigma_spatial(AKMesh &mesh, vector<AKMesh::Point> &face_centroids, AKNumber sigma_c_scalar);

/**
 * @brief computeRadius - computation of radius for radius-based face neighborhood, based on the average of face centroid distances (only between adjacent faces) and multiplied by a scalar
 * @param mesh input mesh
 * @param scalar
 * @return radius value
 */
AKNumber compute_radius(AKMesh &mesh, AKNumber scalar);

/////////////////////////////////////////////////
/// Uniform Laplacian smoothing
/////////////////////////////////////////////////

/**
 * @brief uniformLaplacian - denoise the input mesh using Uniform Laplacian Smoothing Algorithm
 * @param _mesh input mesh
 * @param iteration_number number of iterations
 * @param scale displacement influence for each iteration
 * @return denoised mesh
 */
AKMesh uniform_laplacian_smoothing(AKMesh & _mesh, int iteration_number, AKNumber scale);

/**
 * @brief uniformLaplacian - denoise a subset of vertices of an input mesh using Uniform Laplacian Smoothing Algorithm
 * @param _mesh input mesh
 * @param iteration_number number of iterations
 * @param scale displacement influence for each iteration
 * @param vertex_ids subset of vertex IDs to be denoised
 * @return denoised mesh (taking into account only the subset of IDs)
 */
AKMesh uniform_laplacian_smoothing(AKMesh & _mesh, int iteration_number, AKNumber scale, vector<size_t> & vertex_ids);

/////////////////////////////////////////////////
/// HC Laplacian smoothing (Vollmer et al.)
/////////////////////////////////////////////////

/**
 * @brief HCLaplacian - denoise the input mesh using HC Laplacian Smoothing Algorithm
 * @param _mesh input mesh
 * @param iteration_number number of iterations
 * @param alpha alpha defined in HC Laplacian Smoothing paper
 * @param beta beta defined in HC Laplacian Smoothing paper
 * @return denoised mesh
 */
AKMesh laplacian_smoothing_HC(AKMesh & _mesh,int iteration_number,AKNumber alpha,AKNumber beta);

/**
 * @brief HCLaplacian - denoise a subset of vertices of an input mesh using HC Laplacian Smoothing Algorithm
 * @param _mesh input mesh
 * @param iteration_number number of iterations
 * @param alpha alpha defined in HC Laplacian Smoothing paper
 * @param beta beta defined in HC Laplacian Smoothing paper
 * @param vertex_ids subset of vertex IDs to be denoised
 * @return denoised mesh (taking into account only the subset of IDs)
 */
AKMesh laplacian_smoothing_HC(AKMesh & _mesh, int iteration_number, AKNumber alpha, AKNumber beta, vector<size_t> & vertex_ids);

/////////////////////////////////////////////////
/// Bilateral normal filtering for mesh denoising (Zheng et al.)
/////////////////////////////////////////////////

/**
 * @brief updateFilteredNormals - bilateral filtering of the normal field (face based) of a mesh
 * @param mesh input mesh
 * @param normal_iteration_number number of iterations
 * @param sigma_c_scalar bilateral filtering influence regarding a \f$ sigma_c \f$ based on average edge length (spatial distance influence)
 * @param sigma_s bilateral filtering parameter \f$ sigma_s \f$
 * @param filtered_normals output vector of filtered face normals
 */
void bilateral_filtering(AKMesh &mesh, int normal_iteration_number, AKNumber sigma_c_scalar, AKNumber sigma_s, vector<AKMesh::Normal> &filtered_normals);

/**
 * @brief bilateralNormal - Denoise the input mesh using Bilateral Normal Filtering Algorithm
 * @param _mesh input mesh
 * @param normal_iteration_number number of iterations for normal field filtering
 * @param vertex_iteration_number number of iterations for vertex updating
 * @param sigma_c_scalar bilateral normal field filtering influence regarding a \f$ sigma_c \f$ based on average edge length (spatial distance influence)
 * @param sigma_s bilateral normal field filtering parameter \f$ sigma_s \f$
 * @return denoised mesh
 */
AKMesh bilateral_normal_filtering(AKMesh &_mesh, int normal_iteration_number, int vertex_iteration_number, AKNumber sigma_c_scalar,  AKNumber sigma_s );

/////////////////////////////////////////////////
/// Guided mesh normal filtering (Zhang et al.)
/////////////////////////////////////////////////

/**
 * @brief computeAllFaceNeighborsGMNF - neighborhood computation for all mesh faces (for bilateral filtering)
 * @param mesh input mesh
 * @param face_neighbor_type type of neighborhood to be considered
 * @param radius radius for radius based face neighborhood
 * @param include_central_face if the same face is considered as a neighbor
 * @param all_face_neighbors output vector containing neighboring faces IDs for each face of the mesh
 */
void compute_all_face_neighbors_guided_bilateral_normal_filtering(AKMesh &mesh, FaceNeighborType face_neighbor_type, AKNumber radius, bool include_central_face,
    vector<vector<AKMesh::FaceHandle> > &all_face_neighbors);

/**
 * @brief computeAllGuidedNeighborsGMNF - neighborhood computation for guidance signal
 * @param mesh input mesh
 * @param all_guided_neighbors output vector containing neighboring faces IDs for each face of the mesh
 */
void compute_all_guidance_signal_neighbors(AKMesh &mesh, vector<vector<AKMesh::FaceHandle> > &all_guided_neighbors);

/**
 * @brief computeFaceNeighborsInnerEdges - compute all inner edges of a face neighborhood
 * @param mesh input mesh
 * @param face_neighbors face neighborhood
 * @param inner_edges output vector containing all inner edges
 */
void compute_face_neighbors_inner_edges(AKMesh &mesh, vector<AKMesh::FaceHandle> &face_neighbors, vector<AKMesh::EdgeHandle> &inner_edges);

/**
 * @brief computeConsistenciesAndMeanNormals - consistency and mean normal computation for all patches
 * @param mesh input mesh
 * @param all_guided_neighbors all neighborhoods defining the patches
 * @param face_areas all face areas
 * @param normals all face normals
 * @param consistencies_and_mean_normals output vector containing the consistency and mean normal of all patches
 */
void compute_consistencies_and_mean_normals(AKMesh &mesh, vector<vector<AKMesh::FaceHandle> > &all_guided_neighbors,
    vector<AKNumber> &face_areas, vector<AKMesh::Normal> &normals,
    vector<pair<AKNumber, AKMesh::Normal> > &consistencies_and_mean_normals);

/**
 * @brief computeGuidedNormals - guided normal computation for each face of the mesh
 * @param mesh input mesh
 * @param all_guided_neighbors all neighborhoods defining the patches
 * @param face_areas all face areas
 * @param normals all face normals
 * @param consistencies_and_mean_normals consistencies and mean normals of all patches
 * @param guided_normals output vector containing guided normals for all faces
 */
void compute_guided_normals(AKMesh &mesh, vector<vector<AKMesh::FaceHandle> > &all_guided_neighbors,
    vector<AKNumber> &face_areas, vector<AKMesh::Normal> &normals,
    vector<pair<AKNumber, AKMesh::Normal> > consistencies_and_mean_normals,
	vector<AKMesh::Normal> &guided_normals);

/**
 * @brief updateFilteredNormalsGuided - guided bilateral filtering of a normal field (face based) of a mesh
 * @param mesh input mesh
 * @param filtered_normals guided normals for each face
 * @param radius_scalar average radius ratio for face neighborhood computation
 * @param sigma_c_scalar bilateral filtering influence regarding a \f$ sigma_c \f$ based on average edge length (spatial distance influence)
 * @param normal_iteration_number number of iterations for normal field filtering
 * @param sigma_s bilateral filtering parameter \f$ sigma_s \f$
 * @param vertex_iteration_number number of iterations for vertex updating
 */
void guided_bilateral_normal_filtering(AKMesh &mesh, vector<AKMesh::Normal> &filtered_normals, AKNumber radius_scalar,
	AKNumber sigma_c_scalar, int normal_iteration_number, AKNumber sigma_s, int vertex_iteration_number);

/**
 * @brief guided - Denoise the input mesh using Guided Mesh Normal Filtering Algorithm
 * @param _mesh input mesh
 * @param normal_iteration_number number of iterations for normal field filtering
 * @param vertex_iteration_number number of iterations for vertex updating
 * @param sigma_c_scalar bilateral normal field filtering influence regarding a \f$ sigma_c \f$ based on average edge length (spatial distance influence)
 * @param sigma_s bilateral normal field filtering parameter \f$ sigma_s \f$
 * @param radius_scalar average radius ratio for face neighborhood computation
 * @return denoised mesh
 */
AKMesh guided_denoising(AKMesh _mesh, int normal_iteration_number, int vertex_iteration_number, AKNumber sigma_c_scalar, AKNumber sigma_s, AKNumber radius_scalar);

/** @} */
vector<AKMesh::Normal> adaptive_kernels_normal_filtering(AKMesh & mesh, int normal_iteration_number, int bilateral_iteration_number, int vertex_iteration_number, AKNumber sigma_c_scalar, AKNumber sigma_s, AKNumber alpha, AKNumber beta, AKNumber gamma, AKNumber delta, AKNumber k);

AKMesh adaptive_kernels_denoising(AKMesh _mesh, int normal_iteration_number, int bilateral_iteration_number, int vertex_iteration_number, int external_iteration_number, AKNumber sigma_c_scalar, AKNumber sigma_s, AKNumber alpha, AKNumber beta, AKNumber gamma, AKNumber delta, AKNumber k);




//void computeGeodesics();




#endif // DENOISING_H
