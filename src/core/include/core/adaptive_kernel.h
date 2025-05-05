#ifndef ADAPTIVE_KERNEL_H
#define ADAPTIVE_KERNEL_H


#include "util.h"
#include "solver.h"

void compute_distance_matrix(AKMesh & _mesh, AKNumber _max_dist, vector<vector<pair<size_t, AKNumber>>> & _distance_matrix_sparse, int _max_num_vars);

void sort_distance_matrix(vector<vector<pair<size_t, AKNumber>>> & _distance_matrix_sparse);

AKNumber get_distance_matrix_value_at(vector<vector<pair<size_t, AKNumber>>> & _distance_matrix_sparse, int row, int column, AKNumber default_value);

void compute_all_adaptive_kernels(AKMesh & _mesh, vector<AKNumber> & areas, vector<AKMesh::Point> & centroids, vector<AKMesh::Normal> & normals, vector<vector<pair<size_t, AKNumber>>> & _distance_matrix_sparse, AKNumber alpha, AKNumber beta, AKNumber gamma, AKNumber delta, AKNumber k, vector<vector<pair<size_t, AKNumber>>> & _patches_sparse);

void compute_adaptive_kernel(AKMesh & _mesh, int _face_index, vector<AKNumber> & areas, vector<AKMesh::Point> & centroids, vector<AKMesh::Normal> & normals, vector<pair<size_t, AKNumber>> & _distance_function, AKNumber alpha, AKNumber beta, AKNumber gamma, AKNumber delta, AKNumber k, vector<pair<size_t, AKNumber>> & _patch_function);

#endif // ADAPTIVE_KERNEL_H
