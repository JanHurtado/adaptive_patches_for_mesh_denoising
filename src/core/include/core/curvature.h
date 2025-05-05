#ifndef CURVATURE_H
#define CURVATURE_H

#include "core/util.h"

void compute_vertex_gaussian_curvature(AKMesh &mesh, vector<AKNumber> & curvature);

void compute_vertex_mean_curvature(AKMesh & mesh,vector<AKNumber> & curvature);

void compute_vertex_principal_curvatures_sum(AKMesh & mesh,vector<AKNumber> & mean_curvature,vector<AKNumber> & gaussian_curvature,vector<AKNumber> & principal_curvature);

void compute_all_curvatures(AKMesh & mesh,vector<AKNumber> & mean_curvature,vector<AKNumber> & gaussian_curvature,vector<AKNumber> & principal_curvature);


#endif // CURVATURE_H
