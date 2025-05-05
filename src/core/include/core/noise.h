#ifndef NOISE_H
#define NOISE_H

#include "core/mesh.h"

AKNumber generate_random_gaussian(AKNumber mean, AKNumber StandardDeviation);

void generate_random_gaussian_numbers(AKNumber mean, AKNumber StandardDeviation, int number, vector<AKNumber> &RandomNumbers);

AKMesh add_noise(AKMesh & _mesh);

AKMesh add_noise(AKMesh & _mesh, AKNumber intensity);

AKMesh::Normal generate_random_direction();

void generate_random_directions(int number, vector<AKMesh::Normal> &RandomDirections);

AKMesh add_noise_random_direction(AKMesh & _mesh);

#endif // NOISE_H
