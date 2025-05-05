#include "core/noise.h"

const AKNumber default_noise_level = 0.2;

AKNumber generate_random_gaussian(AKNumber mean, AKNumber StandardDeviation)
{
    static AKNumber v1, v2, s;
    static int phase = 0;
    AKNumber x;

    if(phase == 0)
    {
        do
        {
            v1 = -1 + 2 * (AKNumber)rand() / (AKNumber) RAND_MAX;
            v2 = -1 + 2 * (AKNumber)rand() / (AKNumber) RAND_MAX;
            s = v1 * v1 + v2 * v2;
        }while (s >= 1 || s == 0);

        x = v1 * sqrt(-2 * log(s) / s);
    }
    else
        x = v2 * sqrt(-2 * log(s) / s);

    phase = 1 - phase;

    return x * StandardDeviation + mean;
}

void generate_random_gaussian_numbers(AKNumber mean, AKNumber StandardDeviation, int number, vector<AKNumber> &RandomNumbers)
{
    RandomNumbers.resize(number, 0.0);

    srand((unsigned int)time(NULL));
    for(int i = 0; i < number; i++){
        RandomNumbers[i] = generate_random_gaussian(mean, StandardDeviation);
    }
}

AKMesh add_noise(AKMesh & _mesh)
{
    AKMesh mesh = _mesh;
    // compute average length of mesh
    AKNumber average_length = 0.0;
    for(AKMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
        average_length += mesh.calc_edge_length(*e_it);
    AKNumber edge_numbers = (AKNumber)mesh.n_edges();
    average_length /= edge_numbers;
    AKNumber noise_level = default_noise_level;
    AKNumber standard_deviation = average_length * noise_level;

    mesh.request_face_normals();
    mesh.request_vertex_normals();
    mesh.update_normals();

    vector<AKNumber> GaussianNumbers;

    generate_random_gaussian_numbers(0, standard_deviation, (int)mesh.n_vertices(), GaussianNumbers);

    for(AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++){

        AKMesh::Point p = mesh.point(*v_it) + mesh.normal(*v_it) * GaussianNumbers[v_it->idx()];
        mesh.set_point(*v_it, p);
    }

	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();
    return mesh;
}

AKMesh add_noise(AKMesh & _mesh, AKNumber intensity)
{
	AKMesh mesh = _mesh;
	// compute average length of mesh
	AKNumber average_length = 0.0;
	for (AKMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
		average_length += mesh.calc_edge_length(*e_it);
	AKNumber edge_numbers = (AKNumber)mesh.n_edges();
	average_length /= edge_numbers;
	AKNumber noise_level = intensity;
	AKNumber standard_deviation = average_length * noise_level;

	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();

	vector<AKNumber> GaussianNumbers;

	generate_random_gaussian_numbers(0, standard_deviation, (int)mesh.n_vertices(), GaussianNumbers);

	for (AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++){

		AKMesh::Point p = mesh.point(*v_it) + mesh.normal(*v_it) * GaussianNumbers[v_it->idx()];
		mesh.set_point(*v_it, p);
	}

	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();
	return mesh;
}

AKMesh::Normal generate_random_direction()
{
    AKNumber x, y, z, length;
    do{
        x = -1 + 2 * (AKNumber)rand() / (AKNumber) RAND_MAX;
        y = -1 + 2 * (AKNumber)rand() / (AKNumber) RAND_MAX;
        length = x * x + y * y;
    }while(length>1);

    const AKNumber r = 2 * std::sqrt(1 - length);

    x *= r;
    y *= r;
    z = 1 - 2 * length;

    return AKMesh::Normal(x, y, z);
}

void generate_random_directions(int number, std::vector<AKMesh::Normal> &RandomDirections)
{
    RandomDirections.resize(number, AKMesh::Normal(0.0, 0.0, 0.0));

    srand((unsigned int)time(NULL));
    for(int i = 0; i < number; i++){
        RandomDirections[i] = generate_random_direction();
    }
}

AKMesh add_noise_random_direction(AKMesh & _mesh)
{
    AKMesh mesh = _mesh;
    // compute average length of mesh
    AKNumber average_length = 0.0;
    for(AKMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
        average_length += mesh.calc_edge_length(*e_it);
    AKNumber edge_numbers = (AKNumber)mesh.n_edges();
    average_length /= edge_numbers;
    AKNumber noise_level = default_noise_level;
    AKNumber standard_deviation = average_length * noise_level;

    mesh.request_face_normals();
    mesh.request_vertex_normals();
    mesh.update_normals();

    vector<AKNumber> GaussianNumbers;
    vector<AKMesh::Normal> RandomDirections;

    generate_random_gaussian_numbers(0, standard_deviation, (int)mesh.n_vertices(), GaussianNumbers);
    generate_random_directions((int)mesh.n_vertices(),RandomDirections);

    for(AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++){

        AKMesh::Point p = mesh.point(*v_it) + RandomDirections[v_it->idx()] * GaussianNumbers[v_it->idx()];
        mesh.set_point(*v_it, p);
    }
    return mesh;
}
