#include "core/adaptive_kernel.h"

struct {
	bool operator()(pair<size_t, AKNumber> & p1, pair<size_t, AKNumber> & p2) const
	{
		return (p1.first < p2.first);
	}
} customLess;

struct {
	bool operator()(pair<size_t, AKNumber> & p1, pair<size_t, AKNumber> & p2) const
	{
		return (p1.second < p2.second);
	}
} customLess2;

void sort_distance_matrix(vector<vector<pair<size_t, AKNumber>>> & _sparse_matrix)
{
	for (size_t i = 0; i < _sparse_matrix.size(); i++)
		std::sort(_sparse_matrix[i].begin(), _sparse_matrix[i].end(), customLess);
}

void filter_distance_matrix(vector<vector<pair<size_t, AKNumber>>> & _sparse_matrix, int max_n_var)
{
	for (size_t i = 0; i < _sparse_matrix.size(); i++)
	{
		if (_sparse_matrix[i].size()>max_n_var)
		{
			std::sort(_sparse_matrix[i].begin(), _sparse_matrix[i].end(), customLess2);
			_sparse_matrix[i].resize(max_n_var);
		}
		std::sort(_sparse_matrix[i].begin(), _sparse_matrix[i].end(), customLess);
	}

}

AKNumber get_distance_matrix_value_at(vector<vector<pair<size_t, AKNumber>>> & _sparse_matrix, int row, int column, AKNumber default_value)
{
	int first = 0;
	int last = _sparse_matrix[row].size() - 1;
	int middle = (first + last) / 2;
	while (first <= last)
	{
		if (_sparse_matrix[row][middle].first < column)
		{
			first = middle + 1;

		}
		else if (_sparse_matrix[row][middle].first == column)
		{
			return _sparse_matrix[row][middle].second;
		}
		else
		{
			last = middle - 1;
		}
		middle = (first + last) / 2;
	}
	//if (first > last)
	return default_value;
}

void compute_distance_matrix(AKMesh & _mesh, AKNumber _max_dist, vector<vector<pair<size_t, AKNumber>>> & _distance_matrix_sparse, int _max_num_vars)
{
	clock_t begin = clock();

	_distance_matrix_sparse.clear();
	_distance_matrix_sparse.resize(_mesh.n_faces());
	for (AKMesh::FaceIter f_it = _mesh.faces_begin(); f_it != _mesh.faces_end(); f_it++)
	{
		AKMesh::Point ci = _mesh.calc_face_centroid(*f_it);
		vector<bool> flag((int)_mesh.n_faces(), false);
		flag[f_it->idx()] = true;
		queue<AKMesh::FaceHandle> queue_face_handle;
		queue<AKNumber> queue_face_distance;
		queue_face_handle.push(*f_it);
		queue_face_distance.push(0.0f);
		_distance_matrix_sparse[f_it->idx()].push_back(pair<size_t, AKNumber>(f_it->idx(), 0.0f));
		vector<AKMesh::FaceHandle> temp_face_neighbors;
		while (!queue_face_handle.empty())
		{
			AKMesh::FaceHandle temp_face_handle_queue = queue_face_handle.front();
			AKNumber temp_distance = queue_face_distance.front();
			if (temp_face_handle_queue != *f_it)
				_distance_matrix_sparse[f_it->idx()].push_back(pair<size_t, AKNumber>(temp_face_handle_queue.idx(), temp_distance));
			queue_face_handle.pop();
			queue_face_distance.pop();
			compute_vertex_based_face_neighbors(_mesh, temp_face_handle_queue, temp_face_neighbors);
			for (int i = 0; i < (int)temp_face_neighbors.size(); i++)
			{
				AKMesh::FaceHandle temp_face_handle = temp_face_neighbors[i];
				if (!flag[temp_face_handle.idx()])
				{
					AKMesh::Point cj = _mesh.calc_face_centroid(temp_face_handle);
					AKNumber distance = (ci - cj).length();
					if (distance <= _max_dist)
					{
						queue_face_handle.push(temp_face_handle);
						queue_face_distance.push(distance);
					}
					flag[temp_face_handle.idx()] = true;
				}
			}
		}
		if (_distance_matrix_sparse[f_it->idx()].size()<10)
		{
			_distance_matrix_sparse[f_it->idx()].clear();
			_distance_matrix_sparse[f_it->idx()].push_back(pair<size_t, AKNumber>(f_it->idx(), 0.0f));
			compute_vertex_based_face_neighbors(_mesh, *f_it, temp_face_neighbors);
			for (int i = 0; i < temp_face_neighbors.size(); i++)
			{
				AKMesh::Point cj = _mesh.calc_face_centroid(temp_face_neighbors[i]);
				AKNumber distance = (ci - cj).length();
				_distance_matrix_sparse[f_it->idx()].push_back(pair<size_t, AKNumber>(temp_face_neighbors[i].idx(), distance));
			}

		}
	}
	filter_distance_matrix(_distance_matrix_sparse, _max_num_vars);

	clock_t end = clock();
	AKNumber elapsed_secs = AKNumber(end - begin) / CLOCKS_PER_SEC;
	std::cout << " time for distances: " << elapsed_secs << endl;

	//sort(_distance_matrix_sparse);
}

void writeData(AKMatrix & H, AKMatrix & A, AKNumber area_const)
{
	int num_faces = H.rows();
	ofstream matrix_h;
	ofstream constraints;
	matrix_h.open("mat.txt");
	constraints.open("const.txt");
	constraints << area_const << " ";
	for (int i = 0; i < num_faces; i++)
	{
		for (int j = 0; j < num_faces; j++)
		{
			matrix_h << H(i, j) << " ";
		}
		matrix_h << endl;
		constraints << A(i, i) << " ";
	}
	matrix_h.close();
	constraints.close();
}

void compute_all_adaptive_kernels(AKMesh & _mesh, vector<AKNumber> & areas, vector<AKMesh::Point> & centroids, vector<AKMesh::Normal> & normals, vector<vector<pair<size_t, AKNumber>>> & _distance_matrix_sparse, AKNumber alpha, AKNumber beta, AKNumber gamma, AKNumber delta, AKNumber k, vector<vector<pair<size_t, AKNumber>>> & _patches_sparse)
{
	_patches_sparse.clear();
	_patches_sparse.resize(_mesh.n_faces());

	clock_t begin = clock();

#pragma omp parallel
#pragma omp for
	for (int f_idx = 0; f_idx < _mesh.n_faces(); f_idx++)
	{
		clock_t patch_begin = clock();

		compute_adaptive_kernel(_mesh, f_idx, areas, centroids, normals, _distance_matrix_sparse[f_idx], alpha, beta, gamma, delta, k, _patches_sparse[f_idx]);

		clock_t patch_end = clock();
		AKNumber patch_elapsed_secs = AKNumber(patch_end - patch_begin) / CLOCKS_PER_SEC;
		//std::cout << " time for patch " << f_idx << " : " << patch_elapsed_secs << endl;
		//std::cout << " num vars: " << _patches_sparse[f_idx].size() << endl;

	}
	clock_t end = clock();
	AKNumber elapsed_secs = AKNumber(end - begin) / CLOCKS_PER_SEC;
	std::cout << " time for all patches: " << elapsed_secs << endl;
}

void compute_adaptive_kernel(AKMesh & _mesh, int _face_index, vector<AKNumber> & areas, vector<AKMesh::Point> & centroids, vector<AKMesh::Normal> & normals,
	vector<pair<size_t, AKNumber>> & _distance_function, AKNumber alpha, AKNumber beta, AKNumber gamma, AKNumber delta, AKNumber k, vector<pair<size_t, AKNumber>> & _patch_function)
{
	_patch_function.clear();
	AKNumber t_inf = 999999.0f;

	int num_faces = _distance_function.size();
	AKMatrix error_matrix(num_faces, num_faces);
	error_matrix.fill(0);
	AKMatrix distance_matrix(num_faces, num_faces);
	distance_matrix.fill(0);
	AKMatrix gradient_matrix(num_faces, num_faces);
	gradient_matrix.fill(0);
	AKMatrix area_matrix(num_faces, num_faces);
	area_matrix.fill(0);
	AKMatrix central_face_normal_difference_matrix(num_faces, num_faces);
	central_face_normal_difference_matrix.fill(0);
	AKMatrix linear_term(num_faces, 1);
	linear_term.fill(0);

	AKNumber total_area = 0.0f;

	for (size_t i = 0; i < num_faces; i++)
	{
		size_t idx1 = _distance_function[i].first;
		AKNumber total_edge_length = 0.0f;
		for (size_t j = 0; j < num_faces; j++)
		{
			size_t idx2 = _distance_function[j].first;
			AKNumber normal_dif = (normals[idx1].normalize() - normals[idx2].normalize()).length();
			error_matrix(i, j) = normal_dif;
			error_matrix(j, i) = normal_dif;


			AKNumber t_edge_length = 0.0f;
			for (AKMesh::FaceEdgeIter fe_it = _mesh.fe_iter(AKMesh::FaceHandle(idx1)); fe_it.is_valid(); fe_it++)
			{
				AKMesh::HalfedgeHandle heh = _mesh.halfedge_handle(*fe_it, 0);
				AKMesh::FaceHandle f = _mesh.face_handle(heh);
				AKMesh::HalfedgeHandle heho = _mesh.opposite_halfedge_handle(heh);
				AKMesh::FaceHandle fo = _mesh.face_handle(heho);
				if (f.idx() == idx2 || fo.idx() == idx2)
					t_edge_length = _mesh.calc_edge_length(heh);
			}
			gradient_matrix(i, j) = -t_edge_length;
			gradient_matrix(j, i) = -t_edge_length;
			total_edge_length += t_edge_length;
		}
		AKNumber distance = _distance_function[i].second;
		distance_matrix(i, i) = distance;
		area_matrix(i, i) = areas[idx1];
		gradient_matrix(i, i) = total_edge_length;
		AKNumber central_face_normal_difference = (normals[_face_index].normalize() - normals[idx1].normalize()).length();
		central_face_normal_difference_matrix(i, i) = central_face_normal_difference;
		linear_term(i, 0) = (beta * (distance * areas[_face_index] * areas[idx1])) + (delta* (central_face_normal_difference * areas[_face_index] * areas[idx1]));

		total_area += areas[idx1];
	}


	AKMatrix H(num_faces, num_faces);
	H.fill(0);
	AKMatrix NH(num_faces, num_faces);
	NH.fill(0);
	AKMatrix error_matrix_t = error_matrix.transpose();
	AKMatrix gradient_matrix_t = gradient_matrix.transpose();

	//H = alpha * (area_matrix * (error_matrix_t * error_matrix) * area_matrix) + beta * (distance_matrix * distance_matrix) + gamma * (gradient_matrix_t * gradient_matrix);
	//H = alpha * (area_matrix * (error_matrix) * area_matrix) + beta * (distance_matrix * distance_matrix) + gamma * (gradient_matrix_t * gradient_matrix) + delta * (area_matrix * central_face_normal_difference_matrix * area_matrix);
	H = alpha * (area_matrix * (error_matrix)* area_matrix) + gamma * (gradient_matrix_t * gradient_matrix);

	NH = (H + H.transpose())*0.5f;

	AKNumber area_constraint = total_area * 0.2f;

	//writeData(NH, area_matrix, area_constraint);


	//cout << "area: " << total_area * 0.2f << endl;

	//quadratic_programming_solver(s)
	vector<AKNumber> solution;
	quadratic_programming_solver(num_faces, linear_term, NH, area_matrix, area_constraint, solution);
	for (int j = 0; j < solution.size(); j++)
	{
		_patch_function.push_back(pair<size_t, AKNumber>(_distance_function[j].first, solution[j]));
	}

}
