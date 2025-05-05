#include "core/denoising.h"

/////////////////////////////////////////////////
/// Auxiliar functions
/////////////////////////////////////////////////

bool hasIND(AKMesh::Point & p)
{
	if (p[0] != p[0] || p[1] != p[1] || p[2] != p[2])
		return true;
	else return false;
}

void update_vertex_positions(AKMesh &mesh, vector<AKMesh::Normal> &filtered_normals, int iteration_number, bool fixed_boundary)
{
	clock_t begin = clock();
    vector<AKMesh::Point> new_points(mesh.n_vertices());
    vector<AKMesh::Point> centroids;

    for(int iter = 0; iter < iteration_number; iter++)
    {
        compute_all_face_centroids(mesh, centroids);
        for(AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
        {
            AKMesh::Point p = mesh.point(*v_it);
            if(fixed_boundary && mesh.is_boundary(*v_it))
            {
                new_points.at(v_it->idx()) = p;
            }
            else
            {
                AKNumber face_num = 0.0f;
                AKMesh::Point temp_point(0.0f, 0.0f, 0.0f);
                for(AKMesh::VertexFaceIter vf_it = mesh.vf_iter(*v_it); vf_it.is_valid(); vf_it++)
                {
                    AKMesh::Normal temp_normal = filtered_normals[vf_it->idx()];
                    AKMesh::Point temp_centroid = centroids[vf_it->idx()];
                    temp_point += temp_normal * (temp_normal | (temp_centroid - p));
                    face_num++;
                }
				if (face_num != 0)
					p += temp_point / face_num;
                if (!hasIND(p))
                {
                    new_points.at(v_it->idx()) = p;
                }
                else new_points.at(v_it->idx()) = p;
            }
        }
        for(AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
            mesh.set_point(*v_it, new_points[v_it->idx()]);
    }
	clock_t end = clock();
	AKNumber elapsed_secs = AKNumber(end - begin) / CLOCKS_PER_SEC;
	std::cout << " time for vertex updating: " << elapsed_secs << endl;
}

AKNumber compute_sigma_spatial(AKMesh &mesh, vector<AKMesh::Point> &face_centroids, AKNumber sigma_c_scalar)
{
    AKNumber sigma_c = 0.0f;
    AKNumber num = 0.0f;
    for(AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        AKMesh::Point ci = face_centroids[f_it->idx()];
        for(AKMesh::FaceFaceIter ff_it = mesh.ff_iter(*f_it); ff_it.is_valid(); ff_it++)
        {
            AKMesh::Point cj = face_centroids[ff_it->idx()];
            sigma_c += (ci - cj).length();
            num++;
        }
    }
    sigma_c *= sigma_c_scalar / num;
    return sigma_c;
}

AKNumber compute_radius(AKMesh &mesh , AKNumber scalar )
{
    vector<AKMesh::Point> centroids;
    compute_all_face_centroids(mesh, centroids);

    AKNumber radius = 0.0f;
    AKNumber num = 0.0f;
    for(AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        AKMesh::Point fi = centroids[f_it->idx()];
        for(AKMesh::FaceFaceIter ff_it = mesh.ff_iter(*f_it); ff_it.is_valid(); ff_it++)
        {
            AKMesh::Point fj = centroids[ff_it->idx()];
            radius += (fj - fi).length();
            num++;
        }
    }
	AKNumber res = 0.0f;
	if (num != 0.0f)
		res = radius * scalar / num;
	return res;
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/// Isotropic mesh denoising algorithms
/////////////////////////////////////////////////
/////////////////////////////////////////////////

/////////////////////////////////////////////////
/// Uniform Laplacian smoothing
/////////////////////////////////////////////////

AKMesh uniform_laplacian_smoothing(AKMesh & _mesh, int iteration_number, AKNumber scale)
{
	AKMesh mesh = _mesh;
	vector<AKMesh::Point> displacement_points;
	displacement_points.resize(mesh.n_vertices());
	for (int iter = 0; iter < iteration_number; iter++)
	{
		for (AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
		{
			vector<AKMesh::VertexHandle> vertex_neighbors;
			compute_vertex_neighbors(_mesh, *v_it, 1, vertex_neighbors);
			AKNumber weight = 1.0f / ((AKNumber)vertex_neighbors.size());
			AKMesh::Point sum(0.0f, 0.0f, 0.0f);
			for (int i = 0; i<(int)(vertex_neighbors.size()); i++)
				sum = sum + ((mesh.point(vertex_neighbors[i]) - mesh.point(*v_it))*weight);
			displacement_points[v_it->idx()] = sum;
		}
		for (AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
			mesh.set_point(*v_it, displacement_points[v_it->idx()] * scale + mesh.point(*v_it));
	}
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();
	return mesh;
}

AKMesh uniform_laplacian_smoothing(AKMesh & _mesh, int iteration_number, AKNumber scale, vector<size_t> & vertex_ids)
{
	AKMesh mesh = _mesh;
	vector<AKMesh::Point> displacement_points;
	displacement_points.resize(mesh.n_vertices());
	for (int iter = 0; iter < iteration_number; iter++)
	{
		for (size_t i = 0; i < vertex_ids.size(); i++)
		{
			AKMesh::VertexHandle vh((int)vertex_ids[i]);
			vector<AKMesh::VertexHandle> vertex_neighbors;
			compute_vertex_neighbors(_mesh, vh, 1, vertex_neighbors);
			AKNumber weight = 1.0f / ((AKNumber)vertex_neighbors.size());
			AKMesh::Point sum(0.0f, 0.0f, 0.0f);
			for (int i = 0; i<(int)(vertex_neighbors.size()); i++)
				sum = sum + ((mesh.point(vertex_neighbors[i]) - mesh.point(vh))*weight);
			displacement_points[vh.idx()] = sum;
		}
		for (size_t i = 0; i < vertex_ids.size(); i++)
		{
			AKMesh::VertexHandle vh((int)vertex_ids[i]);
			mesh.set_point(vh, displacement_points[vh.idx()] * scale + mesh.point(vh));
		}
	}
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();
	return mesh;
}

/////////////////////////////////////////////////
/// HC Laplacian smoothing (Vollmer et al.)
/////////////////////////////////////////////////

AKMesh laplacian_smoothing_HC(AKMesh & _mesh,int iteration_number,AKNumber alpha,AKNumber beta)
{
    AKMesh mesh = _mesh;
    vector<AKMesh::Point> original_points;
    compute_all_points(mesh,original_points);
    AKMesh temp_mesh_p;
    vector<AKMesh::Point> temp_points_p(original_points);
	vector<AKMesh::Point> uniform_displacement_points;
	uniform_displacement_points.resize(mesh.n_vertices());
    for(int iter = 0; iter < iteration_number; iter++)
    {
        vector<AKMesh::Point> temp_points_q(temp_points_p);
		for (AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
		{
			vector<AKMesh::VertexHandle> vertex_neighbors;
			compute_vertex_neighbors(_mesh, *v_it, 1, vertex_neighbors);
			if (vertex_neighbors.size()>0)
			{
				AKNumber weight = 1.0f / ((AKNumber)vertex_neighbors.size());
				AKMesh::Point sum(0.0f, 0.0f, 0.0f);
				for (int i = 0; i<(int)(vertex_neighbors.size()); i++)
					sum = sum + ((mesh.point(vertex_neighbors[i]) - mesh.point(*v_it))*weight);
				uniform_displacement_points[v_it->idx()] = sum;
				temp_points_p[v_it->idx()] += sum;
			}
			else
				uniform_displacement_points[v_it->idx()] = AKMesh::Point(0,0,0);
		}
		for (AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
			temp_points_p[v_it->idx()]=uniform_displacement_points[v_it->idx()] + mesh.point(*v_it);
        vector<AKMesh::Point> temp_points_b(temp_points_p.size());
        for(size_t i = 0;i<temp_points_p.size();i++)
            temp_points_b[i]=temp_points_p[i]-(original_points[i]*alpha+(1.0f-alpha)*temp_points_q[i]);
        for(AKMesh::VertexIter v_it = mesh.vertices_begin();v_it!=mesh.vertices_end();v_it++)
        {
            AKMesh::Point t_sum(0,0,0);
            AKNumber num = 0;
            for(AKMesh::VertexVertexIter vv_it = mesh.vv_iter(*v_it);vv_it.is_valid();vv_it++,num+=1.0f)
                t_sum+=temp_points_b[vv_it->idx()];

			AKMesh::Point disp = (beta*temp_points_b[v_it->idx()] + ((1.0f - beta) / num)*t_sum);
			if (!hasIND(disp))
				temp_points_p[v_it->idx()]=temp_points_p[v_it->idx()]-disp;
        }
    }
    for(AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
        mesh.set_point(*v_it, temp_points_p[v_it->idx()]);
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();
    return mesh;
}

AKMesh laplacian_smoothing_HC(AKMesh & _mesh, int iteration_number, AKNumber alpha, AKNumber beta, vector<size_t> & vertex_ids)
{
	AKMesh mesh = _mesh;
	vector<AKMesh::Point> original_points;
	compute_all_points(mesh, original_points);
	AKMesh temp_mesh_p;
	vector<AKMesh::Point> temp_points_p(original_points);
	vector<AKMesh::Point> uniform_displacement_points;
	uniform_displacement_points.resize(mesh.n_vertices());
	for (int iter = 0; iter < iteration_number; iter++)
	{
		vector<AKMesh::Point> temp_points_q(temp_points_p);
		for (size_t i = 0; i < vertex_ids.size(); i++)
		{
			AKMesh::VertexHandle vh((int)vertex_ids[i]);
			vector<AKMesh::VertexHandle> vertex_neighbors;
			compute_vertex_neighbors(_mesh, vh, 1, vertex_neighbors);
			if (vertex_neighbors.size()>0)
			{
				AKNumber weight = 1.0f / ((AKNumber)vertex_neighbors.size());
				AKMesh::Point sum(0.0f, 0.0f, 0.0f);
				for (int i = 0; i<(int)(vertex_neighbors.size()); i++)
					sum = sum + ((mesh.point(vertex_neighbors[i]) - mesh.point(vh))*weight);
				uniform_displacement_points[vh.idx()] = sum;
				temp_points_p[vh.idx()] += sum;
			}
			else
				uniform_displacement_points[vh.idx()] = AKMesh::Point(0, 0, 0);
		}
		for (size_t i = 0; i < vertex_ids.size(); i++)
		{
			AKMesh::VertexHandle vh((int)vertex_ids[i]);
			temp_points_p[vh.idx()] = uniform_displacement_points[vh.idx()] + mesh.point(vh);
		}
		vector<AKMesh::Point> temp_points_b(temp_points_p.size());
		for (size_t i = 0; i<temp_points_p.size(); i++)
			temp_points_b[i] = temp_points_p[i] - (original_points[i] * alpha + (1.0f - alpha)*temp_points_q[i]);
		for (size_t i = 0; i < vertex_ids.size(); i++)
		{
			AKMesh::VertexHandle vh((int)vertex_ids[i]);
			AKMesh::Point t_sum(0, 0, 0);
			AKNumber num = 0;
			for (AKMesh::VertexVertexIter vv_it = mesh.vv_iter(vh); vv_it.is_valid(); vv_it++, num += 1.0f)
				t_sum += temp_points_b[vv_it->idx()];

			AKMesh::Point disp = (beta*temp_points_b[vh.idx()] + ((1.0f - beta) / num)*t_sum);
			if (!hasIND(disp))
				temp_points_p[vh.idx()] = temp_points_p[vh.idx()] - disp;
		}
	}
	for (size_t i = 0; i < vertex_ids.size(); i++)
	{
		AKMesh::VertexHandle vh((int)vertex_ids[i]);
		mesh.set_point(vh, temp_points_p[vh.idx()]);
	}
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();
	return mesh;
}

/////////////////////////////////////////////////
/// Bilateral normal filtering for mesh denoising (Zheng et al.)
/////////////////////////////////////////////////

void bilateral_filtering(AKMesh &mesh, int normal_iteration_number, AKNumber sigma_c_scalar, AKNumber sigma_s, vector<AKMesh::Normal> &filtered_normals)
{
    filtered_normals.resize(mesh.n_faces());

    vector< vector<AKMesh::FaceHandle> > all_face_neighbors;
    compute_all_vertex_based_face_neighbors(mesh,0,all_face_neighbors);
    vector<AKMesh::Normal> previous_normals;
    compute_all_face_normals(mesh, previous_normals);
    vector<AKNumber> face_areas;
    compute_all_face_areas(mesh, face_areas);
    vector<AKMesh::Point> face_centroids;
    compute_all_face_centroids(mesh, face_centroids);

    AKNumber sigma_c = compute_sigma_spatial(mesh, face_centroids, sigma_c_scalar);

    for(int it = 0; it < normal_iteration_number; it++)
    {
        for(AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
        {
            int face_idx = f_it->idx();
			AKMesh::Normal ni = previous_normals[face_idx];
			AKMesh::Point ci = face_centroids[face_idx];
			vector<AKMesh::FaceHandle> face_neighbors = all_face_neighbors[face_idx];
            AKMesh::Normal temp_normal(0.0, 0.0, 0.0);
            AKNumber weight_sum = 0.0;
			for (int i = 0; i < (int)face_neighbors.size(); i++)
            {
                int neighbor_idx = face_neighbors[i].idx();
				AKMesh::Normal nj = previous_normals[neighbor_idx];
				AKMesh::Point cj = face_centroids[neighbor_idx];

                AKNumber spatial_distance = (ci - cj).length();
                AKNumber spatial_weight = compute_gaussian_weight(spatial_distance,sigma_c);
                AKNumber range_distance = (ni - nj).length();
                AKNumber range_weight = compute_gaussian_weight(range_distance,sigma_s);

				AKNumber weight = face_areas[neighbor_idx] * spatial_weight * range_weight;
                weight_sum += weight;
                temp_normal += nj * weight;
            }
            temp_normal /= weight_sum;
            temp_normal.normalize_cond();
			filtered_normals[face_idx] = temp_normal;
        }
        previous_normals = filtered_normals;
    }
}

AKMesh bilateral_normal_filtering(AKMesh &_mesh, int normal_iteration_number, int vertex_iteration_number, AKNumber sigma_c_scalar,  AKNumber sigma_s )
{
    AKMesh mesh = _mesh;
    vector<AKMesh::Normal> filtered_normals;
    bilateral_filtering(mesh, normal_iteration_number, sigma_c_scalar, sigma_s, filtered_normals);
    update_vertex_positions(mesh, filtered_normals, vertex_iteration_number, true);
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();
    return mesh;
}


/////////////////////////////////////////////////
/// Guided mesh normal filtering (Zhang et al.)
/////////////////////////////////////////////////

void compute_all_face_neighbors_guided_bilateral_normal_filtering(AKMesh &mesh, FaceNeighborType face_neighbor_type, AKNumber radius, bool include_central_face,
    vector<vector<AKMesh::FaceHandle> > &all_face_neighbors)
{
    vector<AKMesh::FaceHandle> face_neighbors;
	for (AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		if (face_neighbor_type == kVertexBased)
            compute_vertex_based_face_neighbors(mesh, *f_it, face_neighbors);
		else if (face_neighbor_type == kRadiusBased)
            compute_radius_based_face_neighbors(mesh, *f_it, radius, face_neighbors);

		if (include_central_face)
            face_neighbors.push_back(*f_it);
        all_face_neighbors[f_it->idx()] = face_neighbors;
	}
}

void compute_all_guidance_signal_neighbors(AKMesh &mesh, vector<vector<AKMesh::FaceHandle> > &all_guided_neighbors)
{
    vector<AKMesh::FaceHandle> face_neighbors;
	for (AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
        compute_vertex_based_face_neighbors(mesh, *f_it, face_neighbors);
        face_neighbors.push_back(*f_it);
        all_guided_neighbors[f_it->idx()] = face_neighbors;
	}
}

void compute_face_neighbors_inner_edges(AKMesh &mesh, vector<AKMesh::FaceHandle> &face_neighbors, vector<AKMesh::EdgeHandle> &inner_edges)
{
    inner_edges.clear();
    vector<bool> edge_flags((int)mesh.n_edges(), false);
    vector<bool> face_flags((int)mesh.n_faces(), false);

    for (int i = 0; i < (int)face_neighbors.size(); i++)
        face_flags[face_neighbors[i].idx()] = true;

    for (int i = 0; i < (int)face_neighbors.size(); i++)
	{
        for (AKMesh::FaceEdgeIter fe_it = mesh.fe_iter(face_neighbors[i]); fe_it.is_valid(); fe_it++)
		{
            if ((!edge_flags[fe_it->idx()]) && (!mesh.is_boundary(*fe_it)))
			{
                edge_flags[fe_it->idx()] = true;
				AKMesh::HalfedgeHandle heh = mesh.halfedge_handle(*fe_it, 0);
				AKMesh::FaceHandle f = mesh.face_handle(heh);
				AKMesh::HalfedgeHandle heho = mesh.opposite_halfedge_handle(heh);
				AKMesh::FaceHandle fo = mesh.face_handle(heho);
                if (face_flags[f.idx()] && face_flags[fo.idx()])
                    inner_edges.push_back(*fe_it);
			}
		}
	}
}

void compute_consistencies_and_mean_normals(AKMesh &mesh, vector<vector<AKMesh::FaceHandle> > &all_guided_neighbors,
    vector<AKNumber> &face_areas, vector<AKMesh::Normal> &normals,
    vector<pair<AKNumber, AKMesh::Normal> > &consistencies_and_mean_normals)
{
	const AKNumber epsilon = 1.0e-9f;

	for (AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		int index = f_it->idx();
        vector<AKMesh::FaceHandle> face_neighbors = all_guided_neighbors[index];
		AKNumber metric = 0.0f;
		AKMesh::Normal average_normal(0.0f, 0.0f, 0.0f);
		AKNumber maxdiff = -1.0f;

        for (int i = 0; i < (int)face_neighbors.size(); i++)
		{
            int index_i = face_neighbors[i].idx();
			AKNumber area_weight = face_areas[index_i];
			AKMesh::Normal ni = normals[index_i];
			average_normal += ni * area_weight;

            for (int j = i + 1; j < (int)face_neighbors.size(); j++)
			{
                int index_j = face_neighbors[j].idx();
				AKMesh::Normal nj = normals[index_j];
				AKNumber diff = normal_distance(ni, nj);

				if (diff > maxdiff)
				{
					maxdiff = diff;
				}
			}
		}

        vector<AKMesh::EdgeHandle> inner_edge_handles;
        compute_face_neighbors_inner_edges(mesh, face_neighbors, inner_edge_handles);
		AKNumber sum_tv = 0.0, max_tv = -1.0;
        for (int i = 0; i < (int)inner_edge_handles.size(); i++)
		{
            AKMesh::HalfedgeHandle heh = mesh.halfedge_handle(inner_edge_handles[i], 0);
			AKMesh::FaceHandle f = mesh.face_handle(heh);
			AKMesh::Normal n1 = normals[f.idx()];
			AKMesh::HalfedgeHandle heho = mesh.opposite_halfedge_handle(heh);
			AKMesh::FaceHandle fo = mesh.face_handle(heho);
			AKMesh::Normal n2 = normals[fo.idx()];
			AKNumber current_tv = normal_distance(n1, n2);
			max_tv = (current_tv > max_tv) ? current_tv : max_tv;
			sum_tv += current_tv;
		}
		average_normal.normalize_cond();
		metric = maxdiff * max_tv / (sum_tv + epsilon);

        consistencies_and_mean_normals[index] = make_pair(metric, average_normal);
	}
}

void compute_guided_normals(AKMesh &mesh, vector<vector<AKMesh::FaceHandle> > &all_guided_neighbors,
    vector<AKNumber> &face_areas, vector<AKMesh::Normal> &normals,
    vector<pair<AKNumber, AKMesh::Normal> > consistencies_and_mean_normals,
    vector<AKMesh::Normal> &guided_normals)
{
    compute_consistencies_and_mean_normals(mesh, all_guided_neighbors, face_areas, normals, consistencies_and_mean_normals);

	for (AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
        vector<AKMesh::FaceHandle> face_neighbors = all_guided_neighbors[f_it->idx()];
        AKNumber min_consistency = 1.0e8f;
		int min_idx = 0;
        for (int i = 0; i < (int)face_neighbors.size(); i++)
		{
            AKNumber current_consistency = consistencies_and_mean_normals[face_neighbors[i].idx()].first;
            if (min_consistency > current_consistency){
                min_consistency = current_consistency;
				min_idx = i;
			}
		}
        AKMesh::FaceHandle min_face_handle = face_neighbors[min_idx];
        guided_normals[f_it->idx()] = consistencies_and_mean_normals[min_face_handle.idx()].second;
	}
}

void guided_bilateral_normal_filtering(AKMesh &mesh, vector<AKMesh::Normal> &filtered_normals, AKNumber radius_scalar,
	AKNumber sigma_c_scalar, int normal_iteration_number, AKNumber sigma_s, int vertex_iteration_number)
{
	filtered_normals.resize((int)mesh.n_faces());
	bool include_central_face = 1;
    FaceNeighborType face_neighbor_type = kRadiusBased;
	AKNumber radius;
    radius = compute_radius(mesh, radius_scalar);
    vector<vector<AKMesh::FaceHandle> > all_face_neighbors((int)mesh.n_faces());
    compute_all_face_neighbors_guided_bilateral_normal_filtering(mesh, face_neighbor_type, radius, include_central_face, all_face_neighbors);
    vector<vector<AKMesh::FaceHandle> > all_guided_neighbors((int)mesh.n_faces());
    compute_all_guidance_signal_neighbors(mesh, all_guided_neighbors);
	compute_all_face_normals(mesh, filtered_normals);

	vector<AKNumber> face_areas((int)mesh.n_faces());
	vector<AKMesh::Point> face_centroids((int)mesh.n_faces());
	vector<AKMesh::Normal> previous_normals((int)mesh.n_faces());
	vector<AKMesh::Normal> guided_normals((int)mesh.n_faces());
    vector<pair<AKNumber, AKMesh::Normal> > consistencies_and_mean_normals((int)mesh.n_faces());
	for (int iter = 0; iter < normal_iteration_number; iter++)
	{
		compute_all_face_centroids(mesh, face_centroids);
		AKNumber sigma_c = compute_sigma_spatial(mesh, face_centroids, sigma_c_scalar);
		compute_all_face_areas(mesh, face_areas);
		compute_all_face_normals(mesh, previous_normals);

        compute_guided_normals(mesh, all_guided_neighbors, face_areas, previous_normals, consistencies_and_mean_normals, guided_normals);

		for (AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
		{
			int index = f_it->idx();
            const vector<AKMesh::FaceHandle> face_neighbors = all_face_neighbors[index];
			AKMesh::Normal filtered_normal(0.0f, 0.0f, 0.0f);
            for (int j = 0; j < (int)face_neighbors.size(); j++)
			{
                int current_face_index = face_neighbors[j].idx();

				AKNumber spatial_dis = (face_centroids[index] - face_centroids[current_face_index]).length();
				AKNumber spatial_weight = compute_gaussian_weight(spatial_dis, sigma_c);
				AKNumber range_dis = (guided_normals[index] - guided_normals[current_face_index]).length();
				AKNumber range_weight = compute_gaussian_weight(range_dis, sigma_s);

				filtered_normal += previous_normals[current_face_index] * (face_areas[current_face_index] * spatial_weight * range_weight);
			}
            if (face_neighbors.size())
				filtered_normals[index] = filtered_normal.normalize_cond();
		}
		update_vertex_positions(mesh, filtered_normals, vertex_iteration_number, false);
		//updateVertexPositions(mesh, guided_normals, vertex_iteration_number, false);
	}
	//filtered_normals = guided_normals;
}

AKMesh guided_denoising(AKMesh _mesh, int normal_iteration_number, int vertex_iteration_number, AKNumber sigma_c_scalar, AKNumber sigma_s, AKNumber radius_scalar)
{
	AKMesh mesh = _mesh;
	if (mesh.n_vertices() == 0)
		return mesh;
	vector<AKMesh::Normal> filtered_normals;
	guided_bilateral_normal_filtering(mesh, filtered_normals, radius_scalar, sigma_c_scalar, normal_iteration_number, sigma_s, vertex_iteration_number);
	update_vertex_positions(mesh, filtered_normals, vertex_iteration_number, true);
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();
	return mesh;
}

vector<AKMesh::Normal> adaptive_kernels_normal_filtering(AKMesh & mesh, int normal_iteration_number, int bilateral_iteration_number, int vertex_iteration_number, AKNumber sigma_c_scalar, AKNumber sigma_s, AKNumber alpha, AKNumber beta, AKNumber gamma, AKNumber delta, AKNumber k)
{
	clock_t begin = clock();

	vector<vector<pair<size_t, AKNumber>>> _patches_sparse;
	vector<vector<pair<size_t, AKNumber>>> _distance_matrix_sparse;
	vector<AKNumber> areas;
	vector<AKMesh::Point> centroids;
	vector<AKMesh::Normal> normals;

	AKNumber avg_edge_length = compute_average_edge_length(mesh);
	AKNumber distance_limit = avg_edge_length * 2.0;

	compute_distance_matrix(mesh, distance_limit, _distance_matrix_sparse, k);

	compute_all_face_areas(mesh, areas);
	compute_all_face_centroids(mesh, centroids);
	compute_all_face_normals(mesh, normals);

	compute_all_adaptive_kernels(mesh, areas, centroids, normals, _distance_matrix_sparse, alpha, beta, gamma, delta, k, _patches_sparse);

	AKNumber sigma_c = compute_sigma_spatial(mesh, centroids, sigma_c_scalar);

	vector<AKMesh::Normal> t_normals(_patches_sparse.size());

	//computar normales
	for (int it = 0; it < normal_iteration_number; it++)
	{
		for (int i = 0; i < _patches_sparse.size(); i++)
		{
			AKMesh::Normal avg_normal = AKMesh::Normal(0, 0, 0);
			for (int j = 0; j < _patches_sparse[i].size(); j++)
			{
				avg_normal += normals[_patches_sparse[i][j].first] * _patches_sparse[i][j].second * areas[_patches_sparse[i][j].first];
			}
			avg_normal.normalize_cond();
			t_normals[i] = avg_normal;
		}
		normals = t_normals;
	}
	vector<AKMesh::Normal> adaptive_kernel_normals = t_normals;
	for (int it = 0; it < bilateral_iteration_number; it++)
	{
		for (int i = 0; i < _patches_sparse.size(); i++)
		{
			AKMesh::Normal ni = t_normals[i];
			AKMesh::Point ci = centroids[i];
			AKMesh::Normal temp_normal(0.0, 0.0, 0.0);
			AKNumber weight_sum = 0.0;
			AKMesh::Normal avg_normal = AKMesh::Normal(0, 0, 0);
			for (int j = 0; j < _patches_sparse[i].size(); j++)
			{
				int idx_j = _patches_sparse[i][j].first;
				AKMesh::Normal nj = t_normals[idx_j];
				AKMesh::Point cj = centroids[idx_j];

				AKNumber spatial_distance = (ci - cj).length();
				AKNumber spatial_weight = compute_gaussian_weight(spatial_distance, sigma_c);
				AKNumber range_distance = (ni - nj).length();
				AKNumber range_weight = compute_gaussian_weight(range_distance, sigma_s);
				//cout << _patches_sparse[i][j].second << " " << spatial_weight << " " << range_weight << endl;
				//AKNumber weight = areas[idx_j] * _patches_sparse[i][j].second * spatial_weight * range_weight;
				AKNumber weight = areas[idx_j] * range_weight * spatial_weight;
				//AKNumber weight = _patches_sparse[i][j].second;
				weight_sum += weight;
				temp_normal += nj * weight;


				//avg_normal += normals[_patches_sparse[i][j].first] * _patches_sparse[i][j].second * areas[_patches_sparse[i][j].first];
			}
			temp_normal /= weight_sum;
			temp_normal.normalize_cond();
			adaptive_kernel_normals[i] = temp_normal;
		}
		t_normals = adaptive_kernel_normals;
	}

	clock_t end = clock();
	AKNumber elapsed_secs = AKNumber(end - begin) / CLOCKS_PER_SEC;
	std::cout << " time for normal calculation: " << elapsed_secs << endl;

	return adaptive_kernel_normals;
}


AKMesh adaptive_kernels_denoising(AKMesh _mesh, int normal_iteration_number, int bilateral_iteration_number, int vertex_iteration_number, int external_iteration_number, AKNumber sigma_c_scalar, AKNumber sigma_s, AKNumber alpha, AKNumber beta, AKNumber gamma, AKNumber delta, AKNumber k)
{
	clock_t begin = clock();

	AKMesh mesh = _mesh;
	if (mesh.n_vertices() == 0)
		return mesh;
	vector<AKMesh::Normal> filtered_normals;
	for (int i = 0; i < external_iteration_number; i++)
	{
		filtered_normals = adaptive_kernels_normal_filtering(mesh, normal_iteration_number, bilateral_iteration_number, vertex_iteration_number, sigma_c_scalar, sigma_s, alpha, beta, gamma, delta, k);
		update_vertex_positions(mesh, filtered_normals, vertex_iteration_number, true);
	}
	//filtered_normals = updateRobustPatchNormals4(mesh, normal_iteration_number, vertex_iteration_number, sigma_c_scalar, sigma_s, alpha, beta, gamma, delta, k);
	//updateVertexPositions(mesh, filtered_normals, vertex_iteration_number, true);
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();

	clock_t end = clock();
	AKNumber elapsed_secs = AKNumber(end - begin) / CLOCKS_PER_SEC;
	std::cout << " time for all denoising algorithm: " << elapsed_secs << endl;

	return mesh;
}