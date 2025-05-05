#include "core/metric.h"

const AKNumber average_edge_length_factor = 2.0; //KD-Tree nearest neighbor tolerance radius = factor*average_edge_length

void generate_vertex_point_cloud(AKMesh & mesh, PointCloud &pointCloud)
{
	pointCloud.pts.resize(mesh.n_vertices());
	for (AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		pointCloud.pts[v_it->idx()].x = mesh.point(*v_it)[0];
		pointCloud.pts[v_it->idx()].y = mesh.point(*v_it)[1];
		pointCloud.pts[v_it->idx()].z = mesh.point(*v_it)[2];
	}
}

void generate_face_centroid_point_cloud(AKMesh & mesh, PointCloud &pointCloud)
{
	pointCloud.pts.resize(mesh.n_faces());
	vector<AKMesh::Point> face_centroids;
	compute_all_face_centroids(mesh, face_centroids);

	for (size_t i = 0; i<face_centroids.size(); i++)
	{
		pointCloud.pts[i].x = face_centroids[i][0];
		pointCloud.pts[i].y = face_centroids[i][1];
		pointCloud.pts[i].z = face_centroids[i][2];
	}
}



AKNumber compute_area_error(AKNumber desiredArea, AKMesh & mesh)
{
	AKNumber dif = abs(compute_area(mesh) - desiredArea);
	return dif / desiredArea;
}



AKNumber compute_volume_error(AKNumber desiredVol, AKMesh & mesh)
{
	AKNumber dif = abs(compute_volume(mesh) - desiredVol);
	return dif / desiredVol;
}

AKNumber point_triangle_dist(AKMesh & _mesh, AKMesh::FaceHandle fh, AKMesh::Point _point)
{
	AKMesh::Point _p = _point;
	AKMesh::FaceVertexIter fv_iter = _mesh.fv_iter(fh);
	AKMesh::Point _v0 = _mesh.point(*fv_iter);
	fv_iter++;
	AKMesh::Point _v1 = _mesh.point(*fv_iter);
	fv_iter++;
	AKMesh::Point _v2 = _mesh.point(*fv_iter);
	{
		const AKMesh::Point v0v1 = _v1 - _v0;
		const AKMesh::Point v0v2 = _v2 - _v0;
		const AKMesh::Point n = v0v1 % v0v2; // not normalized !
		const AKNumber d = n.sqrnorm();
		// Check if the triangle is degenerated
		if (d < FLT_MIN && d > -FLT_MIN) {
			return -1.0;
		}
		const AKNumber invD = 1.0 / d;
		// these are not needed for every point, should still perform
		// better with many points against one triangle
		const AKMesh::Point v1v2 = _v2 - _v1;
		const AKNumber inv_v0v2_2 = 1.0 / v0v2.sqrnorm();
		const AKNumber inv_v0v1_2 = 1.0 / v0v1.sqrnorm();
		const AKNumber inv_v1v2_2 = 1.0 / v1v2.sqrnorm();
		AKMesh::Point v0p = _p - _v0;
		AKMesh::Point t = v0p % n;
		AKNumber  s01, s02, s12;
		const AKNumber a = (t | v0v2) * -invD;
		const AKNumber b = (t | v0v1) * invD;
		if (a < 0)
		{
			// Calculate the distance to an edge or a corner vertex
			s02 = (v0v2 | v0p) * inv_v0v2_2;
			if (s02 < 0.0)
			{
				s01 = (v0v1 | v0p) * inv_v0v1_2;
				if (s01 <= 0.0) {
					v0p = _v0;
				}
				else if (s01 >= 1.0) {
					v0p = _v1;
				}
				else {
					v0p = _v0 + v0v1 * s01;
				}
			}
			else if (s02 > 1.0) {
				s12 = (v1v2 | (_p - _v1)) * inv_v1v2_2;
				if (s12 >= 1.0) {
					v0p = _v2;
				}
				else if (s12 <= 0.0) {
					v0p = _v1;
				}
				else {
					v0p = _v1 + v1v2 * s12;
				}
			}
			else {
				v0p = _v0 + v0v2 * s02;
			}
		}
		else if (b < 0.0) {
			// Calculate the distance to an edge or a corner vertex
			s01 = (v0v1 | v0p) * inv_v0v1_2;
			if (s01 < 0.0)
			{
				const AKMesh::Point n = v0v1 % v0v2; // not normalized !
				s02 = (v0v2 | v0p) * inv_v0v2_2;
				if (s02 <= 0.0) {
					v0p = _v0;
				}
				else if (s02 >= 1.0) {
					v0p = _v2;
				}
				else {
					v0p = _v0 + v0v2 * s02;
				}
			}
			else if (s01 > 1.0) {
				s12 = (v1v2 | (_p - _v1)) * inv_v1v2_2;
				if (s12 >= 1.0) {
					v0p = _v2;
				}
				else if (s12 <= 0.0) {
					v0p = _v1;
				}
				else {
					v0p = _v1 + v1v2 * s12;
				}
			}
			else {
				v0p = _v0 + v0v1 * s01;
			}
		}
		else if (a + b > 1.0) {
			// Calculate the distance to an edge or a corner vertex
			s12 = (v1v2 | (_p - _v1)) * inv_v1v2_2;
			if (s12 >= 1.0) {
				s02 = (v0v2 | v0p) * inv_v0v2_2;
				if (s02 <= 0.0) {
					v0p = _v0;
				}
				else if (s02 >= 1.0) {
					v0p = _v2;
				}
				else {
					v0p = _v0 + v0v2*s02;
				}
			}
			else if (s12 <= 0.0) {
				s01 = (v0v1 | v0p) * inv_v0v1_2;
				if (s01 <= 0.0) {
					v0p = _v0;
				}
				else if (s01 >= 1.0) {
					v0p = _v1;
				}
				else {
					v0p = _v0 + v0v1 * s01;
				}
			}
			else {
				v0p = _v1 + v1v2 * s12;
			}
		}
		else {
			// Calculate the distance to an interior point of the triangle
			return ((_p - n*((n | v0p) * invD)) - _p).sqrnorm();
		}
		return (v0p - _p).sqrnorm();
	}
}

AKNumber point_triangle_minDist(AKMesh & _mesh, set<size_t> & _tarcomputeFaces, AKMesh::Point _point, AKMesh::FaceHandle & min_dist_fh)
{
	AKNumber min = 9e10;
	for (set<size_t>::iterator it = _tarcomputeFaces.begin(); it != _tarcomputeFaces.end(); it++)
	{
		AKMesh::FaceHandle fh((int)(*it));
		AKNumber d = point_triangle_dist(_mesh, fh, _point);
		if (d<min)
		{
			min = d;
			min_dist_fh = fh;
		}
	}
	return min;
}

AKNumber compute_mean_vertex_distance_error(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices, AKKDTree & dmKDTreeFaces, AKMesh & mesh)
{
	desired_mesh.request_face_normals();
	desired_mesh.update_normals();
	vector<AKNumber> mesh_vertex_areas;
	compute_all_vertex_areas(mesh, mesh_vertex_areas);

	AKNumber search_radius = average_edge_length_factor*compute_average_edge_length(desired_mesh);
	AKNumber mesh_total_area = 0.0;
	AKNumber error_sum = 0.0;
	AKNumber error_max = -9e10;
	for (AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		AKNumber query_pt[3] = {};

		query_pt[0] = mesh.point(*v_it)[0];
		query_pt[1] = mesh.point(*v_it)[1];
		query_pt[2] = mesh.point(*v_it)[2];
		vector<pair<size_t, AKNumber> >   ret_matches_v;
		vector<pair<size_t, AKNumber> >   ret_matches_f;
		nanoflann::SearchParams params;
		//params.sorted = false;

		dmKDTreeVertices.radiusSearch(&query_pt[0], search_radius, ret_matches_v, params);
		dmKDTreeFaces.radiusSearch(&query_pt[0], search_radius, ret_matches_f, params);

		set<size_t> tarcomputeFaces;
		for (int i = 0; i<ret_matches_f.size(); i++)
			tarcomputeFaces.insert(ret_matches_f[i].first);
		for (int i = 0; i<ret_matches_v.size(); i++)
		{
			AKMesh::VertexHandle vh(ret_matches_v[i].first);
			for (AKMesh::VertexFaceIter vf_it = desired_mesh.vf_iter(vh); vf_it.is_valid(); vf_it++)
			{
				tarcomputeFaces.insert(vf_it->idx());
			}
		}
		if (tarcomputeFaces.size()>0)
		{
			AKMesh::FaceHandle min_dist_fh;
			AKNumber min_dist = point_triangle_minDist(desired_mesh, tarcomputeFaces, mesh.point(*v_it), min_dist_fh);
			mesh_total_area += mesh_vertex_areas[v_it->idx()];
			error_sum += min_dist*mesh_vertex_areas[v_it->idx()];
			if (min_dist>error_max)
				error_max = min_dist;
		}
	}
	//cout<<endl<<error_sum/mesh_total_area<<" "<<error_max;
	return error_sum / mesh_total_area;
}

AKNumber compute_L2_vertex_based_error(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices, AKKDTree & dmKDTreeFaces, AKMesh & mesh)
{
	vector<AKNumber> mesh_vertex_areas;
	compute_all_vertex_areas(mesh, mesh_vertex_areas);

	AKNumber search_radius = average_edge_length_factor*compute_average_edge_length(desired_mesh);
	AKNumber mesh_total_area = 0.0;
	AKNumber error_sum = 0.0;
	for (AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		AKNumber query_pt[3] = {};

		query_pt[0] = mesh.point(*v_it)[0];
		query_pt[1] = mesh.point(*v_it)[1];
		query_pt[2] = mesh.point(*v_it)[2];
		vector<pair<size_t, AKNumber> >   ret_matches_v;
		vector<pair<size_t, AKNumber> >   ret_matches_f;
		nanoflann::SearchParams params;

		dmKDTreeVertices.radiusSearch(&query_pt[0], search_radius, ret_matches_v, params);
		dmKDTreeFaces.radiusSearch(&query_pt[0], search_radius, ret_matches_f, params);

		set<size_t> tarcomputeFaces;
		for (int i = 0; i<ret_matches_f.size(); i++)
			tarcomputeFaces.insert(ret_matches_f[i].first);
		for (int i = 0; i<ret_matches_v.size(); i++)
		{
			AKMesh::VertexHandle vh((int)(ret_matches_v[i].first));
			for (AKMesh::VertexFaceIter vf_it = desired_mesh.vf_iter(vh); vf_it.is_valid(); vf_it++)
			{
				tarcomputeFaces.insert(vf_it->idx());
			}
		}
		if (tarcomputeFaces.size()>0)
		{
			AKMesh::FaceHandle min_dist_fh;
			AKNumber min_dist = point_triangle_minDist(desired_mesh, tarcomputeFaces, mesh.point(*v_it), min_dist_fh);
			mesh_total_area += mesh_vertex_areas[v_it->idx()];
			error_sum += pow(min_dist, 2.0)*mesh_vertex_areas[v_it->idx()];
		}
	}
	return sqrt(error_sum / mesh_total_area);

}

AKNumber compute_mean_quadric_error(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices, AKMesh & mesh)
{
	desired_mesh.request_face_normals();
	desired_mesh.request_vertex_normals();
	desired_mesh.update_normals();
	vector<AKNumber> desired_mesh_face_areas;
	vector<AKNumber> mesh_vertex_areas;
	compute_all_face_areas(desired_mesh, desired_mesh_face_areas);
	compute_all_vertex_areas(mesh, mesh_vertex_areas);
	AKNumber error_sum = 0.0;
	AKNumber total_area = 0.0;
	for (AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		AKNumber query_pt[3] = {};

		query_pt[0] = mesh.point(*v_it)[0];
		query_pt[1] = mesh.point(*v_it)[1];
		query_pt[2] = mesh.point(*v_it)[2];

		const size_t num_results = 1;
		size_t ret_index;
		AKNumber out_dist_sqr;
		nanoflann::KNNResultSet<AKNumber> resultSet(num_results);
		resultSet.init(&ret_index, &out_dist_sqr);
		dmKDTreeVertices.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
		AKNumber vertex_quadric = 0.0;
		AKMesh::VertexHandle nearest_vertex((int)ret_index);
		for (AKMesh::VertexFaceIter vf_it = desired_mesh.vf_iter(nearest_vertex); vf_it.is_valid(); vf_it++)
		{
			AKMesh::Point p = mesh.point(*v_it);
			AKMesh::FaceHandle fh = *vf_it;
			AKMesh::Normal face_normal = desired_mesh.normal(fh);
			face_normal.normalize_cond();

			AKMesh::FaceVertexIter fv_iter = desired_mesh.fv_iter(fh);
			AKMesh::Point p0 = desired_mesh.point(*fv_iter);
			AKNumber d = -(face_normal | p0);

			AKNumber face_quadric = pow((face_normal | p) + d, 2.0);
			vertex_quadric += face_quadric*desired_mesh_face_areas[vf_it->idx()];
		}
		error_sum += mesh_vertex_areas[v_it->idx()] * vertex_quadric;
		total_area += mesh_vertex_areas[v_it->idx()];
	}
	return error_sum / total_area;
}

AKNumber compute_mean_square_angle_error(AKMesh &desired_mesh, AKMesh & mesh)
{
	desired_mesh.request_face_normals();
	desired_mesh.update_face_normals();

	mesh.request_face_normals();
	mesh.update_face_normals();
	AKNumber mean_square_angle_error = 0.0;
	for (AKMesh::FaceIter f_it1 = mesh.faces_begin(), f_it2 = desired_mesh.faces_begin();
		f_it1 != mesh.faces_end(); f_it1++, f_it2++)
	{
		AKMesh::Normal normal1 = mesh.normal(*f_it1);
		AKMesh::Normal normal2 = desired_mesh.normal(*f_it2);
		AKNumber dot_product = normal1 | normal2;
		dot_product = min(1.0, max(dot_product, -1.0));
		mean_square_angle_error += acos(dot_product) * 180.0 / M_PI;
	}
	return mean_square_angle_error / (AKNumber)desired_mesh.n_faces();
}

AKNumber compute_mean_square_angle_error_2(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices, AKKDTree & dmKDTreeFaces, AKMesh & mesh)
{
	desired_mesh.request_face_normals();
	desired_mesh.update_normals();
	mesh.request_face_normals();
	mesh.update_normals();
	vector<AKMesh::Point> centroids;
	compute_all_face_centroids(mesh, centroids);

	AKNumber search_radius = average_edge_length_factor*compute_average_edge_length(desired_mesh);
	AKNumber num_vert = 0.0;
	AKNumber error_sum = 0.0;
	for (size_t i = 0; i<centroids.size(); i++)
	{
		AKNumber query_pt[3] = {};

		query_pt[0] = centroids[i][0];
		query_pt[1] = centroids[i][1];
		query_pt[2] = centroids[i][2];
		vector<pair<size_t, AKNumber> >   ret_matches_v;
		vector<pair<size_t, AKNumber> >   ret_matches_f;
		nanoflann::SearchParams params;

		dmKDTreeVertices.radiusSearch(&query_pt[0], search_radius, ret_matches_v, params);
		dmKDTreeFaces.radiusSearch(&query_pt[0], search_radius, ret_matches_f, params);

		set<size_t> tarcomputeFaces;
		for (int j = 0; j<ret_matches_f.size(); j++)
			tarcomputeFaces.insert(ret_matches_f[j].first);
		for (int j = 0; j<ret_matches_v.size(); j++)
		{
			AKMesh::VertexHandle vh((int)(ret_matches_v[j].first));
			for (AKMesh::VertexFaceIter vf_it = desired_mesh.vf_iter(vh); vf_it.is_valid(); vf_it++)
			{
				tarcomputeFaces.insert(vf_it->idx());
			}
		}
		if (tarcomputeFaces.size()>0)
		{
			AKMesh::FaceHandle min_dist_fh;
			point_triangle_minDist(desired_mesh, tarcomputeFaces, centroids[i], min_dist_fh);
			num_vert += 1.0;
			AKMesh::Normal normal1 = mesh.normal(AKMesh::FaceHandle((int)i));
			AKMesh::Normal normal2 = desired_mesh.normal(min_dist_fh);
			AKNumber dot_product = normal1 | normal2;
			dot_product = min(1.0, max(dot_product, -1.0));
			error_sum += acos(dot_product) * 180.0 / M_PI;
		}
	}
	return error_sum / num_vert;
}

AKNumber compute_L2_normal_based_error(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices, AKKDTree & dmKDTreeFaces, AKMesh & mesh)
{
	desired_mesh.request_face_normals();
	desired_mesh.update_normals();
	mesh.request_face_normals();
	mesh.update_face_normals();
	vector<AKNumber> mesh_face_areas;
	compute_all_face_areas(mesh, mesh_face_areas);
	vector<AKMesh::Point> centroids;
	compute_all_face_centroids(mesh, centroids);

	AKNumber search_radius = average_edge_length_factor*compute_average_edge_length(desired_mesh);
	AKNumber mesh_total_area = 0.0;
	AKNumber error_sum = 0.0;
	for (size_t i = 0; i<centroids.size(); i++)
	{
		AKNumber query_pt[3] = {};

		query_pt[0] = centroids[i][0];
		query_pt[1] = centroids[i][1];
		query_pt[2] = centroids[i][2];
		vector<pair<size_t, AKNumber> >   ret_matches_v;
		vector<pair<size_t, AKNumber> >   ret_matches_f;
		nanoflann::SearchParams params;

		dmKDTreeVertices.radiusSearch(&query_pt[0], search_radius, ret_matches_v, params);
		dmKDTreeFaces.radiusSearch(&query_pt[0], search_radius, ret_matches_f, params);

		set<size_t> tarcomputeFaces;
		for (int j = 0; j<ret_matches_f.size(); j++)
			tarcomputeFaces.insert(ret_matches_f[j].first);
		for (int j = 0; j<ret_matches_v.size(); j++)
		{
			AKMesh::VertexHandle vh((int)(ret_matches_v[j].first));
			for (AKMesh::VertexFaceIter vf_it = desired_mesh.vf_iter(vh); vf_it.is_valid(); vf_it++)
			{
				tarcomputeFaces.insert(vf_it->idx());
			}
		}
		if (tarcomputeFaces.size()>0)
		{
			AKMesh::FaceHandle min_dist_fh;
			point_triangle_minDist(desired_mesh, tarcomputeFaces, centroids[i], min_dist_fh);
			mesh_total_area += mesh_face_areas[i];
			AKMesh::Normal normal1 = mesh.normal(AKMesh::FaceHandle((int)i));
			AKMesh::Normal normal2 = desired_mesh.normal(min_dist_fh);
			AKNumber temp = pow((normal1 - normal2).norm(), 2.0);
			error_sum += temp*mesh_face_areas[i];
			/*if (temp > 5.0 )
			{
				cout << normal1 << endl;
				cout << normal2 << endl;
				cout << temp << endl;
				cout << min_dist_fh.idx() << endl;
				cout << mesh.n_faces() << " " << desired_mesh.n_faces() << endl;
				cout << i << " " << error_sum << endl;
				cout << i << " " << error_sum << endl;
			}*/

		}
	}
	return error_sum / mesh_total_area;
}

AKNumber compute_mean_tangential_error(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices, AKMesh & mesh)
{
	desired_mesh.request_face_normals();
	desired_mesh.request_vertex_normals();
	desired_mesh.update_normals();
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();
	vector<AKNumber> mesh_vertex_areas;
	compute_all_vertex_areas(mesh, mesh_vertex_areas);
	AKNumber error_sum = 0.0;
	AKNumber mesh_total_area = 0.0;
	for (AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		AKNumber query_pt[3] = {};

		query_pt[0] = mesh.point(*v_it)[0];
		query_pt[1] = mesh.point(*v_it)[1];
		query_pt[2] = mesh.point(*v_it)[2];

		const size_t num_results = 1;
		size_t ret_index;
		AKNumber out_dist_sqr;
		nanoflann::KNNResultSet<AKNumber> resultSet(num_results);
		resultSet.init(&ret_index, &out_dist_sqr);
		dmKDTreeVertices.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));


		AKMesh::Normal normal1 = mesh.normal(*v_it);
		AKMesh::Normal normal2 = desired_mesh.normal(AKMesh::VertexHandle((int)ret_index));
		AKNumber temp = pow((normal1 - normal2).norm(), 2.0);

		mesh_total_area += mesh_vertex_areas[v_it->idx()];
		error_sum += temp*mesh_vertex_areas[v_it->idx()];
	}
	return error_sum / mesh_total_area;
}

AKNumber compute_mean_discrete_curvature_error(AKMesh & desired_mesh, AKKDTree &dmKDTreeVertices, vector<AKNumber> & dmCurvature, AKMesh & mesh, vector<AKNumber> & mCurvature)
{
	desired_mesh.request_face_normals();
	desired_mesh.request_vertex_normals();
	desired_mesh.update_normals();
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();
	vector<AKNumber> mesh_vertex_areas;
	compute_all_vertex_areas(mesh, mesh_vertex_areas);
	AKNumber error_sum = 0.0;
	AKNumber mesh_total_area = 0.0;
	for (AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		AKNumber query_pt[3] = {};

		query_pt[0] = mesh.point(*v_it)[0];
		query_pt[1] = mesh.point(*v_it)[1];
		query_pt[2] = mesh.point(*v_it)[2];

		const size_t num_results = 1;
		size_t ret_index;
		AKNumber out_dist_sqr;
		nanoflann::KNNResultSet<AKNumber> resultSet(num_results);
		resultSet.init(&ret_index, &out_dist_sqr);
		dmKDTreeVertices.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

		AKNumber temp = abs(dmCurvature[ret_index] - mCurvature[v_it->idx()]);

		mesh_total_area += mesh_vertex_areas[v_it->idx()];
		error_sum += temp*mesh_vertex_areas[v_it->idx()];
	}
	return error_sum / mesh_total_area;
}

