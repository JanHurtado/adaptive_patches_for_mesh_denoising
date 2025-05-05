#include "core/util.h"

AKNumber compute_area (AKMesh & mesh)
{
    vector<AKNumber> areas;
    compute_all_face_areas(mesh,areas);
    AKNumber sum = 0;
    for (size_t i=0;i<areas.size();i++)
        sum+=areas[i];
    return sum;
}

AKNumber compute_volume ( AKMesh & mesh ) {

    typedef AKMesh::HalfedgeHandle hh_t;
    typedef AKMesh::VertexHandle vh_t;
    typedef AKMesh::Point p_t;

    AKNumber volume = 0.0f;
    AKNumber three_inv = 1.0f / 3.0f;
    AKMesh::FaceIter f_it;
    for (f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
    {
        hh_t hh = mesh.halfedge_handle(*f_it);
        hh_t hhn = mesh.next_halfedge_handle(hh);

        vh_t v0(mesh.from_vertex_handle(hh));
        vh_t v1(mesh.to_vertex_handle(hh));
        vh_t v2(mesh.to_vertex_handle(hhn));

        p_t p0(mesh.point(v0));
        p_t p1(mesh.point(v1));
        p_t p2(mesh.point(v2));

        p_t g = (p0 + p1 + p2) * three_inv;
        p_t n = (p1 - p0) % (p2 - p0);
        volume += (g|n);
    }
    volume *= 1.0f / 6.0f;
    return volume;
}

AKNumber compute_average_edge_length(AKMesh &mesh)
{
    AKNumber average_edge_length = 0.0;
    for(AKMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
        average_edge_length += mesh.calc_edge_length(*e_it);
    AKNumber edgeNum = (AKNumber)mesh.n_edges();
	if (average_edge_length != 0)
		average_edge_length /= edgeNum;

    return average_edge_length;
}

void compute_all_face_areas(AKMesh &mesh, vector<AKNumber> &areas)
{
    areas.resize(mesh.n_faces());

    for(AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        vector<AKMesh::Point> point;
        point.resize(3); int index = 0;
        for(AKMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++)
        {
            point[index] = mesh.point(*fv_it);
            index++;
        }
        AKMesh::Point edge1 = point[1] - point[0];
        AKMesh::Point edge2 = point[1] - point[2];
        AKNumber S = 0.5f * (edge1 % edge2).length();
        areas[(*f_it).idx()] = S;
    }
}

void compute_all_face_centroids(AKMesh &mesh, vector<AKMesh::Point> &centroids)
{
    centroids.resize(mesh.n_faces(), AKMesh::Point(0.0, 0.0, 0.0));
    for(AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        AKMesh::Point pt = mesh.calc_face_centroid(*f_it);
        centroids[(*f_it).idx()] = pt;
    }
}

void compute_all_face_normals(AKMesh &mesh, vector<AKMesh::Normal> &normals)
{
    mesh.request_face_normals();
    mesh.update_face_normals();

    normals.resize(mesh.n_faces());
    for(AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        AKMesh::Normal n = mesh.normal(*f_it);
        normals[f_it->idx()] = n;
    }
}

void computeAllFaceNormals2(AKMesh &mesh, vector<AKMesh::Normal> &normals)
{
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();

	normals.resize(mesh.n_faces());
	for (AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		AKMesh::Normal fn = mesh.normal(*f_it);
		AKMesh::Normal vn = AKMesh::Normal(0, 0, 0);
		for (AKMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++)
		{
			vn += mesh.normal(*fv_it);
		}
		vn.normalize_cond();
		normals[f_it->idx()] = (fn + vn).normalize_cond();
	}
}

void compute_face_vertex_angle(AKMesh &mesh,AKMesh::FaceHandle fh,AKMesh::VertexHandle vh, AKNumber & angle)
{
    AKMesh::FaceVertexIter verts = mesh.fv_iter(fh);
    AKMesh::VertexHandle vh1,vh2;
    if (vh == *verts)
    {
        vh1 = *(++verts);
        vh2 = *(++verts);
    }
    else
    {
        vh1 = *verts;
        vh2 = (vh == *(++verts)) ? *(++verts) : *verts;
    }
    AKMesh::Point v1 = mesh.point(vh) - mesh.point(vh1);
    AKMesh::Point v2 = mesh.point(vh) - mesh.point(vh2);
    AKMesh::Point v3 = (mesh.point(vh1) - mesh.point(vh2));
    v1.normalize_cond(); v2.normalize_cond();
    angle = acos(v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}

AKNumber compute_vertex_area(AKMesh & mesh,AKMesh::VertexHandle vh,vector<AKNumber> & areas)
{
    AKNumber area_sum = 0.0f;
    for(AKMesh::VertexFaceIter vf_iter = mesh.vf_iter(vh);vf_iter.is_valid();vf_iter++)
        area_sum+=areas[vf_iter->idx()];
    return area_sum*(1.0f/3.0f);
}

void compute_all_vertex_areas(AKMesh & mesh,vector<AKNumber> & areas)
{
    areas.resize(mesh.n_vertices());
    vector<AKNumber> face_areas;
    compute_all_face_areas(mesh,face_areas);
    for(AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
        areas[v_it->idx()]=compute_vertex_area(mesh,*v_it,face_areas);
}

void compute_all_points(AKMesh &mesh, vector<AKMesh::Point> & points)
{
    points.resize(mesh.n_vertices());
    for(AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        AKMesh::Point p = mesh.point(*v_it);
        points[v_it->idx()] = p;
    }
}

void generate_data_vectors(AKMesh & _mesh, std::vector<AKNumber> & _points, std::vector<unsigned> & _faces)
{
	_points.resize(_mesh.n_vertices() * 3);
	_faces.resize(_mesh.n_faces() * 3);

	int v_pos = 0;
	for (AKMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); v_it++)
	{
		AKMesh::Point current_point = _mesh.point(*v_it);
		_points[v_pos] = current_point[0];
		v_pos++;
		_points[v_pos] = current_point[1];
		v_pos++;
		_points[v_pos] = current_point[2];
		v_pos++;
	}
	int f_pos = 0;
	for (AKMesh::FaceIter f_it = _mesh.faces_begin(); f_it != _mesh.faces_end(); f_it++)
	{
		for (AKMesh::FaceVertexIter fv_it = _mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++)
		{
			_faces[f_pos] = fv_it->idx();
			f_pos++;
		}
	}
}


void normalize_mesh(AKMesh & _mesh)
{
	AKMesh::Point centroid(0.0f, 0.0f, 0.0f);
	for (AKMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); v_it++)
	{
		centroid += _mesh.point(*v_it);
	}
	centroid = centroid / (static_cast<AKNumber>(_mesh.n_vertices()));
	AKNumber max_dist = 0.0f;
	for (AKMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); v_it++)
	{
		AKNumber temp = (centroid - _mesh.point(*v_it)).length();
		if (temp > max_dist)
			max_dist = temp;
	}
	for (AKMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); v_it++)
		_mesh.set_point(*v_it, (_mesh.point(*v_it) - centroid) / max_dist);
}

void normalize_mesh(AKMesh & _mesh, AKNumber & scale_factor)
{
	AKNumber avg_edge_length = compute_average_edge_length(_mesh);
	scale_factor = 1.0 / avg_edge_length;
	for (AKMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); v_it++)
		_mesh.set_point(*v_it, (_mesh.point(*v_it) * scale_factor));
}

void retrieve_mesh(AKMesh & _mesh, AKNumber scale_factor)
{
	AKNumber temp = 1.0 / scale_factor;
	for (AKMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); v_it++)
		_mesh.set_point(*v_it, (_mesh.point(*v_it) * temp));
}