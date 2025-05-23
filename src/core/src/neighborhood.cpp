#include <core/neighborhood.h>

/////////////////////////////////////////////////
/// Vertex neighborhood
/////////////////////////////////////////////////

void compute_vertex_neighbors(AKMesh &mesh, AKMesh::VertexHandle vh, int k,
                       vector<AKMesh::VertexHandle> &vertex_neighbors)
{
    vertex_neighbors.clear();
    vector<bool> mark(mesh.n_vertices(),false);
    queue<AKMesh::VertexHandle> queue_vertex_handle;
    queue<int> queue_depth_handle;
    mark[vh.idx()] = true;
    queue_vertex_handle.push(vh);
    queue_depth_handle.push(0);
    AKMesh::Point ci = mesh.point(vh);
    while(!queue_vertex_handle.empty())
    {
        AKMesh::VertexHandle vh = queue_vertex_handle.front();
        int current_depth = queue_depth_handle.front();
        vertex_neighbors.push_back(vh);
        queue_vertex_handle.pop();
        queue_depth_handle.pop();
        for(AKMesh::VertexVertexIter vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it)
        {
            AKMesh::VertexHandle vh_neighbor = *vv_it;
            if(mark[vh_neighbor.idx()] == false)
            {
                int new_depth = current_depth + 1;
                if(new_depth <= k)
                {
                    queue_vertex_handle.push(vh_neighbor);
                    queue_depth_handle.push(new_depth);
                }
                mark[vh_neighbor.idx()] = true;
            }
        }
    }
}

void compute_adaptive_vertex_neighbors(AKMesh &mesh, AKMesh::VertexHandle vh, AKNumber radius,
                                                       vector<AKMesh::VertexHandle> &vertex_neighbors)
{
    vertex_neighbors.clear();
    vector<bool> mark(mesh.n_vertices(),false);
    queue<AKMesh::VertexHandle> queue_vertex_handle;
    mark[vh.idx()] = true;
    queue_vertex_handle.push(vh);
    AKMesh::Point ci = mesh.point(vh);

    while(!queue_vertex_handle.empty())
    {
        AKMesh::VertexHandle vh = queue_vertex_handle.front();
        vertex_neighbors.push_back(vh);
        queue_vertex_handle.pop();
        for(AKMesh::VertexVertexIter vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it)
        {
            AKMesh::VertexHandle vh_neighbor = *vv_it;
            if(mark[vh_neighbor.idx()] == false)
            {
                AKMesh::Point cj = mesh.point(vh_neighbor);
                AKNumber length = (cj - ci).length();
                if(length <= radius)
                    queue_vertex_handle.push(vh_neighbor);
                mark[vh_neighbor.idx()] = true;
            }
        }
    }
}

void compute_all_vertex_neighbors(AKMesh &mesh, int k,vector<vector<AKMesh::VertexHandle> > &all_vertex_neighbors)
{
    all_vertex_neighbors.resize(mesh.n_vertices());
    for(AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        vector<AKMesh::VertexHandle> vertex_neighbor;
        compute_vertex_neighbors(mesh, *v_it, k, vertex_neighbor);
        all_vertex_neighbors[v_it->idx()] = vertex_neighbor;
    }
}

void compute_all_adaptive_vertex_neighbors(AKMesh &mesh, AKNumber radius,
                               vector<vector<AKMesh::VertexHandle> > &all_vertex_neighbors)
{
    all_vertex_neighbors.resize(mesh.n_vertices());
    for(AKMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
    {
        vector<AKMesh::VertexHandle> vertex_neighbor;
        compute_adaptive_vertex_neighbors(mesh, *v_it, radius, vertex_neighbor);
        all_vertex_neighbors[v_it->idx()] = vertex_neighbor;
    }
}

/////////////////////////////////////////////////
/// Face neighborhood
/////////////////////////////////////////////////

void compute_edge_based_face_neighbors(AKMesh &mesh, AKMesh::FaceHandle fh, vector<AKMesh::FaceHandle> &face_neighbors)
{
    face_neighbors.clear();
    for(AKMesh::FaceFaceIter ff_it = mesh.ff_iter(fh); ff_it.is_valid(); ff_it++)
        face_neighbors.push_back(*ff_it);
}

void compute_vertex_based_face_neighbors(AKMesh &mesh, AKMesh::FaceHandle fh, vector<AKMesh::FaceHandle> &face_neighbors)
{
    face_neighbors.clear();
    set<int> neighbor_face_index; neighbor_face_index.clear();
    for(AKMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it.is_valid(); fv_it++)
    {
        for(AKMesh::VertexFaceIter vf_it = mesh.vf_iter(*fv_it); vf_it.is_valid(); vf_it++)
        {
            if((*vf_it) != fh)
                neighbor_face_index.insert(vf_it->idx());
        }
    }
    for(set<int>::iterator iter = neighbor_face_index.begin(); iter != neighbor_face_index.end(); ++ iter)
    {
        face_neighbors.push_back(AKMesh::FaceHandle(*iter));
    }
}

void compute_radius_based_face_neighbors(AKMesh &mesh, AKMesh::FaceHandle fh, AKNumber radius, vector<AKMesh::FaceHandle> &face_neighbors)
{
    AKMesh::Point ci = mesh.calc_face_centroid(fh);
    vector<bool> flag((int)mesh.n_faces(), false);

    face_neighbors.clear();
    flag[fh.idx()] = true;
    queue<AKMesh::FaceHandle> queue_face_handle;
    queue_face_handle.push(fh);

    vector<AKMesh::FaceHandle> temp_face_neighbors;
    while(!queue_face_handle.empty())
    {
        AKMesh::FaceHandle temp_face_handle_queue = queue_face_handle.front();
        if(temp_face_handle_queue != fh)
            face_neighbors.push_back(temp_face_handle_queue);
        queue_face_handle.pop();
        compute_vertex_based_face_neighbors(mesh, temp_face_handle_queue, temp_face_neighbors);
        for(int i = 0; i < (int)temp_face_neighbors.size(); i++)
        {
            AKMesh::FaceHandle temp_face_handle = temp_face_neighbors[i];
            if(!flag[temp_face_handle.idx()])
            {
                AKMesh::Point cj = mesh.calc_face_centroid(temp_face_handle);
                AKNumber distance = (ci - cj).length();
                if(distance <= radius)
                    queue_face_handle.push(temp_face_handle);
                flag[temp_face_handle.idx()] = true;
            }
        }
    }
}

void compute_all_edge_based_face_neighbors(AKMesh &mesh, bool include_tarcompute_face, vector<vector<AKMesh::FaceHandle> > &all_face_neighbors)
{
    all_face_neighbors.resize(mesh.n_faces());
    for(AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        vector<AKMesh::FaceHandle> face_neighbors;
        compute_edge_based_face_neighbors(mesh, *f_it, face_neighbors);
        if(include_tarcompute_face) face_neighbors.push_back(*f_it); //include tarcompute face
        all_face_neighbors[f_it->idx()] = face_neighbors;
    }
}

void compute_all_vertex_based_face_neighbors(AKMesh &mesh, bool include_tarcompute_face, vector<vector<AKMesh::FaceHandle> > &all_face_neighbors)
{
    all_face_neighbors.resize(mesh.n_faces());
    for(AKMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
    {
        vector<AKMesh::FaceHandle> face_neighbors;
        compute_vertex_based_face_neighbors(mesh, *f_it, face_neighbors);
        if(include_tarcompute_face) face_neighbors.push_back(*f_it); //include tarcompute face
        all_face_neighbors[f_it->idx()] = face_neighbors;
    }
}