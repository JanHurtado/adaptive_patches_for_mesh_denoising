#include "gui/myDataManager.h"


myDataManager::myDataManager()
{
}


myDataManager::~myDataManager()
{
	input_mesh.clear();
	output_mesh.clear();
	input_mesh_shape.clear();
	output_mesh_shape.clear();
	selection.clear();
}


void myDataManager::loadInputMesh(string & fileName)
{
	input_mesh.clear();
	output_mesh.clear();
	input_mesh_shape.clear();
	output_mesh_shape.clear();
	face_scalar_functions_sparse.clear();
	face_scalar_functions_sparse_2.clear();

	import_mesh(input_mesh, fileName);
	normalize_mesh(input_mesh, scale_factor);

	input_mesh.request_face_normals();
	input_mesh.request_vertex_normals();
	input_mesh.update_normals();

	face_scalar_min = 0.0f;
	face_scalar_max = 1.0f;
	face_scalar_default = 1.0f;
	face_scalar_functions_sparse.resize(input_mesh.n_faces());
	face_scalar_functions_sparse_2.resize(input_mesh.n_faces());

	//normalize_mesh(input_mesh);
	output_mesh = input_mesh;


	
}

void myDataManager::saveOutputMesh(string & fileName)
{
	retrieve_mesh(output_mesh, scale_factor);
	output_mesh.request_face_normals();
	output_mesh.request_vertex_normals();
	output_mesh.update_normals();
	export_mesh(output_mesh, fileName);
}

void myDataManager::updateInputShape()
{
	input_mesh_shape.clear();
	//input_mesh_shape.loadMesh(input_mesh);
	//input_mesh_shape.loadMeshFaceScalar(input_mesh, face_scalar_functions[0], face_scalar_min, face_scalar_max);
	input_mesh_shape.loadMeshFaceScalarSparse(input_mesh, face_scalar_functions_sparse[0], face_scalar_min, face_scalar_max, face_scalar_default);

}

void myDataManager::updateOutputShape()
{
	output_mesh_shape.clear();
	//output_mesh_shape.loadMesh(output_mesh);
	//output_mesh_shape.loadMeshFaceScalar(output_mesh, face_scalar_functions[0], face_scalar_min, face_scalar_max);
	output_mesh_shape.loadMeshFaceScalarSparse(output_mesh, face_scalar_functions_sparse[0], face_scalar_min, face_scalar_max, face_scalar_default);
}

void myDataManager::updateShapes()
{
	updateInputShape();
	updateOutputShape();
}

void myDataManager::updateOutputSelection(vector<size_t> & indices)
{
	selection.clear();
	selection.loadMeshVertexSelection(output_mesh, indices);
}

void myDataManager::setOutputAsInput()
{
	input_mesh = output_mesh;
	updateShapes();
}

void myDataManager::reinitialize()
{
	output_mesh = input_mesh;
	updateShapes();
}

void myDataManager::setFaceScalarFunction(int idx)
{
	//output_mesh_shape.updateFaceScalar(face_scalar_functions[idx], face_scalar_min, face_scalar_max);
	//input_mesh_shape.updateFaceScalar(face_scalar_functions[idx], face_scalar_min, face_scalar_max);
	output_mesh_shape.updateFaceScalarSparse(face_scalar_functions_sparse_2[idx], face_scalar_min, face_scalar_max,face_scalar_default);
	input_mesh_shape.updateFaceScalarSparse(face_scalar_functions_sparse_2[idx], face_scalar_min, face_scalar_max,face_scalar_default);
}

void myDataManager::setSinglePatchFunction()
{
	output_mesh_shape.updateFaceScalarSparse(sp_patch_function, face_scalar_min, face_scalar_max, face_scalar_default);
	input_mesh_shape.updateFaceScalarSparse(sp_patch_function, face_scalar_min, face_scalar_max, face_scalar_default);
}

void myDataManager::genGeodesicData(std::vector<double> & points, std::vector<unsigned> & faces)
{
	points.resize(output_mesh.n_vertices() * 3);
	faces.resize(output_mesh.n_faces() * 3);

	int v_pos = 0;
	for (AKMesh::VertexIter v_it = output_mesh.vertices_begin(); v_it != output_mesh.vertices_end(); v_it++)
	{
		AKMesh::Point current_point = output_mesh.point(*v_it);
		points[v_pos] = current_point[0];
		v_pos++;
		points[v_pos] = current_point[1];
		v_pos++;
		points[v_pos] = current_point[2];
		v_pos++;
	}
	int f_pos = 0;
	for (AKMesh::FaceIter f_it = output_mesh.faces_begin(); f_it != output_mesh.faces_end(); f_it++)
	{
		for (AKMesh::FaceVertexIter fv_it = output_mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++)
		{
			faces[f_pos] = fv_it->idx();
			f_pos++;
		}
	}
}