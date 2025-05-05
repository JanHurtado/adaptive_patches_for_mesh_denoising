#ifndef MYSHAPE_H
#define MYSHAPE_H

#include <glm/glm.hpp>
#include <GL/glew.h>
#include <core/iomesh.h>






/** @addtogroup visualization
  * @brief Shape data.
  *
  * @{
  */


/**
 * @brief The Vertex struct - single vertex with its corresponding attributes
 */
struct Vertex
{
    /**
     * @brief position - vertex position (x,y,z)
     */
    glm::vec3 position;
    /**
     * @brief color - vertex color (r,g,b)
     */
	glm::vec3 color;
    /**
     * @brief normal - vertex normal (nx,ny,nz)
     */
	glm::vec3 normal;
};

/**
 * @brief The ShapeData struct - triangular mesh for OpenGL buffers manipulation
 */
struct ShapeData
{
    /**
     * @brief numVertices - number of vertices of the mesh
     */
    GLuint numVertices;
    /**
     * @brief numIndices - number of indices of all triangles
     */
    GLuint numIndices;

    /**
     * @brief vertices - array containing all vertices
     */
    Vertex* vertices;
    /**
     * @brief indices - array containing all indices
     */
    GLuint* indices;

    /**
     * @brief centroid - shape centroid
     */
	glm::vec3 centroid;

    /**
     * @brief ShapeData - default constructor
     */
    ShapeData():vertices(0), numVertices(0),indices(0), numIndices(0){}
    /**
     * @brief ShapeData - generates shape data using a triangular mesh file supported in OpenMesh library
     * @param fileName - mesh file name
     */
	ShapeData(string fileName)
	{
		loadFromFile(fileName);
	}
    /**
     * @brief loadFromFile - load shape data using a triangular mesh file supported in OpenMesh library
     * @param fileName - mesh file name
     */
	void loadFromFile(string fileName)
	{
		AKMesh _mesh;
		import_mesh(_mesh, fileName);
		_mesh.request_face_normals();
		_mesh.request_vertex_normals();
		_mesh.update_normals();
		loadMesh(_mesh);
	}
    /**
     * @brief loadMesh - load shape data from a data structure defined in OpenMesh library
     * @param _mesh - triangular mesh
     */
	void loadMesh(AKMesh & _mesh)
	{
		numVertices = static_cast<GLuint>(_mesh.n_vertices());
		vertices = new Vertex[numVertices];
		numIndices = static_cast<GLuint>(_mesh.n_faces() * 3);
		indices = new GLuint[numIndices];

		GLuint currentIndex = 0;
		AKMesh::Point c(0, 0, 0);
		float n = 0.0f;
		for (AKMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); v_it++)
		{
			AKMesh::Color current_color = _mesh.color(*v_it);
			AKMesh::Point current_point = _mesh.point(*v_it);
			AKMesh::Normal current_normal = _mesh.normal(*v_it);
			c += current_point;
			n++;
			vertices[currentIndex].position = glm::vec3(current_point[0], current_point[1], current_point[2]);
			vertices[currentIndex].color = glm::vec3(0.5f, 0.5f, 0.5f);
			vertices[currentIndex].normal = glm::vec3(current_normal[0], current_normal[1], current_normal[2]);
			currentIndex++;
		}
		c = c / n;
		centroid[0] = c[0];
		centroid[1] = c[1];
		centroid[2] = c[2];
		currentIndex = 0;
		for (AKMesh::FaceIter f_it = _mesh.faces_begin(); f_it != _mesh.faces_end(); f_it++)
		{
			for (AKMesh::FaceVertexIter fv_it = _mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++)
			{
				indices[currentIndex] = static_cast<GLuint>(fv_it->idx());
				currentIndex++;
			}
		}
	}
    /**
     * @brief loadMeshVertexSelection - load as shape data a subset of vertices (only vertices) from a data structure defined in OpenMesh library
     * @param _mesh - triangular mesh
     * @param selected_vertices - subset of vertices ids of _mesh
     */
	void loadMeshVertexSelection(AKMesh & _mesh,vector<size_t> & selected_vertices)
	{
		numVertices = static_cast<GLuint>(_mesh.n_vertices());
		vertices = new Vertex[numVertices];
		numIndices = 0;
		indices = new GLuint[numIndices];
		AKMesh::Point c(0, 0, 0);
		float n = 0.0f;
		for (size_t i = 0; i < selected_vertices.size(); i++)
		{
			AKMesh::VertexHandle vh((int)selected_vertices[i]);
			AKMesh::Point current_point = _mesh.point(vh);
			AKMesh::Normal current_normal = _mesh.normal(vh);
			c += current_point;
			n++;
			vertices[i].position = glm::vec3(current_point[0], current_point[1], current_point[2]);
			vertices[i].color = glm::vec3(1.0f, 0.0f, 0.0f);
			vertices[i].normal = glm::vec3(current_normal[0], current_normal[1], current_normal[2]);
		}
		c = c / n;
		centroid[0] = c[0];
		centroid[1] = c[1];
		centroid[2] = c[2];
	}
	/**
	*/
	void loadMeshFaceColor(AKMesh & _mesh, vector<glm::vec3> & face_color)
	{
		numVertices = static_cast<GLuint>(_mesh.n_faces() * 3);
		vertices = new Vertex[numVertices];
		numIndices = static_cast<GLuint>(_mesh.n_faces() * 3);
		indices = new GLuint[numIndices];

		GLuint currentIndex = 0;
		AKMesh::Point c(0, 0, 0);
		float n = 0.0f;
		for (AKMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); v_it++)
		{
			AKMesh::Point current_point = _mesh.point(*v_it);
			c += current_point;
			n++;
			currentIndex++;
		}
		c = c / n;
		centroid[0] = c[0];
		centroid[1] = c[1];
		centroid[2] = c[2];
		currentIndex = 0;
		GLuint currentFaceIndex = 0;
		for (AKMesh::FaceIter f_it = _mesh.faces_begin(); f_it != _mesh.faces_end(); f_it++)
		{
			AKMesh::Normal current_face_normal = _mesh.normal(*f_it);
			for (AKMesh::FaceVertexIter fv_it = _mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++)
			{
				indices[currentIndex] = currentIndex;
				AKMesh::Point current_point = _mesh.point(*fv_it);
				vertices[currentIndex].position = glm::vec3(current_point[0], current_point[1], current_point[2]);
				//vertices[currentIndex].color = glm::vec3(0.5f, 0.0f, 0.0f);
				vertices[currentIndex].color = face_color[currentFaceIndex];
				vertices[currentIndex].normal = glm::vec3(current_face_normal[0], current_face_normal[1], current_face_normal[2]);
				currentIndex++;
			}
			currentFaceIndex = 0;
		}
	}



	glm::vec3 GetColorJET(float v, float vmin, float vmax)
	{
		glm::vec3 c = { 1.0, 1.0, 1.0 }; // white
		float dv;

		if (v < vmin)
			v = vmin;
		if (v > vmax)
			v = vmax;
		dv = vmax - vmin;

		if (v < (vmin + 0.25f * dv)) {
			c.r = 0;
			c.g = 4 * (v - vmin) / dv;
		}
		else if (v < (vmin + 0.5f * dv)) {
			c.r = 0;
			c.b = 1 + 4 * (vmin + 0.25f * dv - v) / dv;
		}
		else if (v < (vmin + 0.75f * dv)) {
			c.r = 4 * (v - vmin - 0.5f * dv) / dv;
			c.b = 0;
		}
		else {
			c.g = 1 + 4 * (vmin + 0.75f * dv - v) / dv;
			c.b = 0;
		}

		return(c);
	}

	/**
	*/
	void loadMeshFaceScalar(AKMesh & _mesh, vector<AKNumber> & face_scalar)
	{
		numVertices = static_cast<GLuint>(_mesh.n_faces() * 3);
		vertices = new Vertex[numVertices];
		numIndices = static_cast<GLuint>(_mesh.n_faces() * 3);
		indices = new GLuint[numIndices];

		GLuint currentIndex = 0;
		AKMesh::Point c(0, 0, 0);
		AKNumber n = 0.0f;
		for (AKMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); v_it++)
		{
			AKMesh::Point current_point = _mesh.point(*v_it);
			c += current_point;
			n++;
			currentIndex++;
		}
		c = c / n;
		centroid[0] = c[0];
		centroid[1] = c[1];
		centroid[2] = c[2];
		currentIndex = 0;
		GLuint currentFaceIndex = 0;
		for (AKMesh::FaceIter f_it = _mesh.faces_begin(); f_it != _mesh.faces_end(); f_it++)
		{
			AKMesh::Normal current_face_normal = _mesh.normal(*f_it);
			for (AKMesh::FaceVertexIter fv_it = _mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++)
			{
				indices[currentIndex] = currentIndex;
				AKMesh::Point current_point = _mesh.point(*fv_it);
				vertices[currentIndex].position = glm::vec3(current_point[0], current_point[1], current_point[2]);
				//vertices[currentIndex].color = GetColorJET(((float)currentIndex)/((float)numIndices), 0.0, 1.0);
				vertices[currentIndex].color = GetColorJET(face_scalar[currentFaceIndex],0.0,1.0);
				vertices[currentIndex].normal = glm::vec3(current_face_normal[0], current_face_normal[1], current_face_normal[2]);
				currentIndex++;
			}
			currentFaceIndex++;
		}
	}
	void loadMeshFaceScalarSparse(AKMesh & _mesh, vector<pair<size_t, AKNumber>> & face_scalar_sparse, float _face_scalar_min, float _face_scalar_max, float default_value)
	{
		cout << "start load" << endl;
		numVertices = static_cast<GLuint>(_mesh.n_faces() * 3);
		vertices = new Vertex[numVertices];
		numIndices = static_cast<GLuint>(_mesh.n_faces() * 3);
		indices = new GLuint[numIndices];

		GLuint currentIndex = 0;
		AKMesh::Point c(0, 0, 0);
		float n = 0.0f;
		for (AKMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); v_it++)
		{
			AKMesh::Point current_point = _mesh.point(*v_it);
			c += current_point;
			n++;
			currentIndex++;
		}
		c = c / n;
		centroid[0] = c[0];
		centroid[1] = c[1];
		centroid[2] = c[2];
		currentIndex = 0;
		GLuint currentFaceIndex = 0;
		size_t idx = 0;
		for (AKMesh::FaceIter f_it = _mesh.faces_begin(); f_it != _mesh.faces_end(); f_it++)
		{
			AKMesh::Normal current_face_normal = _mesh.normal(*f_it);
			float val = default_value;
			if (idx < face_scalar_sparse.size() )
			{
				if (face_scalar_sparse[idx].first == currentFaceIndex)
				{
					val = face_scalar_sparse[currentFaceIndex].second;
					idx++;
				}
			}
			for (AKMesh::FaceVertexIter fv_it = _mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++)
			{
				indices[currentIndex] = currentIndex;
				AKMesh::Point current_point = _mesh.point(*fv_it);
				vertices[currentIndex].position = glm::vec3(current_point[0], current_point[1], current_point[2]);
				vertices[currentIndex].color = GetColorJET(val, _face_scalar_min, _face_scalar_max);
				vertices[currentIndex].normal = glm::vec3(current_face_normal[0], current_face_normal[1], current_face_normal[2]);
				currentIndex++;
			}
			currentFaceIndex++;
		}
		cout << "finish load" << endl;
	}

	void loadMeshFaceScalar(AKMesh & _mesh, vector<AKNumber> & face_scalar,float _face_scalar_min,float _face_scalar_max)
	{
		numVertices = static_cast<GLuint>(_mesh.n_faces() * 3);
		vertices = new Vertex[numVertices];
		numIndices = static_cast<GLuint>(_mesh.n_faces() * 3);
		indices = new GLuint[numIndices];

		GLuint currentIndex = 0;
		AKMesh::Point c(0, 0, 0);
		float n = 0.0f;
		for (AKMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); v_it++)
		{
			AKMesh::Point current_point = _mesh.point(*v_it);
			c += current_point;
			n++;
			currentIndex++;
		}
		c = c / n;
		centroid[0] = c[0];
		centroid[1] = c[1];
		centroid[2] = c[2];
		currentIndex = 0;
		GLuint currentFaceIndex = 0;
		for (AKMesh::FaceIter f_it = _mesh.faces_begin(); f_it != _mesh.faces_end(); f_it++)
		{
			AKMesh::Normal current_face_normal = _mesh.normal(*f_it);
			for (AKMesh::FaceVertexIter fv_it = _mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++)
			{
				indices[currentIndex] = currentIndex;
				AKMesh::Point current_point = _mesh.point(*fv_it);
				vertices[currentIndex].position = glm::vec3(current_point[0], current_point[1], current_point[2]);
				//vertices[currentIndex].color = GetColorJET(((float)currentIndex)/((float)numIndices), 0.0, 1.0);
				vertices[currentIndex].color = GetColorJET(face_scalar[currentFaceIndex], _face_scalar_min, _face_scalar_max);
				vertices[currentIndex].normal = glm::vec3(current_face_normal[0], current_face_normal[1], current_face_normal[2]);
				currentIndex++;
			}
			currentFaceIndex = 0;
		}
	}

	void updateFaceScalar(vector<AKNumber> & face_scalar)
	{
		int temp = 0;
		int cont = 0;
		for (GLuint i  = 0; i < numIndices; i++)
		{
			vertices[i].color = GetColorJET(face_scalar[cont], 0.0, 1.0);
			if (temp == 2)
			{
				cont++;
				temp = 0;
			}
			else
				temp++;
		}
	}

	void updateFaceScalar(vector<AKNumber> & face_scalar, float _face_scalar_min, float _face_scalar_max)
	{
		int temp = 0;
		int cont = 0;
		for (GLuint i = 0; i < numIndices; i++)
		{
			vertices[i].color = GetColorJET(face_scalar[cont], _face_scalar_min, _face_scalar_max);
			if (temp == 2)
			{
				cont++;
				temp = 0;
			}
			else
				temp++;
		}
	}

	void updateFaceScalarSparse(vector<pair<size_t, AKNumber>> & face_scalar_sparse, float _face_scalar_min, float _face_scalar_max, float _default_value)
	{
		int temp = 0;
		int cont = 0;
		int idx = 0;
		for (GLuint i = 0; i < numIndices; i++)
		{
			bool flag = false;
			if (idx < face_scalar_sparse.size())
			{
				if (cont == face_scalar_sparse[idx].first)
				{
					vertices[i].color = GetColorJET(face_scalar_sparse[idx].second, _face_scalar_min, _face_scalar_max);
					flag = true;
				}
				else
					vertices[i].color = GetColorJET(_default_value, _face_scalar_min, _face_scalar_max);
			}
			else
			{
				vertices[i].color = GetColorJET(_default_value, _face_scalar_min, _face_scalar_max);
			}
			if (temp == 2)
			{
				if (flag)
					idx++;
				cont++;
				temp = 0;
			}
			else
				temp++;
		}
	}

    /**
     * @brief vertexBufferSize - get size of the array of vertices
     * @return - size of array of vertices
     */
	GLsizeiptr vertexBufferSize() const
	{
		return numVertices * sizeof(Vertex);
	}
    /**
     * @brief indexBufferSize - get size of the array of indices
     * @return - size of array of indices
     */
	GLsizeiptr indexBufferSize() const
	{
		return numIndices * sizeof(GLuint);
	}
    /**
     * @brief clear - erase arrays and reinitialize values
     */
	void clear()
	{
		delete[] vertices;
		delete[] indices;
		vertices = 0;
		indices = 0;
		numVertices = numIndices = 0;
		centroid = glm::vec3(0.0f,0.0f,0.0f);
	}
};

/** @} */

#endif // MYSHAPE_H
