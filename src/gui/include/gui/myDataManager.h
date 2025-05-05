#ifndef MYDATAMANAGER_H
#define MYDATAMANAGER_H

#include <core/iomesh.h>
#include <core/mesh.h>
#include <visualization/myShape.h>
#include <core/util.h>

/**
 * @brief The myDataManager class - data manager for mesh smoothing application
 */
class myDataManager
{
public:
    /**
     * @brief input_mesh - input mesh (left side in the interface)
     */
	AKMesh input_mesh;
    /**
     * @brief output_mesh - output mesh/resulting mesh/edited mesh (right side in the interface)
     */
	AKMesh output_mesh;

    /**
     * @brief input_mesh_shape - input mesh shape data (for visualization)
     */
	ShapeData input_mesh_shape;
    /**
     * @brief output_mesh_shape - output mesh shape data (for visualization)
     */
	ShapeData output_mesh_shape;

    /**
     * @brief selection - shape data for current selection (set of vertices)
     */
	ShapeData selection;

    /**
     * @brief myDataManager - default constructor
     */
	myDataManager();
    /**
     * @brief myDataManager - destructor
     */
	~myDataManager();
    /**
     * @brief loadInputMesh - load a mesh file supported by OpenMesh library
     * @param fileName mesh file name
     */
    void loadInputMesh(string & fileName);
    /**
     * @brief saveOutputMesh - save the output mesh in a file supported bu OpenMesh library
     * @param fileName mesh file name
     */
	void saveOutputMesh(string & fileName);
    /**
     * @brief updateInputShape - generate the corresponding shape data of the input mesh
     */
	void updateInputShape();
    /**
     * @brief updateOutputShape - generate the corresponding shape data of the output mesh
     */
	void updateOutputShape();
    /**
     * @brief updateShapes - generate the corresponding shape data for both meshes (input mesh and output mesh)
     */
	void updateShapes();
    /**
     * @brief updateOutputSelection - generate the selection shape for a set of given vertex indices
     * @param indices target set of vertex indices
     */
	void updateOutputSelection(vector<size_t> & indices);
    /**
     * @brief setOutputAsInput - set output mesh as input mesh (save all changes)
     */
	void setOutputAsInput();
    /**
     * @brief reinitialize - set input mesh as output mesh (undone all changes)
     */
	void reinitialize();


	void setFaceScalarFunction(int idx);

	vector<vector<pair<size_t, AKNumber>>> face_scalar_functions_sparse;
	vector<vector<pair<size_t, AKNumber>>> face_scalar_functions_sparse_2;

	float face_scalar_default;
	float face_scalar_min;
	float face_scalar_max;

	// single patch

	vector<pair<size_t, AKNumber>> sp_patch_function;
	vector<AKNumber> sp_areas;
	vector<AKMesh::Point> sp_centroids;
	vector<AKMesh::Normal> sp_normals;

	AKMesh groundtruth;

	void setSinglePatchFunction();
	void genGeodesicData(std::vector<double> & points, std::vector<unsigned> & faces);

	AKNumber scale_factor;

};
#endif // MYDATAMANAGER_H
