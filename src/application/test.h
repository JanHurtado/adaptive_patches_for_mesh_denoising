#include "core/util.h"
#include "core/iomesh.h"
#include "core/denoising.h"
#include "core/noise.h"
#include "core/metric.h"

typedef AKNumber (*t_testFunctionPtr)(AKMesh &);

AKNumber compare_positions(AKMesh & m1, AKMesh & m2)
{
	AKNumber error = 0;
	for (AKMesh::VertexIter v_it = m1.vertices_begin(); v_it != m1.vertices_end(); v_it++)
	{
		int index = v_it->idx();
		error += (m1.point(AKMesh::VertexHandle(index)) - m2.point(AKMesh::VertexHandle(index))).length();
	}
	return error;
}

AKNumber compare_normals(AKMesh & m1, AKMesh & m2)
{
	AKNumber error = 0;
	for (AKMesh::VertexIter v_it = m1.vertices_begin(); v_it != m1.vertices_end(); v_it++)
	{
		int index = v_it->idx();
		error += (m1.normal(AKMesh::VertexHandle(index)) - m2.normal(AKMesh::VertexHandle(index))).length();
	}
	return error;
}

AKNumber test_update_vertex_positions(AKMesh & mesh)
{
	AKMesh res = mesh;
	vector<AKMesh::Normal> normals;
	compute_all_face_normals(res, normals);
	update_vertex_positions(res,normals,10,false);
	return compare_positions(res, mesh) + compare_normals(res, mesh);
}

AKNumber test_compute_average_edge_length(AKMesh & mesh)
{
	return compute_average_edge_length(mesh);
}

AKNumber test_compute_area(AKMesh & mesh)
{
	return compute_area(mesh);
}

AKNumber test_compute_radius(AKMesh & mesh)
{
	return compute_radius(mesh,1.0);
}

AKNumber test_compute_vertex_area(AKMesh & mesh)
{
	vector<AKNumber> areas;
	compute_all_face_areas(mesh, areas);
	if (mesh.n_vertices() != 0)
		return compute_vertex_area(mesh, AKMesh::VertexHandle(0) , areas);
	else
		return 0.0f;
}


vector<string> testNames = { "updateVertexPositions", "computeAverageEdgeLength", "computeArea", "computeRadius", "computeVertexArea" };
vector<t_testFunctionPtr> testFunctions = { &test_update_vertex_positions, &test_compute_average_edge_length,&test_compute_area, &test_compute_radius, &test_compute_vertex_area};
vector<AKNumber> tolerance = {0.001f , 0.001f, 0.001f, 0.001f, 0.001f};
vector<AKNumber> expectedValuesTestCase0 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
vector<AKNumber> expectedValuesTestCase1 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
vector<AKNumber> expectedValuesTestCase2 = { 0.0f, 1.138f, 0.5f, 0.0f, 0.1667f };
vector<AKNumber> expectedValuesTestCase3 = { 0.0f, 1.0828f, 1.0f, 0.4714f, 0.3333f };

AKMesh testCase0; //empty
AKMesh testCase1; //single point
AKMesh testCase2; //single trianle
AKMesh testCase3; //plane (two triangles)

void initialize_test_cases()
{
	testCase0.clear();
	testCase1.clear();
	import_mesh(testCase1, "TestCase1.off");
	testCase1.request_face_normals();
	testCase1.update_normals();
	import_mesh(testCase2, "TestCase2.off");
	testCase2.request_face_normals();
	testCase2.update_normals();
	import_mesh(testCase3, "TestCase3.off");
	testCase3.request_face_normals();
	testCase3.update_normals();
}

bool OK(AKNumber value, AKNumber expected_value, AKNumber tolerance)
{
	return (value >= (expected_value - tolerance) && value <= (expected_value + tolerance));
}

void run_tests()
{
	for (int i = 0; i < testNames.size(); i++)
	{
		cout << testNames[i] << endl;
		cout << "TestCase0 -- " << "value: " << (*(testFunctions[i]))(testCase0) << " / expected value: " << expectedValuesTestCase0[i] << " / OK: " << OK((*(testFunctions[i]))(testCase0), expectedValuesTestCase0[i], tolerance[i]) << endl;
		cout << "TestCase1 -- " << "value: " << (*(testFunctions[i]))(testCase1) << " / expected value: " << expectedValuesTestCase1[i] << " / OK: " << OK((*(testFunctions[i]))(testCase1), expectedValuesTestCase1[i], tolerance[i]) << endl;
		cout << "TestCase2 -- " << "value: " << (*(testFunctions[i]))(testCase2) << " / expected value: " << expectedValuesTestCase2[i] << " / OK: " << OK((*(testFunctions[i]))(testCase2), expectedValuesTestCase2[i], tolerance[i]) << endl;
		cout << "TestCase3 -- " << "value: " << (*(testFunctions[i]))(testCase3) << " / expected value: " << expectedValuesTestCase3[i] << " / OK: " << OK((*(testFunctions[i]))(testCase3), expectedValuesTestCase3[i], tolerance[i]) << endl;
		cout << endl;
	}
}

vector<AKMesh::Normal> compute_guidance_normals(AKMesh &mesh)
{
	bool include_central_face = 1;
	FaceNeighborType face_neighbor_type = kRadiusBased;
	vector<vector<AKMesh::FaceHandle> > all_guided_neighbors((int)mesh.n_faces());
	compute_all_guidance_signal_neighbors(mesh, all_guided_neighbors);

	vector<AKNumber> face_areas((int)mesh.n_faces());
	vector<AKMesh::Point> face_centroids((int)mesh.n_faces());
	vector<AKMesh::Normal> previous_normals((int)mesh.n_faces());
	vector<AKMesh::Normal> guided_normals((int)mesh.n_faces());
	vector<pair<AKNumber, AKMesh::Normal> > consistencies_and_mean_normals((int)mesh.n_faces());

	compute_all_face_centroids(mesh, face_centroids);
	compute_all_face_areas(mesh, face_areas);
	compute_all_face_normals(mesh, previous_normals);

	compute_guided_normals(mesh, all_guided_neighbors, face_areas, previous_normals, consistencies_and_mean_normals, guided_normals);
	return guided_normals;
}

vector<AKMesh::Normal> compute_adaptive_kernels_normals(AKMesh & mesh)
{
	AKNumber alpha = 1.0; 
	AKNumber beta = 2; 
	AKNumber gamma = 0.8;
	AKNumber delta = 8;
	AKNumber k = 20;
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

	vector<AKMesh::Normal> adaptive_kernel_normals(_patches_sparse.size());
	for (int it = 0; it < 1; it++)
	{
		for (int i = 0; i < _patches_sparse.size(); i++)
		{
			AKMesh::Normal avg_normal = AKMesh::Normal(0, 0, 0);
			for (int j = 0; j < _patches_sparse[i].size(); j++)
			{
				avg_normal += normals[_patches_sparse[i][j].first] * _patches_sparse[i][j].second * areas[_patches_sparse[i][j].first];
			}
			adaptive_kernel_normals[i] = avg_normal.normalize_cond();
		}
		normals = adaptive_kernel_normals;
	}



	return adaptive_kernel_normals;
}


AKNumber angular_error(AKMesh::Normal normal1, AKMesh::Normal normal2)
{
	AKNumber dot_product = normal1 | normal2;
	dot_product = min(1.0, max(dot_product, -1.0));
	return acos(dot_product) * 180.0 / M_PI;
}

void comparison1()
{
	AKMesh original_mesh;
	import_mesh(original_mesh, "fandisk.off");
	AKMesh noisy_mesh;
	import_mesh(noisy_mesh, "fandisk_noisy0.1.off");
	vector<AKMesh::Normal> guided_normals = compute_guidance_normals(noisy_mesh);
	vector<AKMesh::Normal> adaptive_kernel_normals = compute_adaptive_kernels_normals(noisy_mesh);
	vector<AKMesh::Normal> ground_truth;
	vector<AKNumber> areas;
	compute_all_face_areas(original_mesh, areas);
	compute_all_face_normals(original_mesh, ground_truth);
	AKNumber mean_square_angle_error_guided = 0.0;
	AKNumber mean_square_angle_error_robust = 0.0;
	AKNumber mesh_total_area = 0.0;
	AKNumber error_robust = 0.0;
	AKNumber error_guided = 0.0;
	for (int i = 0; i < ground_truth.size(); i++)
	{
		mean_square_angle_error_guided += angular_error(guided_normals[i],ground_truth[i]);
		mean_square_angle_error_robust += angular_error(adaptive_kernel_normals[i], ground_truth[i]);
		AKNumber temp_guided = pow((guided_normals[i] - ground_truth[i]).norm(), 2.0);
		AKNumber temp_robust = pow((adaptive_kernel_normals[i] - ground_truth[i]).norm(), 2.0);
		error_guided += temp_guided*areas[i];
		error_robust += temp_robust*areas[i];
		mesh_total_area += areas[i];
	}
	error_guided /= mesh_total_area;
	error_robust /= mesh_total_area;
	mean_square_angle_error_guided = mean_square_angle_error_guided / (AKNumber)ground_truth.size();
	mean_square_angle_error_robust = mean_square_angle_error_robust / (AKNumber)ground_truth.size();
	cout << "guided: " << mean_square_angle_error_guided << endl << "robust: " << mean_square_angle_error_robust << endl;
	cout << "l2 guided: " << error_guided << endl << "l2 robust: " << error_robust << endl;

	AKMesh mesh_guided = noisy_mesh;
	update_vertex_positions(mesh_guided, guided_normals, 20, true);
	mesh_guided.request_face_normals();
	mesh_guided.request_vertex_normals();
	mesh_guided.update_normals();
	export_mesh(mesh_guided, "t_guided.off");

	AKMesh mesh_robust = noisy_mesh;
	update_vertex_positions(mesh_robust, adaptive_kernel_normals, 20, true);
	mesh_robust.request_face_normals();
	mesh_robust.request_vertex_normals();
	mesh_robust.update_normals();
	export_mesh(mesh_robust, "t_robust.off");


	noisy_mesh.request_face_normals();
	noisy_mesh.request_vertex_normals();
	noisy_mesh.update_normals();
	export_mesh(noisy_mesh, "t_noisy.off");

}


void comparison2()
{
	AKMesh original_mesh;
	import_mesh(original_mesh, "fandisk.off");
	AKMesh bn_mesh;
	import_mesh(bn_mesh, "fandisk_bn.off");
	AKMesh test_mesh;
	import_mesh(test_mesh, "fandisk_test.off");
	AKNumber bn_error = compute_mean_square_angle_error(original_mesh, bn_mesh);
	AKNumber test_error = compute_mean_square_angle_error(original_mesh, test_mesh);
	cout << "ERRORS: " << bn_error << " " << test_error << endl;
}

void comparison3()
{
	//vector<string> algorithm_names = { "bilateral", "nonIterative", "fast", "bilateralNormal", "l0", "guided", "testf" };
	vector<string> algorithm_names = { "testf"};
	string mesh_name = "dragon100000";
	string path = "C:/Users/hurtado/Desktop/results/"+mesh_name+"/";
	string original_file_name = path + mesh_name + ".off";
	AKMesh original_mesh;
	import_mesh(original_mesh, original_file_name);
	string results_file_name = path + "results_"+mesh_name+".txt";
	ofstream results(results_file_name.c_str());

	PointCloud v_cloud, f_cloud;
	generate_vertex_point_cloud(original_mesh, v_cloud);
	generate_face_centroid_point_cloud(original_mesh, f_cloud);
	cout << "point clouds generadas YA" << endl;
	AKKDTree   index_vertices(3 /*dim*/, v_cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index_vertices.buildIndex();

	AKKDTree   index_faces(3 /*dim*/, f_cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index_faces.buildIndex();

	vector<AKNumber> original_mesh_curvature;
	compute_vertex_mean_curvature(original_mesh, original_mesh_curvature);

	AKNumber original_vol = compute_volume(original_mesh);
	AKNumber original_area =  compute_area(original_mesh);

	results << "A:MeanVertexDistance B:L2VertexBased C:MeanQuadric D:MSAE E:L2NormalBased F:Tangential G:MeanDiscreteCurvature H:Area I:Volume" << endl;

	for (int i = 0; i < algorithm_names.size(); i++)
	{
		string currentResult_file_name = path + mesh_name + "_" + algorithm_names[i] + ".off";
		cout << currentResult_file_name << endl;
		AKMesh currentMesh;
		import_mesh(currentMesh, currentResult_file_name);
		vector<AKNumber> current_mesh_curvature;
		compute_vertex_mean_curvature(currentMesh, current_mesh_curvature);

		results << algorithm_names[i] << endl;
		results << "A:" << compute_mean_vertex_distance_error(original_mesh, index_vertices, index_faces, currentMesh) << " ";
		results << "B:" << compute_L2_vertex_based_error(original_mesh, index_vertices, index_faces, currentMesh) << " ";
		results << "C:" << compute_mean_quadric_error(original_mesh, index_vertices, currentMesh) << " ";
		results << "D:" << compute_mean_square_angle_error_2(original_mesh, index_vertices, index_faces, currentMesh) << " ";
		//results << computeMeanSquareAngleError(original_mesh, currentMesh) << " ";
		results << "E:" << compute_L2_normal_based_error(original_mesh, index_vertices, index_faces, currentMesh) << " ";
		results << "F:" << compute_mean_tangential_error(original_mesh, index_vertices, currentMesh) << " ";
		results << "G:" << compute_mean_discrete_curvature_error(original_mesh, index_vertices, original_mesh_curvature, currentMesh, current_mesh_curvature) << " ";
		results << "H:" << compute_area_error(original_area, currentMesh) << " ";
		results << "I:" << compute_volume_error(original_vol, currentMesh) << endl;
	}
}


void latexCodeGenerator(const char * inputFileName, const char * outputFileName)
{
	ifstream iFile(inputFileName);
	ofstream oFile(outputFileName);
	std::string  current_line;
	std::getline(iFile, current_line);
	std::stringstream  linestream(current_line);
	string token;

	vector<string> algorithm_names;
	vector<string> metric_names;

	while (linestream>>token)
	{
		metric_names.push_back(token.substr(2,token.size()-2));
	}
	int num_metrics = metric_names.size();
	vector<vector<string>> values(num_metrics);
	vector<vector<AKNumber>> valuesD(num_metrics);

	int i = 0;
	int num_algorithms = 0;
	while (std::getline(iFile, current_line))
	{
		linestream.clear();
		linestream.str(current_line);
		if (i % 2 == 0)
		{
			linestream >> token;
			algorithm_names.push_back(token);
		}
		else
		{
			while (linestream >> token)
			{
				string val_string = token.substr(2, token.size() - 2);
				AKNumber val = stod(val_string);
				stringstream t_stream;
				t_stream << fixed << setprecision(6) << val;
				string val_res = t_stream.str();
				values[num_algorithms].push_back(val_res);
				valuesD[num_algorithms].push_back(val);
			}
			num_algorithms++;
		}
		i++;
		// Store words
	}


	metric_names = { "Mean Vertex Distance" , "L2 Vertex Based" , "Mean Quadric" , "MSAE" , "L2 Normal Based" , "Tangential" , "Mean Discrete Curvature" , "Area Error" , "Volume Error" };

	for (int i = 0; i < algorithm_names.size(); i++)
	{
		if (algorithm_names[i] == "bilateral") algorithm_names[i] = "\\cite{FDC03}";
		if (algorithm_names[i] == "nonIterative") algorithm_names[i] = "\\cite{JDD03}";
		if (algorithm_names[i] == "fast") algorithm_names[i] = "\\cite{SRML07}";
		if (algorithm_names[i] == "bilateralNormal") algorithm_names[i] = "\\cite{ZFAT11}";
		if (algorithm_names[i] == "l0") algorithm_names[i] = "\\cite{HS13}";
		if (algorithm_names[i] == "guided") algorithm_names[i] = "\\cite{ZDZBL15}";
		if (algorithm_names[i] == "nvt") algorithm_names[i] = "\\cite{YRP16}";
		if (algorithm_names[i] == "robust") algorithm_names[i] = "\\cite{YRP17}";
		if (algorithm_names[i] == "testf") algorithm_names[i] = "Our";

	}

	vector<int> mins_index(num_metrics);
	for (int i = 0; i < num_metrics; i++)
	{
		AKNumber min = 9999.99999;
		int min_index = 0;
		for (int j = 0; j < num_algorithms; j++)
		{
			AKNumber temp = valuesD[j][i];
			if (temp < min)
			{
				min = temp;
				min_index = j;
			}
		}
		mins_index[i] = min_index;
	}

	string allMetricNames = "";
	for (int i = 0; i < metric_names.size(); i++)
	{
		allMetricNames += " & " + metric_names[i];
	}
	oFile << "\\begin{tabular}{ m{1.5cm} m{1cm} m{1cm} m{1cm} m{1cm} m{1cm} m{1cm} m{1cm} m{1cm} m{1cm} } \n";
	oFile << allMetricNames << "\\\\ \n \\hline \n";
	for (int i = 0; i < num_algorithms; i++)
	{
		oFile << algorithm_names[i];
		for (int j = 0; j < num_metrics; j++)
		{
			if (mins_index[j] == i)
			{
				oFile << " & \\textbf{" << values[i][j] << "}";
			}
			else
			{
				oFile << " & " << values[i][j];
			}

		}
		oFile << " \\\\ \n ";
	}
	oFile << "\\end{tabular}\n";
	/*while (std::computeline(iFile, current_line))
	{
		std::stringstream  linestream(current_line);
		string token;
		while (is.compute(c))          // loop computeting single characters
		std::string word1, word2, word3;
		line >> word1 >> word2 >> word3;

		// Store words
	}*/

}