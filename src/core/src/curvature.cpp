#include "core/curvature.h"

void compute_vertex_gaussian_curvature(AKMesh & mesh, vector<AKNumber> & curvature)
{
    curvature.resize(mesh.n_vertices());
    vector<AKNumber> face_areas;
    compute_all_face_areas(mesh, face_areas);
    for (AKMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
    {
        AKNumber total_vertex_area = 0.0;
        AKNumber total_vertex_angle = 0.0;
        for (AKMesh::VertexFaceIter vf_it = mesh.vf_iter(*v_it); vf_it.is_valid() ; ++vf_it)
        {
            AKNumber current_area = face_areas[vf_it->idx()];
            AKNumber current_angle;
            compute_face_vertex_angle(mesh,*vf_it,*v_it,current_angle);
            total_vertex_area += current_area;
            total_vertex_angle += current_angle;
        }
        AKNumber gaussian_curvature = 2.0*PI-total_vertex_angle;
        curvature[v_it->idx()]=gaussian_curvature/(total_vertex_area*0.333);
    }
}

void compute_vertex_mean_curvature(AKMesh & mesh, vector<AKNumber> & curvature)
{
    curvature.resize(mesh.n_vertices());
    for (AKMesh::VertexIter v_it=mesh.vertices_sbegin(); v_it!=(mesh.vertices_end()); ++v_it)
    {
        AKMesh::Point   prev_pt, next_pt, present_pt ;
        OpenMesh::Vec3f v1_beta, v2_beta,v1_alpha,v2_alpha;
        AKMesh::Point lapl_bel(0,0,0);
        AKNumber alpha=0.0, beta=0.0; //Angles
        AKNumber cotan_sum;
        AKNumber total_area=0.0;
        AKNumber area=0.0;
        int is_inside=0;//See if the loop below is entered
        for (AKMesh::VertexVertexIter vv_it=mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
        {
            AKMesh::VertexFaceIter face = mesh.vf_iter(*vv_it);
            if(face.is_valid()){
                is_inside=1;
                AKMesh::VertexVertexIter vv_it_prev=vv_it;
                AKMesh::VertexVertexIter vv_it_next=vv_it;
                prev_pt=mesh.point(*(--vv_it_prev));//Previous point to the present-one ring neighbour
                next_pt=mesh.point(*(++vv_it_next));//Next point to the present-one ring neighbour
                present_pt=mesh.point(*vv_it);

                //Determining alpha
                v1_alpha= mesh.point(*v_it)-prev_pt;
                v2_alpha= present_pt-prev_pt;
                v1_alpha.normalize_cond(); v2_alpha.normalize_cond();
                alpha = acos(v1_alpha[0] * v2_alpha[0] + v1_alpha[1] * v2_alpha[1] + v1_alpha[2] * v2_alpha[2]);//Angle Alpha

                //Determining beta
                v1_beta= mesh.point(*v_it)-next_pt;
                v2_beta= present_pt-next_pt;

                area = 0.5 * v1_beta.norm() * v2_beta.norm();//area= (1/2)*a*b

                v1_beta.normalize_cond(); v2_beta.normalize_cond();
                beta = acos(v1_beta[0] * v2_beta[0] + v1_beta[1] * v2_beta[1] + v1_beta[2] * v2_beta[2]);//Angle Beta

                area=area*sin(beta);//now area=(1/2)*a*b * sin(angle-between-a-and-b)
                total_area+=area;//Total area
                cotan_sum=tan((PI/2)-alpha)+tan((PI/2)-beta);//cot(aplha)+cot(beta)
                lapl_bel+=( (cotan_sum) * (  present_pt -mesh.point(*v_it) ) );
            }
        }
        if(is_inside==1){
            total_area=(1.0/3.0)*total_area;
			AKNumber val = (lapl_bel / (2.0*total_area)).norm()*0.5;
			curvature[v_it->idx()] = min(max(-2.5, val),2.5);
			
        }
    }
}

void compute_vertex_principal_curvatures_sum(AKMesh & mesh,vector<AKNumber> & mean_curvature,vector<AKNumber> & gaussian_curvature,vector<AKNumber> & principal_curvature)
{
    principal_curvature.resize(mesh.n_vertices());
    for (AKMesh::VertexIter v_it=mesh.vertices_sbegin(); v_it!=(mesh.vertices_end()); ++v_it)
    {
        AKNumber temp = pow(mean_curvature[v_it->idx()],2.0)-gaussian_curvature[v_it->idx()];
        if (temp<0.0)
            temp=0;
        AKNumber k1 = mean_curvature[v_it->idx()] + sqrt(temp);
        AKNumber k2 = mean_curvature[v_it->idx()] - sqrt(temp);
        principal_curvature[v_it->idx()] = abs(k1)+abs(k2);
    }
}

void compute_all_curvatures(AKMesh & mesh,vector<AKNumber> & mean_curvature,vector<AKNumber> & gaussian_curvature,vector<AKNumber> & principal_curvature)
{
    compute_vertex_mean_curvature(mesh,mean_curvature);
    compute_vertex_gaussian_curvature(mesh,gaussian_curvature);
    compute_vertex_principal_curvatures_sum(mesh,mean_curvature,gaussian_curvature,principal_curvature);
}
