#include "core/iomesh.h"

bool import_mesh(AKMesh &mesh, string filename)
{
    OpenMesh::IO::Options opt;
    //opt += OpenMesh::IO::Options::VertexColor;
    if(!OpenMesh::IO::read_mesh(mesh, filename, opt))
        return false;
    return true;
}

bool export_mesh(AKMesh &mesh, string filename)
{
    OpenMesh::IO::Options opt;
    //opt += OpenMesh::IO::Options::VertexColor;
    return OpenMesh::IO::write_mesh(mesh, filename,opt);
}
