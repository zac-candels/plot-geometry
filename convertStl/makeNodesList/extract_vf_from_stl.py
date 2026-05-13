import numpy as np
import stl
import os

cwd = os.getcwd()

my_stl_file= cwd + "/Jumping.STL"


def parse_stl(stl_file):

    # Load the STL file
    mesh_data = stl.mesh.Mesh.from_file(stl_file)

    # Extract vertices and faces
    vertices = mesh_data.vectors.reshape((-1, 3))
    faces = np.arange(len(vertices)).reshape((-1, 3))

    return vertices, faces
  
def write_vertex_face_data(vertices,faces):

    np.savetxt('verts.in', vertices)
    np.savetxt('faces.in', faces,fmt='%d')
	   
    return
    
#Main
vertices,faces=parse_stl(my_stl_file)
print(vertices)
print(faces)

write_vertex_face_data(vertices,faces)
