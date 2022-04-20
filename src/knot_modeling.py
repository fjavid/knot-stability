import sys
sys.path.insert(0, "/home/Codim-IPC/Python")
sys.path.insert(0, "/home/Codim-IPC/build")
# print(sys.path)
import Drivers
from Drivers.SimulationBase import make_directory
from JGSL import *
import os
import numpy as np
# import meshio
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# sys.path.insert(0, "/home/Codim-IPC/Python")
from flexible_rod import FlexibleRod

# def blender_to_cipc(file_dir, file_name):
#     fullname = file_dir + file_name + ".obj"
#     with open(fullname, 'r') as f:
#         lines = f.readlines()
#     for idx, line in enumerate(lines):
#         sline = line.split(' ')
#         if sline[0] == 'l':
#             nline = 's' + line[1:]
#             lines[idx] = nline
#     cipc_fullname = file_dir + file_name + "_cipc.obj"
#     with open(cipc_fullname, 'w') as f:
#         for line in lines:
#             # print(line)
#             f.write(line)

# file_dir = "/home/project/input/"
# file_name = "f1x1_1"
# if __name__ == "__main__":
#     blender_to_cipc(file_dir, file_name)



def make_coords_sys(axis, app_n1):
    # axis = vert_normals[0]
    # app_norm = np.array([1, 0, 0], dtype=float)
    axis /= np.linalg.norm(axis)
    n2 = np.cross(axis, app_n1)
    if np.linalg.norm(n2) > 1.0e-6:
        n1 = np.cross(n2, axis)
        n1 /= np.linalg.norm(n1)
        n2 /= np.linalg.norm(n2)
    else:
        n2 = np.random.shuffle(axis)
        n1 = np.cross(n2, axis)
        n1 /= np.linalg.norm(n1)
        n2 /= np.linalg.norm(n2)
    return np.stack([axis, n1, n2])

def sweep_curve(file_dir, file_name, radius):
    fullname = file_dir + file_name + ".obj"
    curv_verts = []
    curv_segs = []
    with open(fullname) as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        sline = line.split(' ')
        if sline[0] == 'v':
            curv_verts.append(np.array([float(sline[1]), float(sline[2]), float(sline[3])], dtype=float))
        if sline[0] == 'f' or sline[0] == 'l':
            curv_segs.append(np.array([sline[1], sline[2]], dtype=int))

    seg_axes = []
    for seg in curv_segs:
        seg_arr = curv_verts[seg[1]-1]-curv_verts[seg[0]-1]
        seg_axes.append(np.array(seg_arr/np.linalg.norm(seg_arr), dtype=float))

    vert_normals = np.zeros((len(curv_verts), 3), dtype=float)
    for idx, seg in enumerate(curv_segs):
        vert_normals[seg[0]-1] += seg_axes[idx]
        vert_normals[seg[1]-1] += seg_axes[idx]

    for idx in range(len(curv_verts)):
        vert_normals[idx] /= np.linalg.norm(vert_normals[idx])

    vert_coords = []
    for idx in range(len(curv_verts)):
        app_n1 = np.array([1.0, 0, 0])
        if idx > 0:
            app_n1 = vert_coords[idx-1][1,:]
        vert_coords.append(make_coords_sys(vert_normals[idx], app_n1))

    surf_verts = []
    cs_res = 8
    for v_id, center_v in enumerate(curv_verts):
        for cs_idx in range(cs_res):
            theta = 2.0 * np.pi * cs_idx /cs_res
            n1 = vert_coords[v_id][1]
            n2 = vert_coords[v_id][2]
            rel_v = radius * (np.cos(theta)*n1 + np.sin(theta)*n2)
            surf_verts.append(center_v+rel_v)

    surf_faces = []
    for v_id in range(len(curv_verts)-1):
        idx_s = v_id * cs_res
        for cs_idx in range(cs_res-1):
            surf_faces.append([idx_s+cs_idx, idx_s+cs_idx+cs_res, idx_s+cs_idx+1])
            surf_faces.append([idx_s+cs_idx+cs_res, idx_s+cs_idx+cs_res+1, idx_s+cs_idx+1])
        surf_faces.append([idx_s+cs_res-1, idx_s+cs_res-1+cs_res, idx_s])
        surf_faces.append([idx_s+cs_res-1+cs_res, idx_s+cs_res, idx_s])

    for idx, face in enumerate(surf_faces):
        surf_faces[idx] = [face[0]+1, face[1]+1, face[2]+1]

    fullname = fullname = file_dir + file_name + "3D.obj"
    with open(fullname, "w") as f:
        for point in surf_verts:
            f.write('v %f %f %f\n' % (point[0], point[1], point[2]) )
        for face in surf_faces:
            f.write('f %i %i %i\n' % (face[0], face[1], face[2]) )
    

def seg_to_curve_Obj(file_dir, file_name):
    fullname = file_dir + file_name + ".obj"
    with open(fullname) as f:
        lines = f.readlines()
    for idx, line in enumerate(lines):
        if line[0] == 'f':
            nline = 'l '
            # print(line)
            sline = line.split(' ')
            for pl in sline[1:-1]:
                nline += pl + ' '
            lines[idx] = nline + '\n'
        # print(lines[idx])
    # name2 = fullname[:-4] + "_1.obj"
    with open(fullname, "w") as f:
        for line in lines:
            f.write(line)

if __name__ == "__main__":
    sim = Drivers.FEMDiscreteShellBase("double", 3)
    sim.output_folder = base_dir+"/output/1x1/" #  + os.path.splitext(os.path.basename(sys.argv[0]))[0] + "/"
    # print(dir(Drivers))
    make_directory(sim.output_folder)
    size = '20'
    if len(sys.argv) > 1:
        size = sys.argv[1]

    # N = 25
    # if len(sys.argv) > 2:
    #     N = int(sys.argv[2])

    sim.mu = 0.0

    # sim.add_shell_with_scale_3D(base_dir+"/input/cylinder3K.obj", Vector3d(0, -0.15, 0), Vector3d(0.05, 0.5, 0.05), \
    #     Vector3d(0.008, -0.099, 0), Vector3d(1, 0, 0), 90)
    # # v, rotCenter, rotAxis, angVelDeg
    # sim.set_DBC(Vector3d(-0.1, -0.1, -0.1), Vector3d(1.1, 1.1, 1.1), 
    #     Vector3d(0, 0, 0), Vector3d(0, 0, 0), Vector3d(0, 0, 0), 0)


    # sim.add_shell_with_scale_3D(base_dir+"/input/f1x1_3d.obj", Vector3d(0.0, 0.0, 0), Vector3d(1, 1, 1),
    #         Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    

    sim.add_rod_3D(base_dir+"/input/f1x1_1_cipc.obj", Vector3d(0.0, 0.0, 0),
            Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, Vector3d(1, 1, 1))

    # top-left
    sim.set_DBC(Vector3d(-0.1, 0.999, -0.1), Vector3d(0.5, 1.1, 1.1), 
        Vector3d(-0.05, -0.025, 0), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    # bottom-left
    sim.set_DBC(Vector3d(-0.1, -0.1, -0.1), Vector3d(0.5, 0.001, 1.1), 
        Vector3d(-0.05, 0.025, 0), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    
    sim.add_rod_3D(base_dir+"/input/f1x1_2_cipc.obj", Vector3d(0.0, -0.0, 0),
            Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, Vector3d(1, 1, 1))
   
    TOL = 0.00001

    # bottom-right
    sim.set_DBC(Vector3d(0.5001, -0.1, -0.1), Vector3d(1.1, 0.0001, 1.1), 
        Vector3d(0.05, 0.025, 0), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    # top-right
    sim.set_DBC(Vector3d(0.5001, 0.9999, -0.1), Vector3d(1.1, 1.1, 1.1), 
        Vector3d(0.05, -0.025, 0), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    
    sim.MDBC_tmin = 0.1
    sim.MDBC_tmax = 0.2
    sim.MDBC_tmin2 = 0.1
    sim.MDBC_tmax2 = 0.2

    sim.gravity = Vector3d(0, 0, 0)
    sim.dt = 0.01
    sim.frame_dt = 0.01
    sim.frame_num = 50
    sim.withCollision = True
    sim.epsv2 = 1e-10

    sim.initialize(sim.cloth_density[0], sim.cloth_Ebase[0], 0.4, \
            sim.cloth_thickness[0], 0)
    Emod = 1.0e3
    nu = 0.3
    rho = 1000.0
    thickness = 1.0e-3
    stiff_mult = 1.0
    h = 1.0e-4
    sim.initialize_rod(rho, Emod, stiff_mult, thickness)
    
    sim.initialize_OIPC(thickness, h)
    # sim.initialize_EIPC(E, nu, thickness, h) # can change offset

    # sim.run()

    sim.write(0)
    # seg_to_curve_Obj(sim.output_folder, "rod"+str(0))
    # sweep_curve(sim.output_folder, "rod"+str(0), 0.5*thickness)
    # sim.write_image(0)
    # do it twice to make sure the image is shown


    
    for f in range(sim.frame_num):
        sim.current_frame = f + 1
        sim.advance_one_frame(f + 1)
        sim.write(f + 1)
        # MeshIO.Convert_Seg_File(sim.output_folder, "rod")
        # seg_to_curve_Obj(sim.output_folder, "rod"+str(f+1))
        # sweep_curve(sim.output_folder, "rod"+str(f+1), thickness)
        if Get_Parameter("Terminate", False):
            break
    




    # self.X, self.segs, self.output_folder + "seg" + str(frame_idx) + ".obj"
    
    # sim.generate_gif()
