from email.mime import base
import sys
sys.path.insert(0, "/home/Codim-IPC/Python")
sys.path.insert(0, "/home/Codim-IPC/build")
# print(sys.path)
import Drivers
from Drivers.SimulationBase import make_directory
from JGSL import *
import os
import numpy as np
import copy
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print(base_dir)
# sys.path.insert(0, "/home/Codim-IPC/Python")
from flexible_rod import FlexibleRod

PI = np.pi


def translation(vertices, trans):
    v = copy.deepcopy(vertices)
    for idx in range(len(v)):
        v[idx] += trans

    return v

def rotation(vertices, rot_axis, rot_angle):
    v = copy.deepcopy(vertices)
    for idx in range(len(v)):
        v[idx] = np.cos(rot_angle) * v[idx] + np.sin(rot_angle) * np.cross(rot_axis, v[idx])\
                        + (1-np.cos(rot_angle)) * np.dot(rot_axis, v[idx]) * rot_axis

    return v

def align_planar_curve(vertices, idset_1, dir_1, idset_2, dir_2, num_iter=3):
    v = copy.deepcopy(vertices)
    for itr in range(num_iter):
        curr_dir1 = v[idset_1[1]] - v[idset_1[0]]
        curr_dir1 /= np.linalg.norm(curr_dir1)
        rot_axis_1 = -np.cross(curr_dir1, dir_1)
        rot_angle_1 = np.linalg.norm(rot_axis_1)
        rot_axis_1 /= rot_angle_1
        v = rotation(v, rot_axis_1, rot_angle_1)

        curr_dir2 = v[idset_2[1]] - v[idset_2[0]]
        curr_dir2 /= np.linalg.norm(curr_dir2)
        rot_axis_2 = -np.cross(curr_dir2, dir_2)
        rot_angle_2 = np.linalg.norm(rot_axis_2)
        rot_axis_2 /= rot_angle_2
        v = rotation(v, rot_axis_2, rot_angle_2)
    
    return v

def make_kont_trefoil(nv=201, k_xy=0.05, k_z=0.1):
    vertices = []
    segs = []
    for i in range(nv):
        theta = i * 2.0*PI / (nv-1)
        vertices.append(np.array([k_xy * (2.+np.cos(3.0*theta))*np.cos(2.0*theta),
                                    k_xy *(2.+np.cos(3.0*theta))*np.sin(2.0*theta),
                                    k_z*np.sin(3.0*theta)]))

    for i in range(nv-1):
        segs.append(np.array([i, i+1]))
    TOL = 0.02 * (vertices[1][1] - vertices[0][1])
    vertices[0][1] += TOL
    vertices[-1][1] -= TOL

    return vertices, segs

def make_kont_fig8foil(nv=201, k_xy=0.05, k_z=0.1):
    vertices = []
    segs = []
    for i in range(nv):
        theta = i * 2*PI / (nv-1)
        vertices.append(np.array([k_xy * np.cos(3.*theta)*(2.0+np.cos(2.*theta)),
                                    k_xy * np.sin(3.*theta)*(2.0+np.cos(2.*theta)),
                                    k_z*np.sin(4.*theta)]))

    for i in range(nv-1):
        segs.append(np.array([i, i+1]))
    TOL = 0.02 * (vertices[1][1] - vertices[0][1])
    vertices[0][1] += TOL
    vertices[-1][1] -= TOL

    return vertices, segs

def make_kont_cinquefoil(nv=201, k_xy=0.05, k_z=0.1):
    vertices = []
    segs = []
    for i in range(nv):
        theta = i * 2*PI / (nv-1)
        vertices.append(np.array([k_xy * np.cos(2.*theta)*(3.+np.cos(5.*theta)),
                                    k_xy * np.sin(2.*theta)*(3.+np.cos(5.*theta)),
                                    k_z*np.sin(5.*theta)]))

    for i in range(nv-1):
        segs.append(np.array([i, i+1]))
    TOL = 0.02 * (vertices[1][1] - vertices[0][1])
    vertices[0][1] += TOL
    vertices[-1][1] -= TOL

    return vertices, segs

def make_kont_sepfoil(nv=201, k_xy=0.05, k_z=0.1):
    vertices = []
    segs = []
    for i in range(nv):
        theta = i * 2*PI / (nv-1)
        vertices.append(np.array([k_xy * np.cos(2.*theta)*(4. + np.cos(7.*theta)),
                                  k_xy * np.sin(2.*theta)*(4. + np.cos(7.*theta)),
                                  k_z*np.sin(7.*theta)]))

    for i in range(nv-1):
        segs.append(np.array([i, i+1]))
    TOL = 0.02 * (vertices[1][1] - vertices[0][1])
    vertices[0][1] += TOL
    vertices[-1][1] -= TOL

    return vertices, segs

def make_cyliner_face_z(cs_res = 17, len_res = 11, r = 0.05, h = 0.1):
    vertices = []
    faces = []
    for ln in range(len_res):
        z = ln * h / (len_res-1)
        for cs in range(cs_res):
            theta = cs * 2.0*PI / (cs_res-1)
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            vertices.append(np.array([x, y, z], dtype=float))

    for ln in range(len_res-1):
            idx_s = ln * cs_res
            for cs in range(cs_res-1):
                faces.append([idx_s+cs, idx_s+cs+cs_res, idx_s+cs+1])
                faces.append([idx_s+cs+cs_res, idx_s+cs+cs_res+1, idx_s+cs+1])
            faces.append([idx_s+cs_res-1, idx_s+cs_res-1+cs_res, idx_s])
            faces.append([idx_s+cs_res-1+cs_res, idx_s+cs_res, idx_s])

    return vertices, faces

def write_obj_face(vertices, faces, file_dir, file_name):
    fullname = file_dir + file_name
    if fullname[-4:] != ".obj":
        fullname += '.obj'
    with open(fullname, 'w') as f:
        for vert in vertices:
            f.write("v %f %f %f\n" % (vert[0], vert[1], vert[2]))
        for face in faces:
            f.write("f %i %i %i\n" % (face[0]+1, face[1]+1, face[2]+1))

def cut_tail_by_length(vertices, segments, cut_l, tail):
    m_l = 0.0
    nvert = len(vertices)
    segs = []
    if tail == 'beginning':
        for idx in range(1, nvert):
            m_l += np.linalg.norm(vertices[idx]- vertices[idx-1])
            if m_l > cut_l:
                edge_v = idx - 1
                break
        v = vertices[edge_v:]
        for s in segments:
            if s[0] >= edge_v and s[1] >= edge_v:
                segs.append(np.array([s[0]-edge_v, s[1]-edge_v]))
    elif tail == 'end':
        for idx in range(nvert-1, -1, -1):
            m_l += np.linalg.norm(vertices[idx-1]- vertices[idx])
            if m_l > cut_l:
                edge_v = idx
                break
        v = vertices[:edge_v+1]
        for s in segments:
            if s[0] <= edge_v and s[1] <= edge_v:
                segs.append(np.array([s[0], s[1]]))
    return v, segs

def mirror_curve(verts, segms):
    x_0 = verts[0][0]
    v = vertices[1:-1]
    v.extend([np.array([2.*x_0-ov[0], ov[1], -ov[2]]) for ov in vertices[-1:0:-1]])
    s = []
    for idx in range(len(v)-1):
        s.append(np.array([idx, idx+1]))
    s.append(np.array([len(v)-1, 0]))
    
    return v, s
    
    

if __name__ == "__main__":
    nv = 402
    r=0.001
    model_name =  'trefoil' # 'trefoil', 'sepfoil', 'cinquefoil', 'fig8foil'
    # k_xy, and k_z values are set for a total length of l = 0.5
    if model_name == 'trefoil':
        k_xy = 0.0167
        k_z = 0.0106
        vertices, segments = make_kont_trefoil(nv, k_xy, k_z)
    elif model_name == 'cinquefoil':
        k_xy=0.01114
        k_z = 0.005
        vertices, segments = make_kont_cinquefoil(nv, k_xy, k_z)
    elif model_name == 'sepfoil':
        k_xy=0.00813
        k_z = 0.0045
        vertices, segments = make_kont_sepfoil(nv, k_xy, k_z)
    elif model_name == 'fig8foil':
        k_xy=0.0125
        k_z = 0.0067
        vertices, segments = make_kont_fig8foil(nv, k_xy, k_z)
    
    input_dir = base_dir + '/input/' + model_name + '_2throw/'
    try:
        os.mkdir(input_dir)
    except OSError:
        pass
    
    # vertices = rotation(vertices, np.array([0.0, 0.0, 1.0]), -PI/12.0)
    vertices, segments = mirror_curve(vertices, segments)
    # print(segments[200])
    segments.pop(200)
    # print(segments[200])
    # print(vertices[200])
    # print(vertices[201])
    init_config = FlexibleRod(vertices, segments, r)
    init_config.cipc_input = input_dir + model_name + '_2throw.obj'
    init_config.write_cipc_input()
    init_config.write_curve(input_dir, model_name + '2t_paraview.obj')
    init_config.sweep_curve()
    init_config.write_surface(input_dir, model_name + '2t_3D.obj')
    init_config.calcul_length()
    print("%s length : %f" % (model_name, init_config.length))

    sim = Drivers.FEMDiscreteShellBase("double", 3)
    sim.output_folder = base_dir + "/output/" + model_name + '_2throw/' #  + os.path.splitext(os.path.basename(sys.argv[0]))[0] + "/"
    make_directory(sim.output_folder)
    sim.mu = 0.05
    t_total = 5.0
    sim.gravity = Vector3d(0, 0, 0)
    sim.dt = 0.2
    sim.frame_dt = 0.2
    sim.frame_num = int(t_total/sim.frame_dt)
    sim.withCollision = True
    sim.epsv2 = 1e-10

    sim.add_rod_3D(init_config.cipc_input, Vector3d(0.0, 0.0, 0),
            Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, Vector3d(1, 1, 1))
    # z-neg-motion vertex
    sim.set_DBC(Vector3d(0.2, 0.5, 0.47), Vector3d(0.45, 0.51, 0.5), 
        Vector3d(-0.0, 0.0, -0.25), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    # z-pos-motion vertex
    sim.set_DBC(Vector3d(0.2, 0.5, 0.45), Vector3d(0.45, 0.53, 0.48), 
        Vector3d(-0.0, 0.0, 0.25), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)

    sim.MDBC_tmin = 0.0
    sim.MDBC_tmax = 0.2*t_total
    # # sim.MDBC_period = 1.0
    sim.DBCPopBackTStart = 0.5*t_total
    sim.DBCPopBackTEnd = t_total
    sim.DBCPopBackStep = 2
    sim.DBCPopBackAmt = 2
    # # sim.PNTol = 5e-4

    # sim.set_DBC2_with_range(Vector3d(-0.1, -0.1, 0.998), Vector3d(1.0, 1.1, 1.1), 
    #     Vector3d(0.0, 0.0, -10.0*pitch), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, Vector4i(-1, 0, 1000, -1))
    # sim.MDBC_tmin2 = 2.0
    # sim.MDBC_tmax2 = 3.0
    # # sim.MDBC_period2 = 1.0

    sim.initialize(sim.cloth_density[0], sim.cloth_Ebase[0], 0.4, \
            sim.cloth_thickness[0], 0)
    Emod = 1.0e7
    nu = 0.3
    rho = 1000.0
    thickness = 0.5*init_config.r
    stiff_mult = 1.0
    h = 0.05*init_config.r
    sim.initialize_rod(rho, Emod, stiff_mult, thickness)
    
    sim.initialize_OIPC(thickness, h)
    # sim.initialize_EIPC(E, nu, thickness, h) # can change offset
    # sim.run()

    sim.write(0)
    for f in range(sim.frame_num):
        sim.current_frame = f + 1
        sim.advance_one_frame(f + 1)
        sim.write(f + 1)
        if Get_Parameter("Terminate", False):
            break













    # # nv = 401
    # # r=0.001
    # # obj_cylinder_r = 0.05
    # # obj_cylinder_h = 0.1
    # # scaled_cylin_r = 0.035
    # # scaled_cylin_h = 0.1
    # # obj_center = np.array([0, 0, 0.5*scaled_cylin_h])
    # # cylin_center = np.array([0, 1.05*scaled_cylin_r, 0.0])
    # # cylin_move = cylin_center - obj_center
    # # cylin_scale = np.array([scaled_cylin_r/obj_cylinder_r, scaled_cylin_r/obj_cylinder_r, scaled_cylin_h/obj_cylinder_h])
    # # model_name =  'trefoil_pull' # 'trefoil_pull', 'sepfoil_pull', 'cinquefoil_pull', 'fig8foil_pull'
    # # # cyl_v, cyl_f = make_cyliner_face_z(33, 11, 0.05, 0.1)
    # # # write_obj_face(cyl_v, cyl_f , base_dir+'/input/', 'cylinder_v363.obj')
    # # if model_name[:-5] == 'trefoil':
    # #     w_c = 3.506
    # # elif model_name[:-5] == 'cinquefoil':
    # #     w_c = 7.640
    # # elif model_name[:-5] == 'sepfoil':
    # #     w_c = 10 # this is a guess no value is reported for.
    # # init_out_dir = base_dir + "/output/" + model_name[:-5] + '/'
    # # init_config = FlexibleRod()
    # # init_config.read_from_file(init_out_dir, model_name[:-5] + '_final_paraview.obj')
    # # init_config.r = r
    # # init_config.calcul_length()
    # # print("%s length : %f" % (model_name, init_config.length))
    # # l_thread = init_config.length
    # # # length of contact
    # # l_cont = 2. * w_c * np.sqrt(2. * r * scaled_cylin_r)
    # # l_loop = 2. * PI * scaled_cylin_r + l_cont
    # # l_tail = 0.5* (l_thread - l_loop)

    # # init_config.v, init_config.seg = cut_tail_by_length(init_config.v, init_config.seg, l_tail, 'beginning')
    # # init_config.v, init_config.seg = cut_tail_by_length(init_config.v, init_config.seg, l_tail, 'end')
    # # init_config.write_curve(self, file_dir, file_name)
    # # # print(l_tail)
    # # # print(init_config.bounds[0])
    # # # print(init_config.bounds[1])
    # # # neg_x_disp = -(l_tail - np.abs(init_config.bounds[0][0]))
    # # # pos_x_disp = l_tail - init_config.bounds[1][0] 
    # # # print(neg_x_disp)
    # # # print(pos_x_disp)

    # # # init_config.cipc_input = init_out_dir + model_name[:-5] + '_final_paraview.obj'
    # # # init_config.calcul_length()
    






    # # sim = Drivers.FEMDiscreteShellBase("double", 3)
    # # sim.output_folder = base_dir + "/output/" + model_name + '/' #  + os.path.splitext(os.path.basename(sys.argv[0]))[0] + "/"
    # # make_directory(sim.output_folder)
    # # sim.mu = 0.2
    # # t_total = 3.0
    # # sim.gravity = Vector3d(0, 0, 0)
    # # sim.dt = 0.05
    # # sim.frame_dt = 0.1
    # # sim.frame_num = int(t_total/sim.frame_dt)
    # # sim.withCollision = True
    # # sim.epsv2 = 1e-10
    
    # # sim.add_shell_with_scale_3D(base_dir+"/input/cylinder_v363.obj", Vector3d(cylin_move[0], cylin_move[1], cylin_move[2]), Vector3d(cylin_scale[0], cylin_scale[1], cylin_scale[2]), \
    # #     Vector3d(0.0, 0.0, 0), Vector3d(1, 0, 0), 0)
    # # # v, rotCenter, rotAxis, angVelDeg
    # # sim.set_DBC(Vector3d(-0.1, -0.1, -0.1), Vector3d(1.1, 1.1, 1.1), 
    # #     Vector3d(0, 0, 0), Vector3d(0, 0, 0), Vector3d(0, 0, 0), 0)

    # # sim.add_rod_3D(init_config.cipc_input, Vector3d(0.0, 0.0, 0),
    # #         Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, Vector3d(1, 1, 1))
    # # # pos-x vertex
    # # sim.set_DBC(Vector3d(0.999, -0.1, -0.1), Vector3d(1.1, 1.1, 1.1), 
    # #     Vector3d(pos_x_disp, 0.0, 0.0), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    # # # neg-x vertex
    # # sim.set_DBC(Vector3d(-0.1, -0.1, -0.1), Vector3d(0.001, 1.1, 1.1), 
    # #     Vector3d(neg_x_disp, 0.0, 0.0), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)

    # # # sim.MDBC_tmin = 0.0
    # # # sim.MDBC_tmax = 1.0
    # # # # sim.MDBC_period = 1.0
    # # sim.DBCPopBackTStart = 1.0
    # # sim.DBCPopBackTEnd = t_total
    # # sim.DBCPopBackStep = 1
    # # sim.DBCPopBackAmt = 1
    # # # # sim.PNTol = 5e-4

    # # # sim.set_DBC2_with_range(Vector3d(-0.1, -0.1, 0.998), Vector3d(1.0, 1.1, 1.1), 
    # # #     Vector3d(0.0, 0.0, -10.0*pitch), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, Vector4i(-1, 0, 1000, -1))
    # # # sim.MDBC_tmin2 = 2.0
    # # # sim.MDBC_tmax2 = 3.0
    # # # # sim.MDBC_period2 = 1.0

    # # sim.initialize(sim.cloth_density[0], sim.cloth_Ebase[0], 0.4, \
    # #         sim.cloth_thickness[0], 0)
    # # Emod = 1.0e7
    # # nu = 0.3
    # # rho = 1000.0
    # # thickness = 0.5*init_config.r
    # # stiff_mult = 1.0
    # # h = 0.05*init_config.r
    # # sim.initialize_rod(rho, Emod, stiff_mult, thickness)
    
    # # sim.initialize_OIPC(thickness, h)
    # # # sim.initialize_EIPC(E, nu, thickness, h) # can change offset
    # # # sim.run()

    # # sim.write(0)
    # # for f in range(sim.frame_num):
    # #     sim.current_frame = f + 1
    # #     sim.advance_one_frame(f + 1)
    # #     sim.write(f + 1)
    # #     if Get_Parameter("Terminate", False):
    # #         break
