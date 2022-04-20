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

def align_planar_curve(vertices, idset_1, dir_1, num_iter=30):
    v = copy.deepcopy(vertices)
    for itr in range(num_iter):
        curr_dir1 = v[idset_1[1]] - v[idset_1[0]]
        curr_dir1 /= np.linalg.norm(curr_dir1)
        rot_axis_1 = -np.cross(curr_dir1, dir_1)
        # print(v[idset_1[1]])
        # print(v[idset_1[0]])
        axis_norm = np.linalg.norm(rot_axis_1)
        # print(axis_norm)
        if np.abs(axis_norm) < 0.00001:
            break
        if axis_norm >= 1.0:
            rot_angle_1 = 0.5*PI
        if axis_norm <= -1.0:
            rot_angle_1 = -0.5*PI
        else:
            rot_angle_1 = np.arcsin(axis_norm)
        rot_axis_1 /= axis_norm
        v = rotation(v, rot_axis_1, rot_angle_1)

        # curr_dir2 = v[idset_2[1]] - v[idset_2[0]]
        # curr_dir2 /= np.linalg.norm(curr_dir2)
        # rot_axis_2 = -np.cross(curr_dir2, dir_2)
        # rot_angle_2 = np.linalg.norm(rot_axis_2)
        # rot_axis_2 /= rot_angle_2
        # v = rotation(v, rot_axis_2, rot_angle_2)
    
    return v

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


def cut_tail_by_idx(vertices, segments, cut_idx, tail):
    m_l = 0.0
    nvert = len(vertices)
    segs = []
    if tail == 'beginning':
        # edge_v = 0
        # for idx in range(1, nvert):
        #     m_l += np.linalg.norm(vertices[idx]- vertices[idx-1])
        #     if m_l > cut_l:
        #         edge_v = idx - 1
        #         break
        v = vertices[cut_idx:]
        for s in segments:
            if s[0] >= cut_idx and s[1] >= cut_idx:
                segs.append(np.array([s[0]-cut_idx, s[1]-cut_idx]))
    elif tail == 'end':
        # edge_v = nvert-1
        # for idx in range(nvert-1, -1, -1):
        #     m_l += np.linalg.norm(vertices[idx-1]- vertices[idx])
        #     if m_l > cut_l:
        #         edge_v = idx
        #         break
        v = vertices[:-cut_idx+1]
        for s in segments:
            if s[0] <= nvert-cut_idx and s[1] <= nvert-cut_idx:
                segs.append(np.array([s[0], s[1]]))
    return v, segs

if __name__ == "__main__":
    # nv = 401
    r=0.001
    # obj_cylinder_r = 0.05
    # obj_cylinder_h = 0.1
    # cylin_r_list = [0.035, 0.03, 0.025, 0.02, 0.015, 0.01]
    # friction_list = [0.1, 0.2, 0.3, 0.4, 0.5]
    # # for scaled_cylin_r in cylin_r_list:
    # fric_coeff = 0.5
    # scaled_cylin_r = 0.03
    # scaled_cylin_h = 0.1
    # obj_center = np.array([0, 0, 0.5*scaled_cylin_h])
    # cylin_center = np.array([0, 1.0*scaled_cylin_r + 10.*r, 0.0])
    # cylin_move = cylin_center - obj_center
    # cylin_scale = np.array([scaled_cylin_r/obj_cylinder_r, scaled_cylin_r/obj_cylinder_r, scaled_cylin_h/obj_cylinder_h])
    # model_base_name = 'sepfoil'
    # model_name =  model_base_name+'_pull_fric' + str(fric_coeff) # 'trefoil_pull', 'sepfoil_pull', 'cinquefoil_pull', 'fig8foil_pull'
    # l_thread = 0.5
    # # cyl_v, cyl_f = make_cyliner_face_z(33, 11, 0.05, 0.1)
    # # write_obj_face(cyl_v, cyl_f , base_dir+'/input/', 'cylinder_v363.obj')
    # if model_base_name == 'trefoil':
    #     w_c = 3.506
    # elif model_base_name == 'cinquefoil':
    #     w_c = 7.640
    # elif model_base_name == 'sepfoil':
    #     w_c = 10 # this is a guess no value is reported for.
    # # length of contact
    # l_cont = 2. * w_c * np.sqrt(2. * r * scaled_cylin_r)
    # l_loop = 2. * PI * scaled_cylin_r
    # l_tail = 0.5* (l_thread - l_loop)


    # final_config = FlexibleRod()
    # final_config.r = init_config.r
    # final_config.read_from_file(sim.output_folder, 'rod' + str(sim.frame_num) + '.obj')
    # final_config.calcul_length()
    # print("%s length : %f" % (model_name, final_config.length))


    loop_name = 'grannyfoil'
    loop_ts = 25
    model_name = loop_name + '_pull'
    init_out_dir = base_dir + "/output/" + loop_name + '/'
    init_config = FlexibleRod()
    init_config.read_from_file(init_out_dir,  'rod' + str(loop_ts) + '.obj')
    init_config.r = r

    cut_l = 830
    tail = 'beginning'
    init_config.v, init_config.seg = cut_tail_by_idx(init_config.v, init_config.seg, cut_l, tail)
    cut_l = 500
    tail = 'end'
    init_config.v, init_config.seg = cut_tail_by_idx(init_config.v, init_config.seg, cut_l, tail)

    init_config = FlexibleRod(init_config.v, init_config.seg, init_config.r)
    mid_of_ends = 0.5 * (init_config.v[0] + init_config.v[init_config.nvert-1])
    dist2_midends = np.array([np.linalg.norm(init_config.v[i]-mid_of_ends) for i in range(init_config.nvert)], dtype=float)
    trans_id = np.argmin(dist2_midends)
    init_config.v = translation(init_config.v, -init_config.v[trans_id])

    end2end_dir = np.array([1., 0., 0.], dtype=float)
    end2end_idset = [0, init_config.nvert-1]
    loop_dir = np.array([0.0, 1.0, 0], dtype=float)
    loop_idset = [int(0.5*init_config.nvert), trans_id]
    init_config.v = align_planar_curve(init_config.v, loop_idset, loop_dir, 10)
    init_config.v = align_planar_curve(init_config.v, end2end_idset, end2end_dir, 10)

    mid_of_ends = 0.5 * (init_config.v[0] + init_config.v[init_config.nvert-1])
    dist2_midends = np.array([np.linalg.norm(init_config.v[i]-mid_of_ends) for i in range(init_config.nvert)], dtype=float)
    trans_id = np.argmin(dist2_midends)
    init_config.v = translation(init_config.v, -init_config.v[trans_id])

    init_config.v = align_planar_curve(init_config.v, loop_idset, loop_dir, 5)
    init_config.v = align_planar_curve(init_config.v, end2end_idset, end2end_dir, 5)
    # print(final_config.v[end2end_idset[1]] - final_config.v[end2end_idset[0]])
    input_dir = base_dir + "/input/" + loop_name + '/'
    init_config.write_curve(input_dir, loop_name+'_loop.obj')
    # final_config.sweep_curve()
    # final_config.write_surface(sim.output_folder+'/', model_name+'_final3D.obj')
    init_config.calcul_length()
    print("%s length : %f" % (loop_name, init_config.length))



    obj_cylinder_r = 0.05
    obj_cylinder_h = 0.1
    cylin_r_list = [0.035, 0.03, 0.025, 0.02, 0.015, 0.01]
    friction_list = [0.1, 0.2, 0.3, 0.4, 0.5]
    # for scaled_cylin_r in cylin_r_list:
    fric_coeff = 0.5
    scaled_cylin_r = 0.03
    scaled_cylin_h = 0.2
    obj_center = np.array([0, 0, 0.5*scaled_cylin_h])
    cylin_center = np.array([0, 1.0*scaled_cylin_r + 10.*r, 0.0])
    cylin_move = cylin_center - obj_center
    cylin_scale = np.array([scaled_cylin_r/obj_cylinder_r, scaled_cylin_r/obj_cylinder_r, scaled_cylin_h/obj_cylinder_h])

    # print(l_tail)
    # print(init_config.bounds[0])
    # print(init_config.bounds[1])
    end_velo = 0.02
    t_tight = 5.0
    t_pull = 3.0*0.003 / end_velo
    t_total = t_tight + t_pull
    scale_small = 0.1
    small_shell_r = scale_small * scaled_cylin_r
    # neg_x_disp = -(l_tail - np.abs(init_config.bounds[0][0]))
    # pos_x_disp = l_tail - init_config.bounds[1][0] 
    # tying_t = max(abs(neg_x_disp/end_velo), abs(pos_x_disp/end_velo))
    # pulling_t = 10.0

    # print(neg_x_disp)
    # print(pos_x_disp)

    init_config.cipc_input = input_dir + loop_name + '_cipc.obj'
    init_config.write_cipc_input()
    # init_config.calcul_length()
    # print("%s length : %f" % (model_name, init_config.length))

    sim = Drivers.FEMDiscreteShellBase("double", 3)
    sim.output_folder = base_dir + "/output/granny_pull/" + model_name + '_sr' + str(small_shell_r) + '/' #  + os.path.splitext(os.path.basename(sys.argv[0]))[0] + "/"
    make_directory(sim.output_folder)
    sim.mu = 0.2
    sim.gravity = Vector3d(0, 0, 0)
    sim.dt = 0.05
    sim.frame_dt = 0.1
    sim.frame_num = int(t_total/sim.frame_dt)
    sim.withCollision = True
    sim.epsv2 = 1e-10
    
    sim.add_shell_with_scale_3D(base_dir+"/input/cylinder_v363.obj", Vector3d(cylin_move[0], cylin_move[1], cylin_move[2]), Vector3d(cylin_scale[0], cylin_scale[1], cylin_scale[2]), \
        Vector3d(0.0, 0.0, 0), Vector3d(1, 0, 0), 0)
    # v, rotCenter, rotAxis, angVelDeg
    sim.set_DBC(Vector3d(-0.1, -0.1, -0.1), Vector3d(1.1, 1.1, 1.1), 
        Vector3d(0, 0, 0), Vector3d(0, 0, 0), Vector3d(0, 0, 0), 0)

    sim.add_shell_with_scale_3D(base_dir+"/input/cylinder_v363.obj", Vector3d(cylin_move[0]+0.0, cylin_move[1]-0.06, cylin_move[2]), Vector3d(scale_small*cylin_scale[0], scale_small*cylin_scale[1], cylin_scale[2]), \
        Vector3d(0.0, 0.0, 0), Vector3d(1, 0, 0), 0)
    # v, rotCenter, rotAxis, angVelDeg
    sim.set_DBC(Vector3d(-0.1, -0.1, -0.1), Vector3d(1.1, 1.1, 1.1), 
        Vector3d(0, 0, 0), Vector3d(0, 0, 0), Vector3d(0, 0, 0), 0)

    # print(init_config.cipc_input)
    sim.add_rod_3D(init_config.cipc_input, Vector3d(0.0, 0.0, 0),
            Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, Vector3d(1, 1, 1))
    # pos-x vertex
    sim.set_DBC(Vector3d(0.99, -0.1, -0.1), Vector3d(1.1, 1.1, 1.1), 
        Vector3d(end_velo, 0.0, 0.0), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    # neg-x vertex
    sim.set_DBC(Vector3d(-0.1, -0.1, -0.1), Vector3d(0.01, 1.1, 1.1), 
        Vector3d(-end_velo, 0.0, 0.0), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)

    # sim.MDBC_tmin = 0.0
    # sim.MDBC_tmax = 1.0
    # # sim.MDBC_period = 1.0
    sim.DBCPopBackTStart = t_tight
    sim.DBCPopBackTEnd = t_total
    sim.DBCPopBackStep = 1
    sim.DBCPopBackAmt = 1
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

    print(t_tight)
    print(t_total)
    
    # del sim