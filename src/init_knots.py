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
    TOL = 0.02
    vertices[0] += TOL * (vertices[1] - vertices[0])
    vertices[-1] += TOL * (vertices[-2] - vertices[-1])

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
    TOL = 0.02
    vertices[0] += TOL * (vertices[1] - vertices[0])
    vertices[-1] += TOL * (vertices[-2] - vertices[-1])

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
    TOL = 0.02
    vertices[0] += TOL * (vertices[1] - vertices[0])
    vertices[-1] += TOL * (vertices[-2] - vertices[-1])

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
    TOL = 0.02
    vertices[0] += TOL * (vertices[1] - vertices[0])
    vertices[-1] += TOL * (vertices[-2] - vertices[-1])

    return vertices, segs

def make_knot_squarefoil(nv=601, k_xy=0.05, k_z=0.05, open=True):
    vertices = []
    segs = []
    for i in range(nv):
        theta = i * 2*PI / (nv-1)
        vertices.append(np.array([k_xy * (3.*np.sin(theta) + 2.0*np.sin(3.*theta)),
                                    k_xy * (1.*np.cos(theta) - 2.0*np.cos(3.*theta)),
                                    k_z*np.cos(5.*theta)]))

    for i in range(nv-1):
        segs.append(np.array([i, i+1]))
    TOL = 0.02
    vertices[0] += TOL * (vertices[1] - vertices[0])
    vertices[-1] += TOL * (vertices[-2] - vertices[-1])
    if not open:
        segs.append(np.array([nv-1, 0]))
    
    return vertices, segs

def make_knot_grannyfoil(nv=601, k_xy=0.05, k_z=0.05, open=True):
    vertices = []
    segs = []
    for i in range(nv):
        theta = i * 2*PI / (nv-1)
        vertices.append(np.array([k_xy * (3.*np.sin(theta) + 2.0*np.sin(3.*theta)),
                                    k_xy * (1.*np.cos(theta) - 2.0*np.cos(3.*theta)),
                                    k_z*np.sin(10.*theta)]))

    for i in range(nv-1):
        segs.append(np.array([i, i+1]))
    TOL = 0.02
    vertices[0] += TOL * (vertices[1] - vertices[0])
    vertices[-1] += TOL * (vertices[-2] - vertices[-1])
    if not open:
        segs.append(np.array([nv-1, 0]))
    
    return vertices, segs


if __name__ == "__main__":
    nv = 401
    r=0.001
    model_name =  'squarefoil_equi' # 'trefoil', 'cinquefoil', 'sepfoil', 'fig8foil', 'squarefoil', 'grannyfoil'
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
    elif model_name == 'squarefoil_equi':
        k_xy=0.01
        k_z = 0.01
        vertices, segments = make_knot_squarefoil(nv, k_xy, k_z, False)
    elif model_name == 'grannyfoil_equi':
        k_xy=0.01
        k_z = 0.01
        vertices, segments = make_knot_grannyfoil(nv, k_xy, k_z, False)
    
    
    input_dir = base_dir + '/input/' + model_name + '/'
    try:
        os.mkdir(input_dir)
    except OSError:
        pass
    
    init_config = FlexibleRod(vertices, segments, r)
    init_config.cipc_input = input_dir + model_name + '.obj'
    init_config.write_cipc_input()
    init_config.write_curve(input_dir, model_name + '_paraview.obj')
    init_config.sweep_curve()
    init_config.write_surface(input_dir, model_name + '3D.obj')
    init_config.calcul_length()
    print("%s length : %f" % (model_name, init_config.length))

    sim = Drivers.FEMDiscreteShellBase("double", 3)
    sim.output_folder = base_dir + "/output/" + model_name + '/' #  + os.path.splitext(os.path.basename(sys.argv[0]))[0] + "/"
    make_directory(sim.output_folder)
    sim.mu = 0.02
    t_total = 5.0
    sim.gravity = Vector3d(0, 0, 0)
    sim.dt = 0.1
    sim.frame_dt = 0.2
    sim.frame_num = int(t_total/sim.frame_dt)
    sim.withCollision = True
    sim.epsv2 = 1e-10

    sim.add_rod_3D(init_config.cipc_input, Vector3d(0.0, 0.0, 0),
            Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, Vector3d(1, 1, 1))
    # end vertex
    sim.set_DBC(Vector3d(0.999, -0.1, -0.1), Vector3d(1.1, 0.5, 1.1), 
        Vector3d(0.0, 0.0, -0.04), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    # start vertex
    sim.set_DBC(Vector3d(0.999, 0.5, -0.1), Vector3d(1.1, 1.1, 1.1), 
        Vector3d(0.0, 0.0, 0.04), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)

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

    final_config = FlexibleRod()
    final_config.r = init_config.r
    final_config.read_from_file(sim.output_folder, 'rod' + str(sim.frame_num) + '.obj')
    final_config.calcul_length()
    print("%s length : %f" % (model_name, final_config.length))

    mid_of_ends = 0.5 * (final_config.v[0] + final_config.v[final_config.nvert-1])
    dist2_midends = np.array([np.linalg.norm(final_config.v[i]-mid_of_ends) for i in range(final_config.nvert)], dtype=float)
    trans_id = np.argmin(dist2_midends)
    final_config.v = translation(final_config.v, -final_config.v[trans_id])

    end2end_dir = np.array([1., 0., 0.], dtype=float)
    end2end_idset = [0, final_config.nvert-1]
    loop_dir = np.array([0.0, 1.0, 0], dtype=float)
    loop_idset = [int(0.5*final_config.nvert), trans_id]
    final_config.v = align_planar_curve(final_config.v, loop_idset, loop_dir, 10)
    final_config.v = align_planar_curve(final_config.v, end2end_idset, end2end_dir, 10)

    mid_of_ends = 0.5 * (final_config.v[0] + final_config.v[final_config.nvert-1])
    dist2_midends = np.array([np.linalg.norm(final_config.v[i]-mid_of_ends) for i in range(final_config.nvert)], dtype=float)
    trans_id = np.argmin(dist2_midends)
    final_config.v = translation(final_config.v, -final_config.v[trans_id])

    final_config.v = align_planar_curve(final_config.v, loop_idset, loop_dir, 5)
    final_config.v = align_planar_curve(final_config.v, end2end_idset, end2end_dir, 5)
    # print(final_config.v[end2end_idset[1]] - final_config.v[end2end_idset[0]])
    
    final_config.write_curve(sim.output_folder+'/', model_name+'_final_paraview.obj')
    final_config.sweep_curve()
    final_config.write_surface(sim.output_folder+'/', model_name+'_final3D.obj')
    final_config.calcul_length()
    print("%s length : %f" % (model_name, final_config.length))