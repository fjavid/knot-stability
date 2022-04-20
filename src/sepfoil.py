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
print(base_dir)
# sys.path.insert(0, "/home/Codim-IPC/Python")
from flexible_rod import FlexibleRod

PI = np.pi

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

if __name__ == "__main__":
    nv = 201
    k_xy=0.00813
    r=0.001
    k_z = 0.0045

    vertices, segments = make_kont_sepfoil(nv, k_xy, k_z)
    fig8 = FlexibleRod(vertices, segments, r)
    fig8.cipc_input = base_dir + '/sepfoil.obj'
    fig8.write_cipc_input()
    fig8.write_curve(base_dir+'/', 'sepfoil_paraview.obj')
    fig8.sweep_curve()
    fig8.write_surface(base_dir+'/', 'sepfoil3D.obj')
    fig8.calcul_length()
    print("sepfoil length : %f" % fig8.length)

    sim = Drivers.FEMDiscreteShellBase("double", 3)
    sim.output_folder = base_dir+"/output/sepfoil/" #  + os.path.splitext(os.path.basename(sys.argv[0]))[0] + "/"
    make_directory(sim.output_folder)
    sim.mu = 0.1
    
    sim.add_rod_3D(fig8.cipc_input, Vector3d(0.0, 0.0, 0),
            Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, Vector3d(1, 1, 1))
    # bot-z end
    sim.set_DBC(Vector3d(0.99, 0.49, -0.1), Vector3d(1.1, 0.5, 0.5), 
        Vector3d(0.0, 0.0, -0.1), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    # top-z end
    sim.set_DBC(Vector3d(0.99, 0.5, -0.1), Vector3d(1.1, 0.51, 0.5), 
        Vector3d(0.0, 0.0, 0.1), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)

    # sim.MDBC_tmin = 0.0
    sim.MDBC_tmax = 1.0
    # sim.MDBC_period = 1.0
    sim.DBCPopBackTStart = 1.0
    sim.DBCPopBackTEnd = 2.0
    sim.DBCPopBackStep = 2
    sim.DBCPopBackAmt = 2
    # sim.PNTol = 5e-4

    # sim.set_DBC2_with_range(Vector3d(0.99, 0.5, -0.1), Vector3d(1.1, 0.51, 0.5), 
    #     Vector3d(0.0, 0.0, 0.1), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, Vector4i(-1, 0, 1000, -1))
    # sim.MDBC_tmin2 = 2.0
    # sim.MDBC_tmax2 = 3.0
    # # sim.MDBC_period2 = 1.0
    t_total = 4.0
    sim.gravity = Vector3d(0, 0, 0)
    sim.dt = 0.2
    sim.frame_dt = 0.2
    sim.frame_num = int(t_total/sim.frame_dt)
    sim.withCollision = True
    sim.epsv2 = 1e-10

    sim.initialize(sim.cloth_density[0], sim.cloth_Ebase[0], 0.4, \
            sim.cloth_thickness[0], 0)
    Emod = 1.0e7
    nu = 0.3
    rho = 1000.0
    thickness = 0.5*fig8.r
    stiff_mult = 1.0
    h = 0.05*fig8.r
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

    cknot = FlexibleRod()
    cknot.r = r
    cknot.read_from_file(base_dir + '/output/sepfoil/', 'rod20.obj') 
    # cknot.write_cipc_input()
    trans = -0.5 * (cknot.v[0] + cknot.v[-1])
    end_axis = cknot.v[1] - cknot.v[0]
    mid_axis = cknot.v[int(cknot.nvert/2)] + trans
    plane_normal = np.cross(end_axis, mid_axis)
    plane_normal /= np.linalg.norm(plane_normal)
    
    # curr_axis /= np.linalg.norm(curr_axis)
    final_axis = np.array([0, 1.0, 0], dtype=float)
    rot_axis = -np.cross(plane_normal, final_axis)
    rot_angle = np.linalg.norm(rot_axis)
    rot_axis /= rot_angle
    print(rot_axis)
    cknot.trans_rot(trans, rot_axis, rot_angle)

    trans = -0.5 * (cknot.v[0] + cknot.v[-1])
    end_axis = cknot.v[1] - cknot.v[0]
    mid_axis = cknot.v[int(cknot.nvert/2)] + trans
    plane_normal = np.cross(end_axis, mid_axis)
    plane_normal /= np.linalg.norm(plane_normal)
    
    # curr_axis /= np.linalg.norm(curr_axis)
    final_axis = np.array([0, 1.0, 0], dtype=float)
    rot_axis = -np.cross(plane_normal, final_axis)
    rot_angle = np.linalg.norm(rot_axis)
    rot_axis /= rot_angle
    print(rot_axis)
    cknot.trans_rot(trans, rot_axis, rot_angle)


    cknot.write_curve(base_dir+'/', 'sepknot_paraview.obj')
    cknot.sweep_curve()
    cknot.write_surface(base_dir+'/', 'sepknot3D.obj')
    cknot.calcul_length()
    print("sepfoil length : %f" % cknot.length)