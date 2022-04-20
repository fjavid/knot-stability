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



# helix = FlexibleRod(vertices, segs, r)
# helix.cipc_input = base_dir + '/helix.obj'
# helix.write_cipc_input()
# helix.write_curve(base_dir+'/', 'helix_paraview.obj')
# helix.sweep_curve()
# helix.write_surface(base_dir+'/', 'helix3D.obj')

# handover = FlexibleRod()
# handover.read_from_file(base_dir+'/input/', 'f1x1.obj')
# handover.r = 0.01
# handover.write_curve(base_dir+'/', 'handover_paraview.obj')
# handover.sweep_curve()
# handover.write_surface(base_dir+'/', 'handover3D.obj')

def make_kont_spiral(nv=101, R_0=0.05, r=0.001, turn=1.0, pitch=0.0025, cone=0.0):
    vertices = []
    segs = []
    for i in range(nv):
        theta = i * 2.*PI* turn/(nv-1)
        R = R_0 - cone * (0.5*turn - theta/2./PI)**4
        vertices.append(np.array([R*np.cos(theta), R*np.sin(theta), pitch*theta/2./PI]))
    for i in range(nv-1):
        segs.append(np.array([i, i+1]))
    return vertices, segs

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

if __name__ == "__main__":
    nv = 201
    R_0=0.05
    r=0.001
    turn=2.5
    pitch=0.0025
    cone=0.005
    k_xy = 0.0167
    k_z = 0.0106

    # vertices, segments = make_kont_spiral(nv, R_0, r, turn, pitch, cone)
    vertices, segments = make_kont_trefoil(nv, k_xy, k_z)
    helix = FlexibleRod(vertices, segments, r)
    helix.cipc_input = base_dir + '/trefoil.obj'
    helix.write_cipc_input()
    helix.write_curve(base_dir+'/', 'trefoil_paraview.obj')
    helix.sweep_curve()
    helix.write_surface(base_dir+'/', 'trefoil3D.obj')
    helix.calcul_length()
    print("helix length : %f" % helix.length)

    sim = Drivers.FEMDiscreteShellBase("double", 3)
    sim.output_folder = base_dir+"/output/handover/" #  + os.path.splitext(os.path.basename(sys.argv[0]))[0] + "/"
    make_directory(sim.output_folder)
    sim.mu = 0.1

    # shell =  sim.add_shell_with_scale_3D(base_dir+"/input/cylinder3K.obj", Vector3d(0, -1.15, 0), Vector3d(0.05, 0.5, 0.05), \
    #     Vector3d(0.008, -0.099, 0), Vector3d(1, 0, 0), 90)
    # # v, rotCenter, rotAxis, angVelDeg
    # sim.set_DBC(Vector3d(-0.1, -0.1, -0.1), Vector3d(1.1, 1.1, 1.1), 
    #     Vector3d(0, 0, 0), Vector3d(0, 0, 0), Vector3d(0, 0, 0), 0)
    
    sim.add_rod_3D(helix.cipc_input, Vector3d(0.0, 0.0, 0),
            Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, Vector3d(1, 1, 1))
    # # bot-z end
    # sim.set_DBC(Vector3d(-0.1, -0.1, -0.1), Vector3d(1.0, 1.1, 0.002), 
    #     Vector3d(0.0, 0.0, 2.0*pitch), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    # # top-z end
    # sim.set_DBC(Vector3d(-0.1, -0.1, 0.998), Vector3d(1.0, 1.1, 1.1), 
    #     Vector3d(0.0, 0.0, -2.0*pitch), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    # bot-z end
    sim.set_DBC(Vector3d(0.999, -0.1, -0.1), Vector3d(1.1, 0.5, 1.1), 
        Vector3d(0.0, -0.1, 0.0), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    # top-z end
    sim.set_DBC(Vector3d(0.999, 0.5, -0.1), Vector3d(1.1, 1.1, 1.1), 
        Vector3d(0.0, 0.1, 0.0), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)

    # sim.MDBC_tmin = 0.0
    # sim.MDBC_tmax = 1.0
    # # # sim.MDBC_period = 1.0
    # # sim.DBCPopBackTStart = 0
    # # sim.DBCPopBackTEnd = 1
    # # sim.DBCPopBackStep = 1
    # # sim.DBCPopBackAmt = 1
    # # # sim.PNTol = 5e-4

    # # sim.set_DBC2_with_range(Vector3d(-0.1, -0.1, 0.998), Vector3d(1.0, 1.1, 1.1), 
    # #     Vector3d(0.0, 0.0, -10.0*pitch), Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, Vector4i(-1, 0, 1000, -1))
    # # sim.MDBC_tmin2 = 2.0
    # # sim.MDBC_tmax2 = 3.0
    # # # sim.MDBC_period2 = 1.0
    # t_total = 4.0
    # sim.gravity = Vector3d(0, 0, 0)
    # sim.dt = 0.2
    # sim.frame_dt = 0.2
    # sim.frame_num = int(t_total/sim.frame_dt)
    # sim.withCollision = True
    # sim.epsv2 = 1e-10

    # sim.initialize(sim.cloth_density[0], sim.cloth_Ebase[0], 0.4, \
    #         sim.cloth_thickness[0], 0)
    # Emod = 1.0e7
    # nu = 0.3
    # rho = 1000.0
    # thickness = 0.5*helix.r
    # stiff_mult = 1.0
    # h = 0.05*helix.r
    # sim.initialize_rod(rho, Emod, stiff_mult, thickness)
    
    # sim.initialize_OIPC(thickness, h)
    # # sim.initialize_EIPC(E, nu, thickness, h) # can change offset
    # # sim.run()

    # sim.write(0)
    # for f in range(sim.frame_num):
    #     sim.current_frame = f + 1
    #     sim.advance_one_frame(f + 1)
    #     sim.write(f + 1)
    #     if Get_Parameter("Terminate", False):
    #         break
