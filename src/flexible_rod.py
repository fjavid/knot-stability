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
# base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class FlexibleRod():
    def __init__(self, vertices=[], segments=[], radius=0):
        self.v = vertices
        self.seg = segments
        self.r = radius
        self.input_file = ''
        self.cipc_input = ''
        self.nvert = len(self.v)
        self.nseg = len(self.seg)
        self.timestep = []
        self.u = []
        self.bounds = []
        if self.nvert > 0:
            self.bounds = [np.min(np.array(self.v), axis=0), np.max(np.array(self.v), axis=0)]
    
    def read_from_file(self, file_dir, file_name):
        # self.full_file = file_dir + file_name
        if file_name[-3:] == 'obj':
            self.read_input_obj(file_dir, file_name)

    def read_input_obj(self, file_dir, file_name):
        self.input_file = file_dir + file_name
        make_cipc_file = False
        with open(self.input_file, 'r') as f:
            lines = f.readlines()
        for idx, line in enumerate(lines):
            if line[0] == 'l' or line[0] == 'f':
                make_cipc_file = True
            sline = line.split(' ')
            if sline[0] == 'v':
                self.v.append(np.array([float(sline[1]), float(sline[2]), float(sline[3])], dtype=float))
            if sline[0] == 'l' or sline[0] == 'f' or sline[0] == 's':
                self.seg.append(np.array([int(sline[1])-1, int(sline[2])-1], dtype=int))
        
        self.nvert = len(self.v)
        self.nseg = len(self.seg)
        if make_cipc_file:
            self.cipc_input = self.input_file[:-4] + '_cipc.obj'
            self.write_cipc_input()
        else: 
            self.cipc_input = self.input_file
        
        if self.nvert > 0:
            self.bounds = [np.min(np.array(self.v), axis=0), np.max(np.array(self.v), axis=0)]

    def add_timestep(self, file_dir, file_name, time):
        self.timestep.append(time)
        self.u.append([])
        fullfile = file_dir + file_name
        with open(self.input_file, 'r') as f:
            lines = f.readlines()
        for idx, line in enumerate(lines):
            sline = line.split(' ')
            if sline[0] == 'v':
                self.u[-1].append(np.array([float(sline[1]), float(sline[2]), float(sline[3])], dtype=float))

    def write_cipc_input(self):
        with open(self.cipc_input, 'w') as f:
            for vert in self.v:
                f.write("v %f %f %f\n" % (vert[0], vert[1], vert[2]))
            for seg in self.seg:
                f.write("s %i %i\n" % (seg[0]+1, seg[1]+1))

    def write_curve(self, file_dir, file_name):
        fullname = file_dir + file_name
        if fullname[-4:] != ".obj":
            fullname += '.obj'
        with open(fullname, 'w') as f:
            for vert in self.v:
                f.write("v %f %f %f\n" % (vert[0], vert[1], vert[2]))
            for seg in self.seg:
                f.write("l %i %i\n" % (seg[0]+1, seg[1]+1))

    def make_seg_axes(self):
        axes = []
        for s in self.seg:
            axes.append(self.v[s[1]]-self.v[s[0]])
            axes[-1] /= np.linalg.norm(axes[-1])
        return np.array(axes)

    def make_v_axes(self):
        seg_axes = self.make_seg_axes()
        axes = np.zeros((self.nvert, 3), dtype=float)
        for idx, s in enumerate(self.seg):
            axes[s[0]] += seg_axes[idx]
            axes[s[1]] += seg_axes[idx]
        for idx in range(self.nvert):
            axes[idx] /= np.linalg.norm(axes[idx])
        return axes


    def make_coordinates(srlf, axis, app_n1):
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

    def make_v_aligned_coords(self):
        v_axes = self.make_v_axes()
        v_coords = []
        for idx in range(self.nvert):
            app_n1 = np.array([1.0, 0, 0])
            if idx > 0:
                app_n1 = v_coords[idx-1][1,:]
            v_coords.append(self.make_coordinates(v_axes[idx], app_n1))
        
        return np.array(v_coords)

    def sweep_curve(self, cs_res=8, cap=True):
        self.surf_v = []
        self.surf_face = []
        v_aligned_coords = self.make_v_aligned_coords()

        for v_id, center_v in enumerate(self.v):
            for cs_idx in range(cs_res):
                theta = 2.0 * np.pi * cs_idx /cs_res
                n1 = v_aligned_coords[v_id][1]
                n2 = v_aligned_coords[v_id][2]
                rel_v = self.r * (np.cos(theta)*n1 + np.sin(theta)*n2)
                self.surf_v.append(center_v+rel_v)

        for v_id in range(self.nvert-1):
            idx_s = v_id * cs_res
            for cs_idx in range(cs_res-1):
                self.surf_face.append([idx_s+cs_idx, idx_s+cs_idx+cs_res, idx_s+cs_idx+1])
                self.surf_face.append([idx_s+cs_idx+cs_res, idx_s+cs_idx+cs_res+1, idx_s+cs_idx+1])
            self.surf_face.append([idx_s+cs_res-1, idx_s+cs_res-1+cs_res, idx_s])
            self.surf_face.append([idx_s+cs_res-1+cs_res, idx_s+cs_res, idx_s])
        
        if cap:
            s_idx = len(self.surf_v)
            self.surf_v.append(self.v[0])
            self.surf_v.append(self.v[-1])
            for cs_idx in range(cs_res-1):
                self.surf_face.append([cs_idx, s_idx, cs_idx+1])
            self.surf_face.append([cs_res-1, s_idx, 0])
            for cs_idx in range(cs_res-1):
                self.surf_face.append([s_idx-1-cs_idx, s_idx+1, s_idx-1-cs_idx-1])
            self.surf_face.append([s_idx-1, s_idx+1, s_idx-cs_res])

    def write_surface(self, file_dir, file_name):
        fullname = file_dir + file_name
        if fullname[-4:] != ".obj":
            fullname += '.obj'
        with open(fullname, 'w') as f:
            for vert in self.surf_v:
                f.write("v %f %f %f\n" % (vert[0], vert[1], vert[2]))
            for face in self.surf_face:
                f.write("f %i %i %i\n" % (face[0]+1, face[1]+1, face[2]+1))
            
    def calcul_length(self):
        self.length = 0.0
        for seg in self.seg:
            self.length += np.linalg.norm(self.v[seg[1]]-self.v[seg[0]])
    
        

