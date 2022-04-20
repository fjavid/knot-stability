from turtle import mode
import numpy as np
import matplotlib.pyplot as plt
import os

PI = np.pi


base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Figs and plots settings
SMALL_SIZE = 14
BIG_SIZE = 18

plt.rc('font', size=BIG_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIG_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIG_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIG_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)  # fontsize of the figure title

pull_velo = 0.005
sim_dt = 0.25/6.0

r = 0.001
area = PI * r * r
I = 0.25 * PI * r * r * r * r
E = 1.E7
nu = 0.3
rho = 1000.

tre_t_start = [6.87705142487, 10.0186440785, 13.1602367321, 16.3018293856, 19.4434220392, 22.5850146928] 
cinque_t_start = [5.63745142487, 8.77904407846, 11.9206367321, 15.0622293856, 18.2038220392, 21.3454146928]
sep_t_start = [4.83965142487, 7.98124407846, 11.1228367321, 14.2644293856, 17.4060220392, 20.5476146928]
small_cylin_r = [0.0018, 0.0027, 0.0036, 0.0045]
sep_pull_s_id = [int(t_s/sim_dt)+10 for t_s in sep_t_start]
cinque_pull_s_id = [int(t_s/sim_dt)+10 for t_s in cinque_t_start]
tre_pull_s_id = [int(t_s/sim_dt)+10 for t_s in tre_t_start]
# print(sep_pull_s_id)
radius_vals = [0.035, 0.03, 0.025, 0.02, 0.015, 0.01]
pull_range = 0.002
pull_range_id = int(pull_range / pull_velo / sim_dt +1)
print(pull_range_id)
model_name = 'squarefoil'

# if model_name == 'trefoil':
#     start_id = tre_pull_s_id
# elif model_name == 'cinquefoil':
#     start_id = cinque_pull_s_id
# elif model_name == 'sepfoil':
#     start_id = sep_pull_s_id
start_id = 44
end_id = 59

pull_stretch = []
# Strain values of simulations
for idx, c_r in enumerate(small_cylin_r):
    stfile = open(base_dir+"/output/square_pull/" + model_name + "_pull_sr" + str(c_r)+"/" + "stretch_rod.txt", 'r')
    stretch_1 = []
    stretch_2 = []
    for line in stfile.readlines():
        stretch_1.append(float(line.split(' ')[0]))
        stretch_2.append(float(line.split(' ')[1]))

    stfile.close() 
    pull_stretch.append(sum(stretch_1[start_id:end_id])/(end_id - start_id))

# print(pull_range_id, len(stretch_2[start_id[1]:start_id[1]+pull_range_id]))

print(pull_stretch)
# area * r * r / I = 4.0
nondim_force = [4.0 * (strtch - 1.0) for strtch in pull_stretch ]
nondim_r = [cyl_r / r for cyl_r in small_cylin_r]

plt.figure(figsize=(10, 6), dpi=200)
plt.plot(nondim_r[1:], nondim_force[1:], 'bs-')
plt.tight_layout()
# plt.title("Nondimensional Tension")
# plt.legend(["dt = 0.01 sec", "dt = 0.001 sec", "dt = 0.0001 sec"])
plt.xlabel("Nondimensional loop radius")
plt.ylabel("Nondimensional tension")
# plt.ylim(0.0, 0.02)
# plt.xlim(0, 0.5)
# plt.show() 
plt.savefig(model_name+'_pulling.png')
plt.savefig(model_name+'_pulling.eps', format='eps')
# plt.close()