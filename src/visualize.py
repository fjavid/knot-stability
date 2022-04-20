from turtle import mode
import numpy as np
import matplotlib.pyplot as plt
import os


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
tre_t_start = [5.63745142487, 8.77904407846, 11.9206367321, 15.0622293856, 18.2038220392, 21.3454146928]
cinque_t_start = [6.87705142487, 10.0186440785, 13.1602367321, 16.3018293856, 19.4434220392, 22.5850146928] 
sep_t_start = [4.83965142487, 7.98124407846, 11.1228367321, 14.2644293856, 17.4060220392, 20.5476146928]
sep_pull_s_id = [int(t_s/sim_dt)+1 for t_s in sep_t_start]
cinque_pull_s_id = [int(t_s/sim_dt)+1 for t_s in cinque_t_start]
tre_pull_s_id = [int(t_s/sim_dt)+1 for t_s in tre_t_start]
# print(sep_pull_s_id)
radius_vals = [0.035, 0.03, 0.025, 0.02, 0.015, 0.01]
pull_range = 0.005
pull_range_id = int(pull_range / pull_velo / sim_dt +1)
print(pull_range_id)
model_name = 'sepfoil'
if model_name == 'trefoil':
    start_id = tre_pull_s_id
elif model_name == 'cinquefoil':
    start_id = cinque_pull_s_id
elif model_name == 'sepfoil':
    start_id = sep_pull_s_id

pull_stretch = []
# Strain values of simulations
for idx, c_r in enumerate(radius_vals):
    stfile = open(base_dir+"/output/pulling_cylin_radius/" + model_name + "_pull_cr" + str(c_r)+"/" + "stretch_rod.txt", 'r')
    stretch_1 = []
    stretch_2 = []
    for line in stfile.readlines():
        stretch_1.append(float(line.split(' ')[0]))
        stretch_2.append(float(line.split(' ')[1]))

    stfile.close() 
    pull_stretch.append(sum(stretch_1[start_id[idx]:start_id[idx]+pull_range_id])/pull_range_id)

print(pull_stretch)

plt.figure(figsize=(11, 4), dpi=200)
plt.plot(radius_vals, pull_stretch, 'b')
plt.tight_layout()
plt.title("rod stretches")
# plt.legend(["dt = 0.01 sec", "dt = 0.001 sec", "dt = 0.0001 sec"])
plt.xlabel("Time step")
plt.ylabel("Stretch")
# plt.ylim(1.0, 1.00005)
# plt.xlim(100, 150)
plt.show() 
# plt.savefig('rod_stretches.png')
# plt.close()