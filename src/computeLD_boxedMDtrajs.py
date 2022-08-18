## Processing trajectory data

# Computing LD on (y,py) sections from the trajectory data on the box boundaries. We know the singular values of LD (at maxima because of variable time integration) should show reactive islands.


import numpy as np
import pandas as pd

import solutesolventLJ2dof
import importlib
importlib.reload(solutesolventLJ2dof)

from scipy.integrate import trapz



def compute_ld(datapath, tFlag, filesID, system_mass, bath_mass, time_step, total_int_time):

    ld_sampling = np.array([])
    for index in np.arange(filesID[0],filesID[1],1):
        

        filename = "trj_M" + str(bath_mass) + "T" + str(total_int_time) + "_" + str(tFlag) + \
            "_n" + str(index) + "_0-0.txt"
        
        # boxed_traj = np.loadtxt(datapath + filename)

        boxed_traj = np.array([])
        traj_lengths = np.array([])
        traj_count = 0
        
        
        with open(datapath + filename) as fp:
            for line in fp:
                if line != '\n':
                    boxed_traj = np.append(boxed_traj,pd.to_numeric(line.strip().split()))
                else:
                    boxed_traj = np.reshape(boxed_traj,(int(len(boxed_traj)/4),4),'C')
                    traj_lengths = np.append(traj_lengths, np.size(boxed_traj,0))


                    # compute forward LD for the initial condition at the mid point of the trajectory
                    ke_vec = solutesolventLJ2dof.kinetic_energy(boxed_traj[:,2], boxed_traj[:,3], [system_mass, bath_mass])
                    if tFlag == -1:
                        dt_vec = np.linspace(len(ke_vec)*time_step, 0, len(ke_vec))
                        ld = trapz(ke_vec, dt_vec)
                    elif tFlag == 1:
                        dt_vec = np.linspace(0, len(ke_vec)*time_step, len(ke_vec))
                        ld = trapz(ke_vec, dt_vec)
                    
                    ld_sampling = np.append(ld_sampling, np.append(boxed_traj[0,:], abs(ld)))


                    boxed_traj = np.array([])            
                    traj_count = traj_count + 1
                    
                    # if traj_count > 1e4:
                    #     break

        print(['Finished processing file ' + filename + ' and got %5d'%(traj_count) + ' trajectories'])


    # reshape the post-processed LD data
    ld_sampling = np.reshape(ld_sampling, (int(len(ld_sampling)/5),5))
    # np.savetxt('ld_combined_' + filename, ld_sampling)


    np.savetxt('ld_combined_' + filename, ld_sampling)

if __name__ == '__main__':
    datapath = '/Users/OptimusPrime/Google Drive/systembath2dof-LJrepulsion/Boxed_Trajectories_M0.1/'
    tFlag = 1
    filesID = np.array([0, 2], dtype=int)
    system_mass = 1
    bath_mass = 1
    time_step = 1e-3
    total_int_time = 500
    compute_ld(datapath, tFlag, filesID, system_mass, bath_mass, time_step, total_int_time)
    
    
    