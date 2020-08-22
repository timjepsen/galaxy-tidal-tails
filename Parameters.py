import numpy as np

#This file contains all variable parameters

#G, the gravitational constant
G = 1

#Parameters corresponding to all heavy masses. 
number_of_heavy_masses = 2 #TODO: add a check to make sure this is all correct
heavy_mass_masses = (1,1)
heavy_mass_starting_positions = [[0,0,0],[18,-35,0]]
heavy_mass_starting_velocity_directions= [np.array([0,-1,0]),np.array([0,1,0])]
heavy_mass_radii = [0.5,0.5] #these radii are the radii below which light particles are considered to collide with the heavy particle

light_particle_radii = (2,3,4,5,6,7)
points_at_each_radius = (10,10,500,1000,1000,1000) #(12,18,24,30,36,40)


max_time = 200

save_file_name = "data.csv"

def check_parameters():

    if (number_of_heavy_masses != len(heavy_mass_masses)):
        raise Exception("The number of heavy masses should be the same as the length of the heavy mass masses tuple.")
    elif (number_of_heavy_masses != len(heavy_mass_starting_positions)):
        raise Exception("The number of heavy masses should be the same as the length of the heavy mass starting positions list.")
    elif (number_of_heavy_masses != len(heavy_mass_starting_velocity_directions)):
        raise Exception("The number of heavy masses should be the same as the length of the heavy mass starting velocities list")
    elif (number_of_heavy_masses != len(heavy_mass_radii)):
        raise Exception("The number of heavy masses should be the same as the length of the heavy mass radii list")
    elif (len(light_particle_radii) != len(points_at_each_radius)):
        raise Exception("The length of the radii and points_at_each_radius lists should be the same length.")
    else:
        return True