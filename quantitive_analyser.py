#This file was made to get qauantative results from the results of a computation of the system.
import numpy as np
import csv
from Parameters import *
from Helper_Functions import Particle

ts =(0,140,199)

def read_csv_file(t):
    """Reads the CSV file and returns the array of particles at time t"""
    with open (save_file_name, newline = '') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')# quotechar = '')
        column_index_at_time_t = -1
        for row in reader:
            for column_index in range(len(row)):
                cell_value = float(row[column_index])
                if (cell_value > t):
                    column_index_at_time_t = column_index
                    break
            break

        if (column_index_at_time_t==-1):
            raise Exception("No results at a time above t were found.")

        particles = []
        for row in reader:
            particles.append(float(row[column_index_at_time_t]))

        return particles

def is_particle_gravitationally_bound(particles, particle_index, heavy_particle_index):
    """Determines if a particle (at particle_index) is gravitationally bound by the heavy particle at 
    index heavy_particle_index. particles should be an array-like of Particle objects. The returned value is a boolean
    corresponding to whether or not the particle is gravitationally bound. """
    relative_position = np.subtract(particles[heavy_particle_index].position, particles[particle_index].position)
    relative_velocity = np.subtract(particles[heavy_particle_index].velocity, particles[particle_index].velocity)

    kinetic_energy = 0.5 * np.linalg.norm(relative_velocity)**2

    gravitational_potential_energy = G * heavy_mass_masses[heavy_particle_index] / (np.linalg.norm(relative_position)**2)

    return (kinetic_energy < gravitational_potential_energy)

def distance(particle_1, particle_2):
    return np.linalg.norm(np.subtract(particle_1.position,particle_2.position))

def print_particle_data(orbited_particle_index):
    for t in ts:
        print("Reseting t")
        print("t =", t)
        all_particles = read_csv_file(t)

        light_particles = []
        heavy_particles = [] 
        for i in range(0,len(all_particles),6):
            if i >= 6 * number_of_heavy_masses:
                light_particles.append(Particle([all_particles[i],all_particles[i+1],all_particles[i+2]],[all_particles[i+3],all_particles[i+4],all_particles[i+5]]))
            else:
                heavy_particles.append(Particle([all_particles[i],all_particles[i+1],all_particles[i+2]],[all_particles[i+3],all_particles[i+4],all_particles[i+5]]))

        for i in range(len(light_particle_radii)):
            initial_radius = light_particle_radii[i]
            start_index = sum(points_at_each_radius[0: i]) 
            for j in range(orbited_particle_index):
                start_index+= sum(points_at_each_radius)
            particles_at_initial_radius = light_particles[start_index:start_index+points_at_each_radius[i]]
            offset_factor = 1.5

            total_number_of_particles = len(particles_at_initial_radius)
            def  is_particle_unpeturbed(particle, orbited_particle_index):
                distance = np.linalg.norm(np.subtract(particle.position,heavy_particles[orbited_particle_index].position))
                if (((distance > initial_radius /offset_factor) and (distance < initial_radius * offset_factor))):
                    return True
                else:
                    return False

            peturbed_particles = [particle for particle in particles_at_initial_radius if is_particle_unpeturbed(particle, orbited_particle_index)]
            number_of_unpeturbed_particles = len(peturbed_particles)
            print("Initial radius: ", initial_radius)
            print("Number unpeturbed: ", number_of_unpeturbed_particles)
            print("Total particles: ", total_number_of_particles)
            print()

print_particle_data(1)