import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import time
import sys

from Parameters import *
from Helper_Functions import get_particle_position, get_particle_velocity
from particles_generator import set_up_initial_conditions
#heavy_mass_starting_velocity_directions= [np.array([0,-1,0]),np.array([0,1,0])]
#heavy_mass_radii = [0.5,0.5]
squares_of_heavy_mass_radii = np.square(heavy_mass_radii)

def gravity_force(t,particles):   
    """Calculate gravitational force on all particles at some point in time t. particles should be an array of all the particles
    in the format [x,y,z, xdot, ydot, zdot]. An array containing the time derivatives of each of these variables is returned."""              
    
    #First, find the positions of all heavy masses
    heavy_mass_positions = []
    for i in range(number_of_heavy_masses):
        heavy_mass_positions.append(get_particle_position(particles,6*i))
    
    answer = np.empty(len(particles))

    for particle_start_index in range(0, len(particles), 6):
        particle_position = get_particle_position(particles, particle_start_index)

        #find the distances to heavy particles
        #fast_subtract replaces np.subtract - it's slightly faster at the expense of being less generalised and elegant.
        def fast_subtract(list_1, list_2):
            return np.array([list_1[0]-list_2[0], list_1[1] - list_2[1], list_1[2] - list_2[2]])
        heavy_particle_distances = [fast_subtract(particle_position, heavy_mass_positions[0]), 
        fast_subtract(particle_position, heavy_mass_positions[1])]
        
        def calculate_acceleration(heavy_particle_distances, particle_start_index):
            """Calculate the total acceleration by adding all acceleration_term values for every heavy particle in the system."""
            def acceleration_term(heavy_particle_distances, heavy_mass_index):
                """Calculate (and return) the acceleration on a particle due to the heavy mass with index heavy_mass_index and at position r (relative to the particle
                whose acceleration is being calculated"""

                #r squared is stored because finding the square root was found to be significantly more computationally expensive
                r_squared = heavy_particle_distances[heavy_mass_index][0]**2 + heavy_particle_distances[heavy_mass_index][1]**2 + heavy_particle_distances[heavy_mass_index][2]**2

                magnitude = (-1 * G * heavy_mass_masses[heavy_mass_index] / r_squared**(3/2))
                if (r_squared < squares_of_heavy_mass_radii[heavy_mass_index]):
                    magnitude = 0
                direction = heavy_particle_distances[heavy_mass_index]

                return (magnitude * direction) #todo: rename the variables in this function
            
            acceleration = 0
            for heavy_mass_index in range(number_of_heavy_masses):  
                if (particle_start_index != 6*heavy_mass_index): #Don't consider the gravitational force on a particle due to itself
                    acceleration += acceleration_term(heavy_particle_distances, heavy_mass_index)
            return acceleration
        #calculate_acceleration = np.vectorize(calculate_acceleration)

        acceleration = calculate_acceleration(heavy_particle_distances, particle_start_index)

        for j in range(3):
            answer[particle_start_index + j] = particles[particle_start_index + j+3] #The time derivative of a particle's position is just its velocity, which is stored three elements along in the array
            answer[particle_start_index + j + 3] = acceleration[j] #The time derivative of velocity is the acceleration, which has been calculated.
    return answer

def heavy_particles_closest_approach_event(t, particles):
    """Returns the dot product of the relative position and velocity of the (first) two heavy particles. This will equal 
    zero at a extremal point of particle distance."""
    relative_position = np.subtract(get_particle_position(particles,6), get_particle_position(particles,0))
    relative_velocity = np.subtract(get_particle_velocity(particles,6), get_particle_velocity(particles, 0))

    return np.dot(relative_position, relative_velocity)

def get_closest_approach_values(integral_answer):
    """Returns the time and distance of closest approach of the (first) two heavy particles. The argument, integral_answer, should be the returned value from a call
    to solve_ivp, where solve_ivp's events arguments some function that crosses 0 at the particle's closest approach."""
    closest_approach_y_events = integral_answer.y_events
    if (len(closest_approach_y_events[0]) == 0):
        raise Exception("find_distance_of_closest_approach failed; no extremal point was reached.")

    elif (len(closest_approach_y_events[0]) > 1):
        raise Exception("find_distance_of_closest_approach failed; multiple extremal points were reached.")
    else:
        closest_approach_time = integral_answer.t_events[0][0]
        particles_at_closest_approach = closest_approach_y_events[0][0]
        closest_approach_distance = np.linalg.norm(np.subtract(get_particle_position(particles_at_closest_approach,6), get_particle_position(particles_at_closest_approach, 0)))#make this line more legible

        return closest_approach_time, closest_approach_distance

def solve_ode(initial_condition):
    """Given an array of particles and velocities initial_condition, this uses SciPy solve_ivp to calculate the motion of the particles.
    The return value from the call to SciPy's solve_ivp is returned."""
    tsToPlot = np.linspace(0,max_time, 1000)
    integral_answer = solve_ivp(gravity_force, [0,max_time], initial_condition, events = heavy_particles_closest_approach_event, method = "RK45", t_eval = tsToPlot, first_step = 0.01)
    if (integral_answer.success == False):
        print("ODE Solver failed")
    try:
        closest_approach_time, closest_approach_distance = get_closest_approach_values(integral_answer)
        print("Closest approach: time  = ",closest_approach_time," and distance =",closest_approach_distance)
    except Exception as err:
        print(err)
    
    return integral_answer

def export_computation_to_csv(ode_solution):
    """Formats the result from solve_ivp, ode_solution, such that it can be saved to a csv file then saves it to a csv file."""
    data = [ode_solution.t,]
    for i in range(len(ode_solution.y)):
        data.append(ode_solution.y[i])
    print("Saving to file, time =", time.time())
    try:
        np.savetxt(save_file_name, data, delimiter=',')
    except:
        input("Saving to file failed. Check that it isn't open.")
        np.savetxt(save_file_name, data, delimiter=',')

def solve_ode_with_timestamps(particles):
    """Runs solve_ode (with argument particles) and prints timestamps at the start and end."""
    start_time = time.time()
    print("Start time: ",start_time)
    answer = solve_ode(particles)
    export_computation_to_csv(answer)
    finish_time = time.time()
    print("Finish time: ",finish_time, ", Difference: ", finish_time - start_time)

    return answer

def run():
    return solve_ode_with_timestamps(set_up_initial_conditions())