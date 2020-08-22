#Using the parameters from Parameters, this file generates a list of particles formatted such that they can be solved by the ODE solver.

import numpy as np
from Parameters import *
from Helper_Functions import Particle

def circular_orbit_initial_condition(particles, relative_position, orbitted_particle_index, clockwise = True):
    """Returns a Particle object set to circularly orbit a heavy particle. particles should be a list of Particle objects corresponding to all the
    particles in the system (where the first number_of_heavy_masses particles are heavy). relative_position is the relative position of the particle being set
    from the particle it is orbiting, and orbitted_particle_index is the index of the particle being orbitted. clockwise determines if the orbit is clockwise
    or anticlockwise."""
    speed = np.sqrt(G * heavy_mass_masses[orbitted_particle_index] / np.linalg.norm(relative_position))

    relative_velocity = np.cross(relative_position,[0,0,1]) * (speed / np.linalg.norm(relative_position))
    if not(clockwise):
        relative_velocity = relative_velocity * -1
    return Particle(np.add(heavy_mass_starting_positions[orbitted_particle_index],relative_position),np.add(particles[orbitted_particle_index].velocity, relative_velocity))

def circular_polar_to_cartesian(radius, theta):
    """Converts a position from circular polar coordinates to 3d cartesian coordinates. radius is the distance from 0,0,0 and theta is the 
    angle to the x axis."""
    return [radius * np.cos(theta), radius * np.sin(theta), 0]

def calculate_parabolic_orbit_particle_speed(particles, particle_index):
    '''Calculates the speed required for a particle to have a total energy of 0, given the particles list (of Particle objects) and the index of the particle.'''
    gravitational_potential_energy = 0
    for i in range(number_of_heavy_masses):
        if(i != particle_index):
            gravitational_potential_energy -= G *  heavy_mass_masses[i] / np.linalg.norm(np.subtract(particles[particle_index].position, particles[i].position))

    speed = np.sqrt(-1 * 2 * gravitational_potential_energy)

    return speed

def format_particles_for_solving(particles):
    """Reformats the array particles so that every Particle object in the initial array is replaced with three elements corresponding to its position then
    three elements corresponding to its velocity."""
    formatted_particles = []
    for i in range(len(particles)):
        particle = particles[i]
        formatted_particles.extend(particle.position)
        formatted_particles.extend(particle.velocity)

    return formatted_particles

def set_up_initial_conditions():
    """Set up the positions of all heavy masses. Then, using these, calculate the velocity of each heavy mass such that its kinetic energy plus its gravitational potenetial energy is zero.
    This ensures that the total energy is zero, i.e. the orbit is parabolic. Next, generate all light particles, and finally reformat them so they can be sent into
    solve_ode."""

    check_parameters()

    particles = []
    for j in range(number_of_heavy_masses):
        particles.append(Particle(heavy_mass_starting_positions[j], [0,0,0]))

    for j in range(number_of_heavy_masses):
        particles[j].velocity = (heavy_mass_starting_velocity_directions[j]) * calculate_parabolic_orbit_particle_speed(particles, j)

    #Set up light particles
    for k in range(1):
        for i in range(len(light_particle_radii)):
            r = light_particle_radii[i]
            for j in range(points_at_each_radius[i]):     
                    particles.append(circular_orbit_initial_condition(particles, circular_polar_to_cartesian(r, 2*np.pi * j/points_at_each_radius[i]), k, clockwise = (k==1)))

    particles = format_particles_for_solving(particles)

    return particles