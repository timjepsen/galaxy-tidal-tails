
def get_particle_attribute(particles, particle_start_index, is_velocity = False):
    """Gets some parameter - either the position (if is_velocity is false) or velocity (if is_velocity is true) - of a particle"""
    velocity_offset = 3 if is_velocity else 0
    answer = []
    for i in range(3):
        answer.append(particles[particle_start_index + i + velocity_offset])
    return answer

def get_particle_position(particles, particle_start_index):
    """Gets the position of a particle at particle_start_index. Requires the array particles to be formatted as an array of particles,
    first position then velocity."""
    return get_particle_attribute(particles, particle_start_index, False)

def get_particle_velocity(particles, particle_start_index):
    """Gets the velocity of a particle at particle_start_index. Requires the array particles to be formatted as an array of particles,
    first position then velocity."""
    return get_particle_attribute(particles, particle_start_index, True)

#The Particle class stores Particles in a tidy, readable way. It can't be used in the ODE solver (which
# requires a 1d array-like of numbers) but is useful for particle generation.
class Particle:
    """The Particle class stores Particles in a tidy, readable way. It can't be used in the ODE solver (which
    requires a 1d array-like of numbers) but is useful for particle generation."""
    def __init__(self, position, velocity):
        self.position = position
        self.velocity = velocity

    def getPosition(self):
        return self.position

    def getVelocity(self):
        return self.velocity