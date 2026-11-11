# we can add GR effects through post Newtonian corrections 
# for this we use a different gravity solver: Hermite_GRX
import numpy as np
import matplotlib.pyplot as plt
from amuse.lab import *
from amuse.ext.orbital_elements import new_binary_from_orbital_elements, generate_binaries
from amuse.community.hermite_grx.interface import Hermitegrx
from amuse.plot import scatter


"""
This file governs the gravity of the gravity part of the simulation. This should create the two black holes, 
"""

class grav():
    def __init__(self, Mcentral):
        """
        This file should create empty containers for point masses for use in the gravity. These will be populated by the following functions and 
        used to create a gravitational simulation in further phases.
        """

        # set the mass of the central black hole
        self.Mass_central_bh = Mcentral

        # empty lists for information to be recorded
        self.times = []
        self.separations = []
        self.relative_velocities = []
        
        self.kinetic_energies = []
        self.potential_energies = []
        self.total_energies = []

    def create_black_hole(self, r_a, r_p, mass_2bh):
        """
        Creates a binary black hole system with properties specified by the arguments:

        Args:
            r_a: the farthest distance in the orbit from the central black hole
            r_p: the closest distance in the orbit from the central black hole
            mass_2bh: the mass of the incoming black hole
        """

        # calculate the orbital elements (semi major axis and eccentricity) from the furthest and closests approach
        semi_major = (r_p + r_a) / 2
        eccentri = (r_a - r_p) / (r_a + r_p)
        
        self.bodies = new_binary_from_orbital_elements(self.Mass_central_bh, mass_2bh, semi_major, eccentri, 240|units.deg)
        self.bodies.name = ["Central black hole", "Incoming black hole"]

        # Calculate Schwarzschild radii based on their mass
        r_s = 2 * constants.G * self.bodies.mass / constants.c**2
        self.bodies.radius = r_s


    def setup_grav_simulation(self):
        # largely taken from book section 4.2.3.3
        converter = nbody_system.nbody_to_si(self.bodies.mass.sum(), 1|units.pc)
    
        gravity = Hermitegrx(converter)
    
        # Include post-Newtonian corrections up to 2.5 PN
        gravity.parameters.perturbation = "2.5PN_EIH"
        gravity.parameters.integrator = "RegularizedHermite"
        gravity.parameters.dt_param = 0.01
        gravity.parameters.light_speed = constants.c
    
        gravity.particles.add_particles(self.bodies)
        self.channel_from_gravity = (gravity.particles.new_channel_to(self.bodies))
        self.gravity = gravity

    def record_information(self):
        """
        Here we record the information that we can obtain during the simulation. None have units yet
        """

        # time
        self.times.append(self.gravity.model_time)
    
        # distance between the black holes
        separation = (self.bodies.position[1] - self.bodies.position[0]).length()
        self.separations.append(separation)
    
        # relative velocity
        relative_velocity = (self.bodies.velocity[1] - self.bodies.velocity[0]).length()
        self.relative_velocities.append(relative_velocity)
    
        # kinetic and potential energy
        self.kinetic_energies.append(self.gravity.kinetic_energy)
        self.potential_energies.append(self.gravity.potential_energy)
    
        # total energy
        total_energy = (self.gravity.kinetic_energy + self.gravity.potential_energy)
        self.total_energies.append(total_energy)

    def propagate_grav_simulation(self, dt, tend):
        """
        Do a part of the gravity simulation
        """
        time = 0 | units.yr
    
        # record the initial conditions
        self.record_information()
        
        time = 0|units.yr
        while self.gravity.model_time < tend:
            time += dt
            self.gravity.evolve_model(time)
            self.channel_from_gravity.copy()
            scatter(self.bodies.x,self.bodies.y,s=np.sqrt(self.bodies.mass/self.bodies.mass.sum()))

            # record current information
            self.record_information()

        plt.show()
        
        self.gravity.stop()
        

    def plot_information(self):
        """
        Function, that when called, returns plots of the recorded information against time.
        """

        # remove units from the lists, matplotlib cant handle them
        times = [t.value_in(units.yr) for t in self.times]
        separations = [r.value_in(units.pc) for r in self.separations]
        velocities = [v.value_in(units.km / units.s)for v in self.relative_velocities]
        kinetic = [e.value_in(units.J) for e in self.kinetic_energies]
        potential = [e.value_in(units.J) for e in self.potential_energies]
        total = [e.value_in(units.J) for e in self.total_energies]

        # 1. Separation
        plt.figure()
        plt.plot(times, separations)
        plt.scatter(times[np.argmin(separations)], separations[np.argmin(separations)], label=f"Closest Separation (t={times[np.argmin(separations)]:.0f} yr) = {separations[np.argmin(separations)]:.4f} pc", c="yellowgreen")
        plt.scatter(times[0], separations[0], label=f"Initial Separation (t=0) = {separations[0]:.1f} pc", c="orange")
        plt.xlabel("Time (yr)")
        plt.ylabel("Separation (pc)")
        plt.title("Black Hole Separation")
        plt.legend()
        plt.show()
    
        # 2. Relative velocity
        plt.figure()
        plt.plot(times, velocities)
        plt.scatter(times[np.argmin(separations)], velocities[np.argmin(separations)], c="yellowgreen", label=f"Rel velocity at closest separation = {velocities[np.argmin(separations)]:.6f} km/s")
        plt.xlabel("Time (yr)")
        plt.ylabel("Relative velocity (km/s)")
        plt.title("Relative Velocity")
        plt.legend()
        plt.show()
    
        # 3. Energies
        plt.figure()
        plt.plot(times, kinetic, label="Kinetic energy")
        plt.plot(times, potential, label="Potential energy")
        plt.plot(times, total, label="Total energy")
        plt.xlabel("Time (yr)")
        plt.ylabel("Energy (J)")
        plt.title("Energy Evolution")
        plt.legend()
        plt.show()



if __name__ == "__main__":
    """
    This block of code only runs if this file specifically is ran, not if this code is imported from another module.
    Use this to test whether everything is working properly and test your code independently of the rest.
    """

    # set the orbital parameters and simulate 
    r_a = 40 | units.pc # farthest distance between incoming and central black hole, is not what actually is simulated
    r_p = 0.6 | units.pc # closests separation
    
    grav_sim = grav(1e7|units.Msun)
    grav_sim.create_black_hole(r_a, r_p, 1e4|units.Msun)

    grav_sim.setup_grav_simulation()
    grav_sim.propagate_grav_simulation(10|units.yr, 40000|units.yr)

    # plot information
    grav_sim.plot_information()

    # NOTE: it's a little bit confusing but the simulation starts from bottom left and ends up top left