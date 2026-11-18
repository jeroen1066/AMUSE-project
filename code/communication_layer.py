from grav_setup import grav
from hydro_setup import hydro_sim
from amuse.community.gadget2 import Gadget2
from amuse.community.phantom import Phantom
from amuse.units import units
import numpy as np
"""
This file handles the integration of the gravity simulation and the hydro simulation. It is where the check for star formation is done, 
and the simulation is managed. This is also the layer the main program will interact with directly. 
"""



class simulation:
    def __init__(self, center_bh_mass):
        self.hydro = hydro_sim()
        self.grav = grav(center_bh_mass)
        pass

    def initial_setup(self,apoapsis,periapsis,masses_bh,masses_disks,inner_gaps,outer_radii,particles_per_disk):
        """
        This creates the black holes and disks around them from the gravity and hydro codes.

        Args:
            location_bh: location of the second black hole. The first black hole is always generated at the origin.
            vel_bh: the initial velocity of the second black hole. The first black hole is always at rest.
            masses_bh: tuple of the masses of the black holes
            masses_disks: tuple of the masses of the two disks
            inner_gaps: tuple of the gap before the inner disks starts
            outer_radii: tuple of the outer radii of the disks
        """

        self.grav.create_black_hole(apoapsis,periapsis,masses_bh[1])
        
        self.hydro.create_disk(masses_disks,inner_gaps,outer_radii,masses_bh,particles_per_disk,self.grav.bodies)
        

        particles = self.grav.bodies

        self.hydro.setup_hydro_simulation(Gadget2,particles)
        self.grav.setup_grav_simulation()

        self.update_particles_to_disks = self.grav.gravity.particles.new_channel_to(self.hydro.hydro.dm_particles)

        return

    def propagate_simulation(self,t_end,timestep):
        """
        iterate the simulations.
        """
        t = 0 |units.yr
        treshold= 100
        filename = "disk_overview"
        index = 0

        while t < t_end:
            t+= timestep
            self.hydro.propagate_hydro(timestep)

            #new_stars = self.hydro.check_high_density(1000 | units.kg/units.m**3)
            new_stars = self.hydro.check_jeans_mass()

            self.grav.bodies.add_particles(new_stars)
            self.update_particles_to_disks.copy()
            self.grav.propagate_grav_simulation_single_step(timestep)
            self.hydro.plot_disks(filename,index)
            index += 1
            print(t)

        self.grav.gravity.stop()
        self.hydro.hydro.stop()
        


if __name__ == "__main__":
    """
    This block of code only runs if this file specifically is ran, not if this code is imported from another module.
    Use this to test whether everything is working properly and test your code independently of the rest.
    """

    sim = simulation(10**6|units.MSun)
    r_p = 2|units.pc
    r_a = 10000|units.pc

    masses_bh = np.array([10**6|units.MSun,10**4|units.MSun])
    masses_disks= masses_bh*0.0001
    inner_gaps = [1|units.pc,0.1|units.pc]
    outer_radii = [15|units.pc,0.5|units.pc]

    sim.initial_setup(r_a,r_p,masses_bh,masses_disks,inner_gaps,outer_radii,20000)
    sim.propagate_simulation(1000000|units.yr,1000|units.yr)