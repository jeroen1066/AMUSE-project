import grav_setup as grav
import hydro_setup as hydro_sim
from amuse.community.gadget2 import Gadget2
from amuse.units import units
"""
This file handles the integration of the gravity simulation and the hydro simulation. It is where the check for star formation is done, 
and the simulation is managed. This is also the layer the main program will interact with directly. 
"""



class simulation:
    def __init__(self, center_bh_mass):
        self.hydro = hydro_sim()
        self.grav = grav(center_bh_mass)
        pass

    def initial_setup(self,orbelem_bh,location_bh,vel_bh,masses_bh,masses_disks,inner_gaps,outer_radii,particles_per_disk):
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

        self.grav.create_black_holes(orbelem_bh,masses_bh[1])
        
        
        self.hydro.create_disk(masses_disks,inner_gaps,outer_radii,masses_bh,particles_per_disk,location_bh,vel_bh)
        

        particles = self.grav.particles

        self.hydro.setup_hydro_simulation(Gadget2,particles)
        self.grav.setup_grav_simulation()

        self.update_particles_to_disks = self.grav.grav.particles.new_channel_to(self.hydro.hydro.dm_particles)

        return

    def propagate_simulation(self,t_end,timestep):
        """
        iterate the simulations.
        """
        t = 0 |units.yr
        treshold= 100

        while t < t_end:
            t+= timestep
            self.hydro.propagate_hydro(timestep)

            new_stars = self.hydro_1.check_high_density()

            self.grav.particles.add_particles(new_stars)
            self.update_particles_to_disks.update()
            self.propagate_simulation()

        self.grav.gravity.stop()
        self.hydro.hydro.stop()


if __name__ == "__main__":
    """
    This block of code only runs if this file specifically is ran, not if this code is imported from another module.
    Use this to test whether everything is working properly and test your code independently of the rest.
    """