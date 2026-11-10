import grav_setup as grav
import hydro_setup as hydro_sim
from amuse.community.gadget2 import Gadget2
from amuse.units import units
"""
This file handles the integration of the gravity simulation and the hydro simulation. It is where the check for star formation is done, 
and the simulation is managed. This is also the layer the main program will interact with directly. 
"""



class simulation:
    def __init__(self):
        self.hydro_1 = hydro_sim()
        self.hydro_2 = hydro_sim()
        self.grav = grav()
        pass

    def initial_setup(self,location_bh,vel_bh,masses_bh,masses_disks,inner_gaps,outer_radii,particles_per_disk):
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


        self.hydro_1.create_disk(masses_disks[0],inner_gaps[0],outer_radii,masses_bh[0],particles_per_disk)
        self.hydro_2.create_disk(masses_disks[1],inner_gaps[1],outer_radii,masses_bh[1],particles_per_disk)

        self.hydro_1.setup_hydro_simulation(Gadget2)
        self.hydro_2.setup_hydro_simulation(Gadget2)


        return

    def propagate_simulation(self,t_end,timestep):
        """
        iterate the simulations.
        """
        t = 0 |units.yr

        while t < t_end:

            self.hydro_1.propagate_hydro()
            self.hydro_2.propagate_hydro()




if __name__ == "__main__":
    """
    This block of code only runs if this file specifically is ran, not if this code is imported from another module.
    Use this to test whether everything is working properly and test your code independently of the rest.
    """