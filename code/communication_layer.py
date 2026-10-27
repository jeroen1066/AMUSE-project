import grav_setup as grav
import hydro_setup as hydro

"""
This file handles the integration of the gravity simulation and the hydro simulation. It is where the check for star formation is done, 
and the simulation is managed. This is also the layer the main program will interact with directly. 
"""

class simulation:
    def __init__(self):
        self.hydro = hydro()
        self.grav = grav()
        pass

    def initial_setup(self,location_bh,vel_bh,masses_bh,masses_disks,inner_gaps,outer_radii):
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

        return

    def propagate_simulation(self):
        """
        iterate the simulations.
        """



if __name__ == "__main__":
    """
    This block of code only runs if this file specifically is ran, not if this code is imported from another module.
    Use this to test whether everything is working properly and test your code independently of the rest.
    """