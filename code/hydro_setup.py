
"""
This program handles the hydrodynamics simulation with the disks. It should contain parts to create the disks, propagate the disks
"""

from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particles
from amuse.community.gadget2 import Gadget2
from amuse.ext.protodisk import ProtoPlanetaryDisk

class hydro_sim():
    def __init__(self):
        pass

    def create_disk(self, disk_mass, disk_inner_radius, disk_outer_radius, bh_mass, n_particles):
        """
        Creates a hydrodynamics disk around the specified location. Stores it as a property of the disk.

        Args:
            disk_location: center of the disk 
            disk_velocity: initial velocity of the disk
            disk_mass: total mass in the disk
            disk_inner_radius: gap between the center and the start of the disk
            disk_outer_radius: furthest radial extent of the disk 
        """
        

        converter = nbody_system.nbody_to_si(bh_mass, disk_inner_radius)
        disk = ProtoPlanetaryDisk(
            n_particles,
            convert_nbody=converter,
            densitypower=1.5,
            Rmin=1,
            Rmax=disk_outer_radius /disk_inner_radius,
            q_out=1.0,
            discfraction=disk_mass / bh_mass,
        ).result
        com = disk.center_of_mass()

        self.disk = disk

        bh = Particles()
        bh.mass = bh_mass
        self.bh = bh
        return disk



    def setup_hydro_simulation(self,hydro):
        """
        Set up a simulation 
        """

        
        self.time = 0|units.yr
        self.disk.r2 = self.disk.x**2 +self.disk.y**2
        self.disk = self.disk.sorted_by_attributes("r2")
        radius_max = self.disk.r2.max().sqrt()
        self.hydro = hydro(
            nbody_system.nbody_to_si(self.disk.mass.sum(), radius_max), number_of_workers=4,
                            )

        self.hydro.parameters.epsilon_squared = (10 | units.au) ** 2

        self.hydro.gas_particles.add_particles(self.disk)
        self.hydro.gas_particles.new_channel_to(self.disk)
        self.hydro.dm_particles.add_particles(self.bh)
        #self.hydro.dm_particles.new_channel_to(self.star)

        self.hydro_to_disk = self.hydro.gas_particles.new_channel_to(self.disk)
        self.disk_to_hydro = self.disk.new_channel_to(self.hydro.gas_particles)

        self.hydro.evolve_model(1 | units.s)
        self.hydro_to_disk.copy()

        return

    def propagate_hydro(self,dt):
        """
        todo figure out what the disk propagation should ask and do
        
        """


        self.old_time = self.time
        self.time += dt 

        self.hydro.evolve_model(self.time)
        self.hydro_to_disk.copy()
        return


if __name__ == "__main__":
    """
    This block of code only runs if this file specifically is ran, not if this code is imported from another module.
    Use this to test whether everything is working propersly and test your code independently of the rest.
    """
    sim = hydro_sim()

    sim.create_disk(10|units.MSun, 1|units.au, 40|units.au,10**5 |units.MSun,200)
    sim.setup_hydro_simulation(Gadget2)
    sim.propagate_hydro(200|units.yr,1000)
    print("Done!")



#to use later
def mu(X=None, Y=0.25, Z=0.02, x_ion=0.1):
    """
    Compute the mean molecular weight in kg (the average weight of
    particles in a gas) X, Y, and Z are the mass fractions of
    Hydrogen, of Helium, and of metals, respectively.  x_ion is the
    ionisation fraction (0 < x_ion < 1), 1 means fully ionised.
    """
    if X is None:
        X = 1.0 - Y - Z
    elif abs(X + Y + Z - 1.0) > 1e-6:
        raise Exception(
            "Error in calculating mu: mass " + "fractions do not sum to 1.0"
        )
    return constants.proton_mass / (
        X * (1.0 + x_ion) + Y * (1.0 + 2.0 * x_ion) / 4.0 + Z * x_ion / 2.0
    )