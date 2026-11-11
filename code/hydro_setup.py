
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

    def create_disk(self, disk_mass, disk_inner_radius, disk_outer_radius, bh_mass, n_particles, bh_bodies):
        """
        Creates a hydrodynamics disk around the specified location. Stores it as a property of the disk.

        Args:
            disk_location: center of the disk 
            disk_velocity: initial velocity of the disk
            disk_mass: total mass in the disk
            disk_inner_radius: gap between the center and the start of the disk
            disk_outer_radius: furthest radial extent of the disk 
        """
        density_power_exp=1.5
        q_out = 1.

        central_bh = bh_bodies[0]
        infalling_bh = bh_bodies[1]

        central_loc = central_bh.position
        central_vel = central_bh.velocity

        infalling_loc = infalling_bh.position
        infalling_vel = infalling_bh.velocity

        converter1 = nbody_system.nbody_to_si(bh_mass[0], disk_inner_radius[0])
        disk = ProtoPlanetaryDisk(
            n_particles,
            convert_nbody=converter1,
            densitypower=density_power_exp,
            Rmin=1,
            Rmax=disk_outer_radius[0] /disk_inner_radius[0],
            q_out=q_out,
            discfraction=disk_mass[0] / bh_mass[0],
        ).result
        #com = disk.center_of_mass()

        converter2 = nbody_system.nbody_to_si(bh_mass[1], disk_inner_radius[1])
        disk_2 = ProtoPlanetaryDisk(
            n_particles,
            convert_nbody=converter2,
            densitypower=density_power_exp,
            Rmin=1,
            Rmax=disk_outer_radius[1] /disk_inner_radius[1],
            q_out=q_out,
            discfraction=disk_mass[1] / bh_mass[1],
        ).result

        
        disk_2.position += infalling_loc
        disk_2.velocity += infalling_vel
        

        disk.position += central_loc
        disk.velocity += central_vel
        
        disk.add_particles(disk_2)

        self.disk = disk

        bh = Particles()
        bh.mass = bh_mass
        self.bh = bh
        return disk



    def setup_hydro_simulation(self,hydro,particles:List[Particles]):
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
        #self.hydro.dm_particles.add_particles(self.bh)
        self.hydro.dm_particles.add_particles(particles)
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

    def check_high_density(self,threshold):
        """
        this is a primitive version of the simulation, an actual jeans test would be much better.
        """

        highdens = self.disk.select_array(lambda rho:rho>threshold, ["rho"])
        candidate_sinks = highdens.copy()
        self.disk.remove_particles(highdens)

        #newsinks_in_code = hydro.dm_particles.add_particles(candidate_sinks)
        #newsinks = Particles()
        #for nsi in newsinks_in_code:
        #    if nsi not in self.sink_particles:
        #        newsinks.add_particle(nsi)
        #sink_particles.add_sinks(newsinks)
        return highdens


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