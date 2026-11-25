
import numpy as np
from code.communication_layer import *
import amuse
from amuse.lab import *
#%matplotlib inline 
import numpy as np
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.orbital_elements import get_orbital_elements_from_binary
from amuse.plot import scatter 
from amuse.community.ph4.interface import Ph4


"""
initial_conditions
"""

r_p = 2|units.pc                                            #Closest approach between the black holes
r_a = 10000|units.pc                                        #furthest extent of the black holes

masses_bh = np.array([10**6|units.MSun,10**4|units.MSun])   #Central bh first
masses_disks= masses_bh*0.001                               #masses of the disks as a fraction of 
inner_gaps = [1|units.pc,0.01|units.pc]                     #gap on the inside of the disk, between black hole and disk. Too low and it will give arbitrarily small or zero timesteps
outer_radii = [5|units.pc,0.5|units.pc]                     #Furthest extent of the disks
num_particles = 2000                                        #number of particles for the sph dynamics
end_time = 100000|units.yr                                  #time to run the simulation for
timestep = 1000|units.yr                                    #Timesteps for plots to be made.


"""
Running code
"""


sim = simulation(masses_bh[0])
sim.initial_setup(r_a,r_p,masses_bh,masses_disks,inner_gaps,outer_radii,num_particles)
sim.propagate_simulation(end_time,timestep)