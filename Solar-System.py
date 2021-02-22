"""
Alasdair Pedley s1717209 and Daniel Smith s1751483.

This program models the solar system interacting via gravity, using a Velocity Verlet integration.

"""

import sys              
import math
import numpy as np
import matplotlib.pyplot as pyplot
import statistics as stats
from Particle3D import Particle3D

# Constants

G_metres = 6.67408e-11       #m^3 kg^-1 s^-2, SI units.
AU = 1.495978707e11          # 1 AU =  1.49...e11 metres
G = (G_metres/(AU**3))*86400**2  # This gives G in AU^3 kg^-1 days^-2 (86400 seconds in a day)

############## Force and Potential Energy Functions ################

"""
This function calculates the pairwise gravitational force between two particles.
param G: Constant of gravitation in units of AU, kg, days
params particle1, particle2: Particle3D instances

returns: Force on particle 1 as a result of particle 2. (This is easily flipped to give the force on particle 2)
"""

def force_pair(G,particle1,particle2):

    vec_r12 = Particle3D.separation(particle1,particle2)
    r12 = np.linalg.norm(vec_r12)
    r12_hat = vec_r12/r12
    F = (G*particle1.mass*particle2.mass/r12**2)*r12_hat
    return F                         

"""
This function contains two loops. 
The inner loop takes one body, and builds up the total force on that body as a result of all other bodies in the simulation.
The outer loop does this whole process for each of the n bodies in the simulation.
param bodies_list: This is the list of all the bodies, all Particle3D instances.
param n: The number of bodies in the simulation.
returns: List of the forces on each body.
"""

def update_forces(bodies_list,n):
    
    forces_list = []
    for i in range(n):  
        forces_on_i = 0
        for j in range(n):
            if i!=j:
                fi = force_pair(G,bodies_list[i],bodies_list[j])
                forces_on_i += (fi)

        forces_list.append(forces_on_i)

    return forces_list


"""      
This function calculates the potential energy of the pairwise interaction between two bodies.              
param G: Constant of gravitation in units of AU, kg, days.
params particle1, particle2: Particle3D instances.
returns: potential energy    
"""

def pot_energy_pair(G,particle1,particle2):
    
    vec_r12 = Particle3D.separation(particle1,particle2) 
    r12 = np.linalg.norm(vec_r12)
    U = -(G*particle1.mass*particle2.mass)/(r12)
    return U

"""
This function contains two loops. 
The inner loop takes one body, and builds up the total potential energy of that body as a result of all other bodies in the simulation.
The outer loop does this whole process for each of the n bodies in the simulation and adds it to the total.
Note that interactions should not be double counted; the if statement 'if i>j' prevents this.
param bodies_list: This is the list of all the bodies, all Particle3D instances.
param n: The number of bodies in the simulation.
returns: the total potential energy of the system.
"""

def update_pot_energy(bodies_list,n):
    
    pot_energy = 0
    for i in range(n):
        for j in range(n):
            if i>j:
                pot_ij = pot_energy_pair(G,bodies_list[i],bodies_list[j])
                pot_energy += pot_ij

    return pot_energy

def main():
    
    """   
    The following files are read in from from the command line:
    An input file containing particle information.
    An imput file containing the simulation parameters.
    A xyz file to write the trajectory data to for future use in VMD.
    A file to write the post-simulation orbital data to.
    """
    
    if len(sys.argv)!=5:
        
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <input file (particle info)> <input file (simulation parameters)> <output file (trajectory)> <output file (orbit info)>")
        quit()
        
    else:
        input_filename = sys.argv[1]
        params_filename = sys.argv[2]
        trajectory_file = sys.argv[3]
        orbit_file = sys.argv[4]
    

    # Open input files, these variables are the file handles
        
    bodies_infile = open(input_filename,"r") 
    params_infile = open(params_filename,"r")

    params_line = params_infile.read()
    params = params_line.split(" ")   # The params are separated by spaces in the text file; this splits them into a list.

    # Static method of Particle3D is called, which will return a list of particles.
    bodies_list = Particle3D.particles_from_file(bodies_infile) 

    # Simulation parameters
    
    dt = float(params[0])        # This is the timestep, smaller timestep => more accurate.
    t_end = float(params[1])     # This is the total time the simulation runs for.
    numstep = int(t_end/dt)      # Needs to be an integer as you cannot do "half an iteration" in a loop.
    time = 0.0
    
    bodies_infile.close()
    params_infile.close()
    
    n = len(bodies_list)         # Total number of bodies in the simulation.

    #Get initial forces
    
    forces_list = update_forces(bodies_list,n)
    
    #Initialise data lists for plotting later

    time_list = []
    separation_list = []
    
    ######### Centre of mass correction #########

    mass = 0
    for i in range(n):
        mass += bodies_list[i].mass
    
    momentum = 0
    for i in range(n):
        momentum += bodies_list[i].velocity*bodies_list[i].mass
        
    COM_velocity = momentum/mass
    correction = -COM_velocity
    
    for i in range(n):
        bodies_list[i].velocity += correction
     
    ######### End centre of mass correction #########

    # Open output files

    traj_outfile = open("trajectory_file.xyz", "w")
    orbit_outfile = open("orbit_file", "w")
    
    # Create separation list with n-1 body elements (not Sun)

    separ_list = []
    for i in range(n-1): 
        separ_list.append([])     # Needs to be a list of lists

    # Create total energy list for plotting fluctuations

    total_energy_list = []

    ######### Time period calculation set up #########

    """
    The boolean variable "sign" is whether the x component of velocity is +ve or -ve. 
    It will change exactly twice per orbit.
    The time of the changes is recorded, in a list of "special times". 
    Then the time period is found by subtracting every other element in the list.
    This is carried out in the "Time period calculation method" within the integration loop.
    The section here is just setting up the initial lists at time = 0.
    """

    sign_old_list = []    
    for i in range(1,n):
        if bodies_list[i].velocity[0]>0:
            sign_old_list.append(True)
        else:
            sign_old_list.append(False)

    special_t_list = []
    for i in range(n-1):            # n-1 since no Sun
        special_t_list.append([]) 
        
    ############## Start the Verlet time integration loop ##############

    for i in range(numstep): 

        # Update positions
        for i in range(n):
            bodies_list[i].leap_pos2nd(dt,forces_list[i])  
        
        # Note: this relies on forces list having the same number of elements as bodies list

        # Update forces
        forces_list_new = update_forces(bodies_list,n)

        # Update particle velocities by averaging current and new forces
        for i in range(n):
            bodies_list[i].leap_velocity(dt,0.5*(forces_list[i]+forces_list_new[i]))

        # Re-define force value
        forces_list = forces_list_new
        
        # Update total potential energy
        
        pot_energy = update_pot_energy(bodies_list,n)
        
        # Calculate total energy and add to list (used to plot total energy fluctuations)
        
        total_energy = pot_energy
        for i in range(n):
            total_energy += bodies_list[i].kinetic_energy()    
        
        total_energy_list.append(total_energy)
        
        # Add updated position to the Trajectory File
        # The __str__ method from Particle3D is used to output the position in the required format for VMD.

        traj_outfile.write(str(n) + "\n")
        traj_outfile.write("Point = " + str(time) + "\n")
        for i in range(n):
            traj_outfile.write(bodies_list[i].__str__())

        ######### Calculating separations, so that apo- and periapses can be found #########
        
        # Seperation between Sun and bodies ((i-1) so not Moon)  
  
        for i in range(1,(n-1)):
            sep_vec = Particle3D.separation(bodies_list[0],bodies_list[i])
            sep = np.linalg.norm(sep_vec)
            separ_list[i-1].append(sep) 
                    
        # Seperation between Moon (list element 11) and Earth (element 3)

        moon_sep_vec = Particle3D.separation(bodies_list[3],bodies_list[11])
        moon_sep = np.linalg.norm(moon_sep_vec)
        separ_list[10].append(moon_sep)        
        
        ######### Time period calculation method #########
        
        sign_new_list = []
        for i in range(n-1):
            sign_new_list.append(0)

        # Moon velocity needs to be relative to Earth, so this is adjusted, and then switched back after this section.

        bodies_list[11].velocity = bodies_list[11].velocity - bodies_list[3].velocity

        for i in range(1,n):
        
            if bodies_list[i].velocity[0]>0:
                sign_new_list[i-1] = True
            else:
                sign_new_list[i-1] = False

            if sign_new_list[i-1] != sign_old_list[i-1]:
                special_t_list[i-1].append(time)

        # Switching Moon velocity back

        bodies_list[11].velocity = bodies_list[11].velocity + bodies_list[3].velocity

        # Redefine "signs"

        sign_old_list = sign_new_list

        ######### End time period calculation method #########

        # Record and then increase time
        
        time_list.append(time)
        time += dt

    ############### Post-simulation ###############
    """
    Finding the Apoapses(Max), Periapses(Min) and Periods of the orbits.
    For each of these, the data is calculated and then immediately written to the orbit output file.
    """

    for i in range(1,n):

        orbit_outfile.write(bodies_list[i].label)
        orbit_outfile.write("\n")

        max_position = max(separ_list[i-1])
        orbit_outfile.write("Apoapsis of " + bodies_list[i].label + " = " + str(max_position) + " AU ")
        orbit_outfile.write("\n")
        
        min_position = min(separ_list[i-1])
        orbit_outfile.write("Periapses of " + bodies_list[i].label + " = " + str(min_position) + " AU ")
        orbit_outfile.write("\n")
        
        """
        Orbital period calculations.
        Since the x velocity changes exactly twice per rotation, odd and even elements of the list are subtracted.
        The statistics module is used to take the mean of all these individual time periods.
        The error is calculated as the standard error on the mean, which scales with sqrt(N). 
        """

        period_list = []
        for j in range(2,len(special_t_list[i-1])):
            period = special_t_list[i-1][j] - special_t_list[i-1][j-2]
            period_list.append(period)

        period_mean = stats.mean(period_list)
        period_stdev = stats.stdev(period_list)
        period_error = period_stdev / math.sqrt(float(len(period_list)))
        
        orbit_outfile.write("Orbital Period of " + bodies_list[i].label + " = " + str(period_mean) + " earth days  ±" + str(period_error) +" earth days")
        orbit_outfile.write("\n" + "\n")

    ####### Closing the output files ########

    traj_outfile.close()
    orbit_outfile.close()
    
    ####### Plotting the fluctuations of the total energy of the system #######
    pyplot.plot(time_list,total_energy_list)
    pyplot.title("Total energy fluctuations")
    pyplot.xlabel("Time (Days)")
    pyplot.ylabel("Total energy (AU^3 kg^-1 Days^-2)")
    pyplot.show()

main()
