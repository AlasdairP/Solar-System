"""
 Alasdair Pedley, s1717209
 A class to describe 3D particles
"""

import numpy as np


class Particle3D(object): 
    """
    Class to describe 3D particles.

    Properties:
    position(array) - position in 3D space stored as numpy array [x,y,z]
    velocity(array) - velocity in 3D space stored as numpy array [vx,vy,vz]
    mass(float) - particle mass
    label(string) - particle label

    Instance Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first and second order position updates

    Static Methods:
    * Creating particles when given a text file containing particle information
    * Calculating separation of two particles
    """
    
    def __init__(self,label,x,y,z,vx,vy,vz,mass):
        """
        Initialise a Particle3D instance
        
        :params x,y,z: positions as floats
        :params vx,vy,vz: velocities as floats
        :param mass: mass as float
        :param label: label as string
        """
        self.label = label
        self.position = np.array([x,y,z],float)
        self.velocity = np.array([vx,vy,vz],float)
        self.mass = mass

    def __str__(self):
        """
        Define output format.
        For particle p=(particle1, [2.0,2.0,2.0], [0.5,0.5,0.5], 1.0) this will print as
        "particle1 x y z "
        """
        return("{0:s} {1:f} {2:f} {3:f}\n".format(self.label,self.position[0],self.position[1],self.position[2]))
        
    
    def kinetic_energy(self):
        """
        Return kinetic energy as
        1/2*mass*vel^2
        """
        return 0.5*self.mass*np.dot(self.velocity,self.velocity)      # Velocity is a numpy array so use numpy's dot product to calculate square
        

    # Time integration methods
    def leap_velocity(self,dt,F1):
        """
        First-order velocity update,
        vector_v(t+dt) = vector_v(t) + dt*vector_f(t)/mass

        :param dt: timestep as float
        :param force: 3D force on particle as numpy array
        """
        self.velocity += (dt*F1/self.mass)


    def leap_pos1st(self,dt):
        """
        First-order position update,
        r(t+dt) = r(t) + dt*v(t)    

        :param dt: timestep as float
        """
        self.position += dt*self.velocity


    def leap_pos2nd(self,dt,force):
        """
        Second-order position update,
        r(t+dt) = r(t) + dt*v(t) + 1/2*dt^2*F(t)/mass

        :param dt: timestep as float
        :param force: current force as 3D numpy array
        """
        self.position += dt*self.velocity + 0.5*dt**2*force/self.mass


    @staticmethod
    def particles_from_file(particles_infile):
        
        lines = particles_infile.readlines()  # lines is a list where each line is held as a string element

        particles = []       # Opening an empty list to store particles in, which will be Particle3D instances

        for i in range(len(lines)):           # Loops through all the lines in the file, one particle per line

            line = lines[i]         # The line of interest is the i-th line
            pinfo = line.split(" ") # This splits the line at the spaces, creating a list (of strings) of the elements
            
            label = pinfo[0]
            x = float(pinfo[1])     # Converts strings to floats
            y = float(pinfo[2])
            z = float(pinfo[3])
            vx = float(pinfo[4])
            vy = float(pinfo[5])
            vz = float(pinfo[6])
            mass = float(pinfo[7])
            
            particle = Particle3D(label,x,y,z,vx,vy,vz,mass) # Creates a Particle3D instance

            particles.append(particle)  # Adds the particle to the list
        
        return particles         # "particles" is a list of Particle3D instances 
        
        
    @staticmethod
    def separation(particle1,particle2):
        """
        Returns the separation of two Particle3D objects.
        
        param particle1: The first particle, as a Particle3D instance
        param particle2: The second particle, as a Particle3D instance
        """

        separation = particle2.position - particle1.position
        return separation
        


    




