# Solar-System
A Spring 2021 project with Daniel Smith to model the Solar System.

Input Files 

	Solar-System.py
	- Main code, with time integration loop & force function that creates the output files

	Particle3D.py
	- File for defining the objects in the simulation (used for positions, velocities etc.)

	bodies.txt
	- list of initial conditions for the bodies in the simulation
	- Starting conditions at 01/01/2000

	simparams.txt 
	- Timestep size ‘dt’ (Days)
	- Total simulation time Duration (Days)
	

Output Files - Names are arbitrary 

	Trajectory Outfile
	- e.g. trajectory.xyz 

	Orbit Outfile
	- e.g. orbit.txt

#########################################################################

To run the simulation, these are the arguments needed in the command line, in order.
-
Argument 1 - Main file

Argument 2 - Input file for the initial conditions of each planet and celestial body (position, velocity & mass)

Argument 3 - Input file for the conditions of the simulation (time-step and duration)

Argument 4 - Output file name for the trajectories of each body (user defined) 

Argument 5 - Output file name for the orbital data of the bodies (user defined)

-
1) To run the simulation, type - python3 Solar-System.py bodies.txt simparams.txt trajectory.xyz orbit.txt - into the terminal whilst in the root folder containing all the files.

2) To visualise the data from the trajectory.xyz file, type: vmd trajectory.xyz - into the terminal (whilst in the root folder). This should open VMD with the trajectory file. The visualisation method by default is not useful therefore going to settings - Select “points” - Set Size ~ 10. This should display the bodies as points that can be seen on screen. To see the orbital motion of the planets sketched out smooth out the motion for the data, which will require at least a full orbit of pluto to see its motion. 

3) The orbit file should be a text file that contains all the information for each planet in the form;

	*Body Name*
	Apihelion - number AU
	Perihelion - number AU
	Orbital Period - number ± Error on mean days

(where AU is Astronomical Units and the orbital period is calculated as the average Period for the duration of the simulation, error on the mean is calculated by the Standard Deviation of the period divided by the square root of the amount of periods in the simulation (For longer periods, less values are found))


Common Errors

 i) If the total duration of the code is less than the period of a full orbit for any of the planets an error will be returned :

##################################################################################
		
		Traceback (most recent call last):
  		File "Solar-System.py", line 350, in <module>
   		 main()
  		File "Solar-System.py", line 331, in main
    			period_mean = stats.mean(period_list)
  		File "/Users/user/opt/anaconda3/lib/python3.7/statistics.py", line 310, in mean
    			raise StatisticsError('mean requires at least one data point')
		statistics.StatisticsError: mean requires at least one data point
	
##################################################################################

To fix this increase the total simulation duration time. 

ii) having a smaller timestep makes the positions of the bodies less accurate. The code may still run but the data will not be accurate. 

    - For timesteps greater then 10 days planets can have unusual orbits.
    - For timesteps greater than 4 days the Moon is no longer in a stable orbit.
    
Suggested simulation parameters:
    
    - dt = 2 days, total time = 200,000 days (given) will give 2 full orbits of Pluto and a stable orbit of the Moon.
    - dt = 4 days, total time = 200,000 days will take half thetime to run  but the orbits (especially the Moon) will be less accurate.
    - dt = 1 day will give a more accurate orbit of the Moon, but will take 4 times longer to run and a larger trajectory file.
    - total time less than 90,000 will not give a full orbit of Pluto, so it will return the error above.
