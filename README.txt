PARTICLE DYNAMICS SIMULATION PROJECT

MEMBERS: ANTONIO CUADRADO COBO; IXAKA LABADIE GARCIA; ASIER ZULUETA BARBADILLO

ANTONIO CUADRADO COBO: InputData, Init (doInit), Neigbours, Force (calcForce, evalPotential, drawPairDistributionFunction), 
		       Velocity, Output and PDS. BoUML project. Power Point. README

ASIER ZULUETA BARBADILLO: InputData, ResultsFiles (class), R3Point (class) , Output. BoUML project. Power Point.

IXAKA LABADIE GARCIA: Init (velocityDistribution), Force (drawPotential, drawPairDistributionFunction, drawForce) and PDS. BoUML project
		      README


This project consists in the simulation of the movement of the atoms in a fluid starting with the structure of a BCC or FCC type crystal. 
The atoms move under the influence of the Lennard-Jones potential and the Verlet Algorithm has been used to calculate such movement.
The final goal is to obtain the positions of all the atoms of the crystal in each step that the simulation goes through. 
This way, using the program xCrySDen, it is possible to obtain an animation of the movement of the atoms inside the crystal.

**CONTENTS**
This folder contains all the material needed to run the simulation including the Eclipse code 
and the BoUML project to generate that code. The BoUML project is able to generate the complete Eclipse code correctly 
(BoUML version 7.9, Eclipse IDE version 4.13.0; not tested with other versions).
The folder "ws" contains the code that runs the simulation. It also contains another folder, called "Output", in which ".txt" files 
with the energy, positions, velocities... data will be generated.
Furthermore there is another folder called "Input" which contains the file "InputFile" with all the data needed to run the simulation.
In the folder Simulation a runnable "simulation.jar" file of the simulation, the same input example, and an output folder can be found.
Finally, a powerpoint of the presetation we will do is also included. (It will be uploaded apart in the Ixaka Moodle)

**INPUT DATA**
The "InputFile" file includes the data needed for a simulation of the dynamics of Argon as an example. 
It is recommended to use this file as a model for simulations with other elements or different data, 
this way the probability of mistakes in the reading of the file will be minimised. 
The order of the parameters does not affect the reading, but changing names of the data might. It can also be a text file.

**HOW TO RUN THE SIMULATION**
The simulation can be executed by running the "main" operation in the PDS module using Eclipse. Here the path of the input file has to be changed.
The fastest and easiest way is to use the "simulation.jar" file and execute it via the command line, introducing the input file at the same time. 
The results will be created in the "Output" folder.

**********************************************WARNING**********************************************
* When the code from the BoUML model is generated in the Ouput module you will have to delete the *
* import java.awt.*; because if it is not deleted there will be an ambigous for the class List    *
* because the package java.util.*; also import the class List. We have not been able to delete it * 
* directly from the BoUML model.                                                                  *
*************************************************************************************************** 
          
