
package pms.initModule;

import pds.common.*;
import java.util.*;
public class Init implements pds.common.IInit {
  /**
   * This is a graph with the histogram of velocities as they are 
   * assigned at the initialization of the lattice. It gets its
   * value at the end of the doInit operation.
   */
  protected pds.common.Graph initialVelocityDistribution;

  /**
   * This static method returns an InitModule object that implements the interface IInit.
   * 
   * @return an object that implements the interface IInit
   */
  public static pds.common.IInit createInit()
  {
    return new Init();
  }

  /**
   * Create the initial state for the given experiment and returns it with the 
   * initial positions, velocities and neighbours lists for all particles.
   * 
   * @param experiment is the experiment object with the data needed to configure and 
   * run the simulation
   * @return a fresh new State ready to start the simulation
   */
  public pds.common.State doInit(pds.common.InputData experiment) {
     //Calculate the size of the unit cell, assuming that the length of the supercell
    	  // is equal to one.
    	  double cell=1.0/(experiment.numberCells);
    	  
    	  //Build the unit cell 
    	  List<Particle> particles=new ArrayList<Particle>();
    	  
    	  switch (experiment.initialLattice) {
    	  
    	  case FCC:
    		  //Construct the FCC lattice from the unit cell 
    		  for (int i=0;i<experiment.numberCells;i++) {
    			  for (int j=0;j<experiment.numberCells;j++) {
    				  for (int z=0;z<experiment.numberCells;z++) {  
    					  particles.add(new Particle(i*cell,j*cell,z*cell));
    					  particles.add(new Particle(0.5*cell+i*cell,0.5*cell+j*cell,z*cell));
    					  particles.add(new Particle(0.5*cell+i*cell,j*cell,0.5*cell+z*cell));
    					  particles.add(new Particle(i*cell,0.5*cell+j*cell,0.5*cell+z*cell));
    					  	
    				  }
    			  }
    		  }
    		  break;
    	  case BCC:
    		//Construct the BCC lattice from the unit cell 
    		  for (int i=0;i<experiment.numberCells;i++) {
    			  for (int j=0;j<experiment.numberCells;j++) {
    				  for (int z=0;z<experiment.numberCells;z++) {  
    					  particles.add(new Particle(i*cell,j*cell,z*cell));
    					  particles.add(new Particle(0.5*cell+i*cell,0.5*cell+j*cell,0.5*cell+z*cell));  	
    				  }
    			  }
    		  }
    		  break;
    	  default: 
    		  System.out.println("Error: Type of laticce not recognised");
    	  }
     
    
    	   //Now, the center will be shifted to be in the center of the supercell
    	  //Due to the size of the supercell is one, the shift can be done easily subtracting 0.5 to the
    	  //positions of all the atoms.
    	  
    	  for (int i=0;i<particles.size();i++) {
    		  particles.get(i).r.x=particles.get(i).r.x-0.5;
    		  particles.get(i).r.y=particles.get(i).r.y-0.5;
    		  particles.get(i).r.z=particles.get(i).r.z-0.5;
    
    	  }
    	  // Now the units will be change, so the fundamental unit of length is made equal to sigma.
    	  // To change the units, multiply it by the real length (number_fcc_units*latcon)
    	  // and divide by sigma.
    	  for (int i=0;i<particles.size();i++) {
    		  particles.get(i).r.x=particles.get(i).r.x*experiment.numberCells*experiment.latcon/experiment.sigma;
    		  particles.get(i).r.y=particles.get(i).r.y*experiment.numberCells*experiment.latcon/experiment.sigma;
    		  particles.get(i).r.z=particles.get(i).r.z*experiment.numberCells*experiment.latcon/experiment.sigma;
    
    	  }
    	  
    	 
    	  
    	  // Next step: set-up the initial velocities
    	  // The initial velocities will be random velocities, with magnitudes conforming 
    	  // to the Maxwell-Boltzman distribution at the required temperature.
    	  // The Maxwell-Boltzmann distribution is a normal distribution with mean = 0 and
    	  // variance = \sigma^{2} = \sqrt{k_{B} T / m} 
    	  // In this program, since the mass of the atoms is taken as one, and
    	  // k_{B} is taken also as 1, then
    	  // variance = \sqrt{T}
    	  
    	  //Initialisation of the random numbers
    	  @SuppressWarnings("unused")
    	  Random first=new Random(22567);
    	  @SuppressWarnings("unused")
    	  Random second=new Random(88);
    	  Random third=new Random(8);
    	  
    	  
    	  for (int i=0;i<particles.size();i++) {
    		  particles.get(i).v.x=Math.sqrt(experiment.temperature)*third.nextGaussian();
    		  particles.get(i).v.y=Math.sqrt(experiment.temperature)*third.nextGaussian();
    		  particles.get(i).v.z=Math.sqrt(experiment.temperature)*third.nextGaussian();
    		 
    	  }
    	  
    	  //The velocities will be corrected, so that there will be no overall momentum.
    	  // To compute the momentum, we profit of the fact that we have
    	  // taken the mass of the atom as a fundamental unit and,
    	  // as a consequence, the velocity of each atom is numerically identical to its
    	  // momentum
    	  double netMomentumSystemX=0;
    	  double netMomentumSystemY=0;
    	  double netMomentumSystemZ=0;
    	  
    	  for (int i=0;i<particles.size();i++) {
    		  netMomentumSystemX=netMomentumSystemX+particles.get(i).v.x;
    		  netMomentumSystemY=netMomentumSystemY+particles.get(i).v.y;
    		  netMomentumSystemZ=netMomentumSystemZ+particles.get(i).v.z; 
    	  }
    	  netMomentumSystemX=netMomentumSystemX/particles.size(); // The size of the ArrayList particles
    	  netMomentumSystemY=netMomentumSystemY/particles.size(); //  is equal to the total number 
    	  netMomentumSystemZ=netMomentumSystemZ/particles.size(); // of atoms in the supercell
    	  
    	  // Correct the velocities so that there is no net momentum
    	  
    	  for (int i=0;i<particles.size();i++) {
    		  particles.get(i).v.x=particles.get(i).v.x-netMomentumSystemX;
    		  particles.get(i).v.y=particles.get(i).v.y-netMomentumSystemY;
    		  particles.get(i).v.z=particles.get(i).v.z-netMomentumSystemZ; 
    	  }
    	  
    	  
    
    	   //In the initial state the potential energy of the atoms will set equal to zero
    	  // and the kinetic energy will be the sum of the cuadratic velocitie in each direction 
    	  // for each particle.
    	  
    	  //Create an object of the class EnergyLog, which holds the kinetic and potential energy of an state
    	 
    	  EnergyLog initialEnergies=new EnergyLog();
    	  initialEnergies.u=0; //potential energy
    	  
    	  //Calculation of the kinetic energy
    	  double kinetic=0;//
    	  
    	  for (int i=0;i<particles.size();i++) {
    		  kinetic=kinetic+particles.get(i).v.kinetic();
    	  }
    	  initialEnergies.k=kinetic;
    	  initialEnergies.total=0;
    	  
    			
    	  State initialState= new State(experiment, particles, initialEnergies);
    	  
    	  
    	  initialVelocityDistribution = velocityDistribution(initialState);
    	  
    	  
    	  
    	  
        return initialState; 
  }

  /**
   * Calculates a Graph showing the statistical distribution of velocities for 
   * all the particles of the system at the given state.
   * 
   * @param s areference to the State object for which the distribution will be calculated
   * @return a reference to the Graph object that holds the velocities distribution function 
   */
  public pds.common.Graph velocityDistribution(pds.common.State s) {
    //The limits of the histogram will be
    	  int left=(int) (-6*Math.sqrt(s.e.temperature));
    	  int right=(int) (6*Math.sqrt(s.e.temperature));
    	  int[] histogram=new int[s.e.histogramResolution];
    	  int casilla;  //Division of the histogram where a given component of the velocity is given
    	  double delta=1.0*(right-left)/s.e.histogramResolution;
    	  for (int i=0;i<s.particles.size();i++) {
    		  casilla=(int)(Math.round((s.particles.get(i).v.x/delta+s.e.histogramResolution/2)));
    		  if (casilla>=s.e.histogramResolution) {
    				casilla=s.e.histogramResolution-1;
    			}
    			if (casilla<0) {
    				casilla=0;
    			}
    		  histogram[casilla]++;		  
    	  }
    	  for (int i=0;i<s.particles.size();i++) {
    		  casilla=(int)(Math.round((s.particles.get(i).v.y/delta+s.e.histogramResolution/2)));
    		  if (casilla>=s.e.histogramResolution) {
    				casilla=s.e.histogramResolution-1;
    			}
    			if (casilla<0) {
    				casilla=0;
    			}
    		  histogram[casilla]++;		  
    	  }
    	  for (int i=0;i<s.particles.size();i++) {
    		  casilla=(int)(Math.round((s.particles.get(i).v.z/delta+s.e.histogramResolution/2)));
    		  if (casilla>=s.e.histogramResolution) {
    				casilla=s.e.histogramResolution-1;
    			}
    			if (casilla<0) {
    				casilla=0;
    			}
    		  histogram[casilla]++;		  
    	  }
    	  
    	  List<R2Point> velocitiesReps=new ArrayList<R2Point>(s.e.histogramResolution);
    	  
    	  for (int i=0;i<s.e.histogramResolution;i++) {
    		  velocitiesReps.add(new R2Point((i-s.e.histogramResolution/2.0)/(1.0*s.e.histogramResolution)*2*right,histogram[i]));
    		  
    		 	
    		  
    	  }
    	  Graph initialMaxwellDistribution=new Graph(velocitiesReps);
    	  
        return initialMaxwellDistribution;
  }

  /**
   * Returns a Graph showing the statistical distribution of velocities for 
   * all the particles of the system as they were assigned at the initial state.
   * 
   * @return a reference to the Graph object that holds the velocities distribution function 
   */
  public pds.common.Graph getInitialVelocityDistribution() {
        return initialVelocityDistribution;
  }

}
