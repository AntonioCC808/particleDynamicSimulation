
package pds.forceModule;

import pds.common.*;
import java.util.*;
/**
 * This class calculates the current accelerations for all particles.
 */
public class Force implements pds.common.IForce {
  /**
   * The state for which forces and hence accelerations are calculated.
   * This value is initialized by the constructor.
   */
  protected pds.common.State s;

  /**
   * This static method creates an object of a class that implements the IForce interface
   * 
   * @param s the reference to the state object that holds the system for which the 
   * acceleration and potentials will be calculated.
   * @return an object of a class that implements the IForce interface
   */
  public static pds.common.IForce createForce(pds.common.State s)
  {
    	  return new Force(s);
  }

  /**
   * This constructor receives the state for which the Force object
   * will be created.
   * 
   * @param s the state for whose particles the accelerations will be 
   * calculated when requested
   */
  public  Force(pds.common.State s) {
      this.s = s;
  }

  /**
   * Calculates the forces that apply to all particles and update the 
   * corresponding accelerations in the current State. The operation 
   * returns the current potential energy of the system. 
   * The positions of every particle are taken from State, as well as the 
   * neighbour list of each particle. For a particle, just the particles 
   * stated in the neighbour list are taken into account.
   * The force acting over a particle because of the other ones has three components. 
   * First, it calculates the distance between each neighbour particle and 
   * the force is calculated using the equation for every particle against this one,
   * therefore the total force in a component is the sum of all these pair forces.
   * The force acting on a particle operates in the same way for the 3 components.
   * 
   * @return a double real value with the current potential energy of the system.
   */
  public double calcForce() {
    	// Set up the forces equal to zero in each iteration of the simulation to
    		// calculate the new ones
    		for (Particle p : s.particles) {
    			p.a.x = 0;
    			p.a.y = 0;
    			p.a.z = 0;
    		}
    		// In reduced units, epsilon is taken as 1.0_dp,
    		// Length of the supercell (in units of sigma)
    		double boxSize = s.e.numberCells * s.e.latcon / s.e.sigma;
    		// Initialize the potential for the current coordinates
    		double potential = 0;
    		// TEST
    		double maxAcelX = 0;
    		double maxAcelY = 0;
    		double maxAcelZ = 0;
    		double newAcelerationX;
    		double newAcelerationY;
    		double newAcelerationZ;
    
    		for (Particle p : s.particles) {
    			if (p.numberNeighbors > 0) {
    				for (Integer jneig : p.vicinity) {
    
    					// Find the relative position between the two atoms
    					double relDisX = p.r.x - s.particles.get(jneig).r.x;
    					double relDisY = p.r.y - s.particles.get(jneig).r.y;
    					double relDisZ = p.r.z - s.particles.get(jneig).r.z;
    					// Consider the periodic boundary conditions. The calculation of minimum
    					// image distances is simplified by the use of reduced units:
    					// the length of the sigma parameter is taken to define the
    					// fundamental unit of length in the simulation.
    					relDisX = relDisX - boxSize * (Math.round(relDisX / boxSize));
    					relDisY = relDisY - boxSize * (Math.round(relDisY / boxSize));
    					relDisZ = relDisZ - boxSize * (Math.round(relDisZ / boxSize));
    					// Distance module
    					double modRelDistance = Math
    							.sqrt(Math.pow(relDisX, 2) + Math.pow(relDisY, 2) + Math.pow(relDisZ, 2));
    
    					// Calculation of the potential
    					potential = potential + evalPotential(modRelDistance);
    					
    					// Value of the derivative of the force at the cutOff Radius
    					double dcutOffPotential = 24 / s.e.cutoffRadius
    							* (2 * Math.pow(sigmaReduced / s.e.cutoffRadius, 12)
    									- Math.pow(sigmaReduced / s.e.cutoffRadius, 6));
    					
    					// Calculation of the force for each particle
    					switch (s.e.potential) {
    					case LENNARDJONESUNSHIFTED:
    						
    						if (modRelDistance < s.e.cutoffRadius) {
    							// Acceleration/Force X
    							newAcelerationX = (24 / Math.pow(modRelDistance, 2)
    									* (2 * Math.pow(sigmaReduced / modRelDistance, 12)
    											- Math.pow(sigmaReduced / modRelDistance, 6))) * relDisX;
    
    							p.a.x = p.a.x + newAcelerationX;
    
    							// Acceleration/Force Y
    							newAcelerationY = (24 / Math.pow(modRelDistance, 2)
    									* (2 * Math.pow(sigmaReduced / modRelDistance, 12)
    											- Math.pow(sigmaReduced / modRelDistance, 6))) * relDisY;
    
    							p.a.y = p.a.y + newAcelerationY;
    
    							// Acceleration/Force Z
    							newAcelerationZ = (24 / Math.pow(modRelDistance, 2)* (2 * Math.pow(sigmaReduced / modRelDistance, 12)- 
    									Math.pow(sigmaReduced / modRelDistance, 6))) * relDisZ;
    
    							p.a.z = p.a.z + newAcelerationZ;
    
    							///// TEST////////////
    							if (maxAcelX < Math.abs(p.a.x)) {
    								maxAcelX = p.a.x;
    							}
    							if (maxAcelY < Math.abs(p.a.y)) {
    								maxAcelY = p.a.y;
    							}
    							if (maxAcelZ < Math.abs(p.a.z)) {
    								maxAcelZ = p.a.z;
    							}
    							//////////////////////////
    							// Applying the third Newton law, the force that the atom jneig acts on atom p
    							// is minus the force that atom p acts on atom jneig.
    							s.particles.get(jneig).a.x = s.particles.get(jneig).a.x - newAcelerationX;
    							s.particles.get(jneig).a.y = s.particles.get(jneig).a.y - newAcelerationY;
    							s.particles.get(jneig).a.z = s.particles.get(jneig).a.z - newAcelerationZ;
    
    						}
    						break;
    					case LENNARDJONESSHIFTED:
    						if (modRelDistance < s.e.cutoffRadius) {
    							// Acceleration/Force X
    							newAcelerationX = (24 / Math.pow(modRelDistance, 2)
    									* (2 * Math.pow(sigmaReduced / modRelDistance, 12)
    											- Math.pow(sigmaReduced / modRelDistance, 6))) * relDisX;
    
    							p.a.x = p.a.x + newAcelerationX;
    
    							// Acceleration/Force Y
    							newAcelerationY = (24 / Math.pow(modRelDistance, 2)
    									* (2 * Math.pow(sigmaReduced / modRelDistance, 12)
    											- Math.pow(sigmaReduced / modRelDistance, 6))) * relDisY;
    
    							p.a.y = p.a.y + newAcelerationY;
    
    							// Acceleration/Force Z
    							newAcelerationZ = (24 / Math.pow(modRelDistance, 2)* (2 * Math.pow(sigmaReduced / modRelDistance, 12)- 
    									Math.pow(sigmaReduced / modRelDistance, 6))) * relDisZ;
    
    							p.a.z = p.a.z + newAcelerationZ;
    
    							// Applying the third Newton law, the force that the atom jneig acts on atom p
    							// is minus the force that atom p acts on atom jneig.
    							s.particles.get(jneig).a.x = s.particles.get(jneig).a.x - newAcelerationX;
    							s.particles.get(jneig).a.y = s.particles.get(jneig).a.y - newAcelerationY;
    							s.particles.get(jneig).a.z = s.particles.get(jneig).a.z - newAcelerationZ;
    
    						}
    						break;
    					case LENNARDJONESTOTALSHIFTED:
    						// Value of the derivative of the force at the cutOff Radius
    						if (modRelDistance < s.e.cutoffRadius) {
    							// Acceleration/Force X
    							newAcelerationX = (24 / Math.pow(modRelDistance, 2)
    									* (2 * Math.pow(sigmaReduced / modRelDistance, 12)
    											- Math.pow(sigmaReduced / modRelDistance, 6))
    									- dcutOffPotential / modRelDistance) * relDisX;
    
    							p.a.x = p.a.x + newAcelerationX;
    
    							// Acceleration/Force Y
    							newAcelerationY = (24 / Math.pow(modRelDistance, 2)
    									* (2 * Math.pow(sigmaReduced / modRelDistance, 12)
    											- Math.pow(sigmaReduced / modRelDistance, 6))
    									- dcutOffPotential / modRelDistance) * relDisY;
    
    							p.a.y = p.a.y + newAcelerationY;
    
    							// Acceleration/Force Z
    							newAcelerationZ = (24 / Math.pow(modRelDistance, 2)
    									* (2 * Math.pow(sigmaReduced / modRelDistance, 12)
    											- Math.pow(sigmaReduced / modRelDistance, 6))
    									- dcutOffPotential / modRelDistance) * relDisZ;
    
    							p.a.z = p.a.z + newAcelerationZ;
    							
    							// Applying the third Newton law, the force that the atom jneig acts on atom p
    							// is minus the force that atom p acts on atom jneig.
    							s.particles.get(jneig).a.x = s.particles.get(jneig).a.x - newAcelerationX;
    							s.particles.get(jneig).a.y = s.particles.get(jneig).a.y - newAcelerationY;
    							s.particles.get(jneig).a.z = s.particles.get(jneig).a.z - newAcelerationZ;
    
    						}
    						break;
    					default:
    						System.out.println("Error: Type of potential not recognised");
    					}
    
    				}
    			}
    		}
    
    		return potential;
  }

  /**
   * This function returrns the potential between two particles of the 
   * current experiment evaluated for a concrete distance given as input
   * parameter.
   * 
   * @param distance value in meters that indicates the distance between 
   * the particles for which the potential is to be calculated
   * @return value of the potential calculated for two particles of the current 
   * experiment when they are at the distance given as invocation argument.  
   */
  public double evalPotential(double distance) {
    	double pot = 0;
    		double cutOffPotential = 4
    				* (Math.pow(sigmaReduced / s.e.cutoffRadius, 12) - Math.pow(sigmaReduced / s.e.cutoffRadius, 6));
    		double dcutOffPotential = 24 / s.e.cutoffRadius
    				* (2 * Math.pow(sigmaReduced / s.e.cutoffRadius, 12) - Math.pow(sigmaReduced / s.e.cutoffRadius, 6));
    		switch (s.e.potential) {
    		case LENNARDJONESUNSHIFTED:
    
    			if (distance < s.e.cutoffRadius) {
    				pot = 4 * (Math.pow(sigmaReduced / distance, 12) - Math.pow(sigmaReduced / distance, 6));
    			}
    			break;
    		case LENNARDJONESSHIFTED:
    			if (distance < s.e.cutoffRadius) {
    				pot = 4 * (Math.pow(sigmaReduced / distance, 12) - Math.pow(sigmaReduced / distance, 6))- cutOffPotential ;
    			}
    
    			break;
    		case LENNARDJONESTOTALSHIFTED:
    			if (distance < s.e.cutoffRadius) {
    				pot = 4 * (Math.pow(sigmaReduced / distance, 12) - Math.pow(sigmaReduced / distance, 6))
    						- cutOffPotential - dcutOffPotential * (distance - s.e.cutoffRadius);
    			}
    
    			break;
    		default:
    			System.out.println("Error: Type of potential not recognised");
    		}
    		return pot;
  }

  /**
   * This method returns a graph in the form of a list of points (x,y) that 
   * describe the potential for the range of distance vales observed in the 
   * current experiment.
   * The number of points is taken from the method numberOfPoints that comes from Input.
   * The objective of this method is to provide (review) graph which contains two list of numbers
   * In the x coordinate we take as a first value, the minimum distance mind and as a maximum
   * value the maximum distance maxd. The intermediate points are values separate (maxd-mind)/numberOfPoints
   * In the Y coordinate, the potential for each value of the distances is calculated using the method evalPotential
   * 
   * @return a graph that is a list of pairs of numbers (x,y). These values 
   * are of the type double. x represents the distance in meters and y potentials
   * in energy units. The first and last points in the graph correspond to the minimum and 
   * maximum distances observed among all particles in the current experiment.
   */
  public pds.common.Graph drawPotential() {
    	// Length of the supercell (in units of sigma)
    		double boxSize = s.e.numberCells * s.e.latcon / s.e.sigma;
    		double maxDistance = 0;
    		double minDistance = 100;
    		for (Particle p : s.particles) {
    			if (p.numberNeighbors > 0) {
    				for (Integer jneig : p.vicinity) {
    					// Find the relative position between the two atoms
    					double relDisX = p.r.x - s.particles.get(jneig).r.x;
    					double relDisY = p.r.y - s.particles.get(jneig).r.y;
    					double relDisZ = p.r.z - s.particles.get(jneig).r.z;
    					// Consider the periodic boundary conditions. The calculation of minimum
    					// image distances is simplified by the use of reduced units:
    					// the length of the sigma parameter is taken to define the
    					// fundamental unit of length in the simulation.
    					relDisX = relDisX - boxSize * (Math.round(relDisX / boxSize));
    					relDisY = relDisY - boxSize * (Math.round(relDisY / boxSize));
    					relDisZ = relDisZ - boxSize * (Math.round(relDisZ / boxSize));
    					// Distance module
    					double modRelDistance = Math
    							.sqrt(Math.pow(relDisX, 2) + Math.pow(relDisY, 2) + Math.pow(relDisZ, 2));
    					if (modRelDistance > maxDistance) {
    						maxDistance = modRelDistance;
    					}
    					if (modRelDistance < minDistance) {
    						minDistance = modRelDistance;
    					}
    				}
    			}
    		}
    		double cutOffPotential = 4
    				* (Math.pow(sigmaReduced / s.e.cutoffRadius, 12) - Math.pow(sigmaReduced / s.e.cutoffRadius, 6));
    		List<Double> distance = new ArrayList<Double>();
    		List<Double> pot = new ArrayList<Double>();
    		List<R2Point> distancePot = new ArrayList<R2Point>();
    		double dcutOffPotential = 24 / s.e.cutoffRadius
    				* (2 * Math.pow(1 / s.e.cutoffRadius, 12) - Math.pow(1 / s.e.cutoffRadius, 6));
    		for (int i = 0; i < s.e.histogramResolution; i++) {
    			distance.add(minDistance + i * (maxDistance - minDistance) / 100);
    			pot.add(evalPotential(distance.get(i)));
    			distancePot.add(new R2Point(distance.get(i), pot.get(i)));
    		}
    		Graph potentialTest = new Graph(distancePot);
    
    		return potentialTest;
  }

  /**
   * This method takes from the class Force maxd and mind from the current state.
   * It will divide the interval between mind and maxd in some divisions which 
   * are collected in X. It will calculate the number of particles in each division
   * and collect these values on the Y axis.
   */
  public pds.common.Graph drawPairDistributionFunction() {
    			// First, initialize the maximum distance between two atoms that will be
    		// considered in the histogram (in units of sigma)
    		// Here, we take the box length
    		double maxDistance = s.e.numberCells * s.e.latcon / s.e.sigma;
    		double boxSize = maxDistance;
    		double delta = maxDistance / s.e.histogramResolution;
    		int[] histogramPair = new int[s.e.histogramResolution];
    		double[] pairDistFunct = new double[s.e.histogramResolution];
    		for (int iatom = 0; iatom < s.particles.size() - 1; iatom++) {// iatom Counter for loops on atoms
    			for (int jneig = iatom + 1; jneig < s.particles.size(); jneig++) { // Counter for loops on possible
    																				// neighbors
    				// Relative distance between two atoms
    				double revDisX = s.particles.get(iatom).r.x - s.particles.get(jneig).r.x;
    				double revDisY = s.particles.get(iatom).r.y - s.particles.get(jneig).r.y;
    				double revDisZ = s.particles.get(iatom).r.z - s.particles.get(jneig).r.z;
    				// Consider the periodic boundary conditions. The calculation of minimum
    				// image distances is simplified by the use of reduced units:
    				// the length of the sigma parameter is taken to define the
    				// fundamental unit of length in the simulation.
    				revDisX = revDisX - boxSize * (Math.round(revDisX / boxSize));
    				revDisY = revDisY - boxSize * (Math.round(revDisY / boxSize));
    				revDisZ = revDisZ - boxSize * (Math.round(revDisZ / boxSize));
    				// Distance module
    				double modRelDistance = Math.sqrt(Math.pow(revDisX, 2) + Math.pow(revDisY, 2) + Math.pow(revDisZ, 2));
    
    				// Compute the bin in the histogram where the previous distance lies in
    				int bin = (int) (modRelDistance / delta);// +1?
    				// Update the histogram.
    				// We add "2" below because we have to consider the pair iatom -> jatom,
    				// and jatom -> iatom
    				if (bin <= s.e.histogramResolution) {
    					histogramPair[bin] = histogramPair[bin] + 2;
    				}
    			}
    		}
    		// Normalize the histogram by the average number of atoms in the
    		// same interval in an ideal gas at the same density, nideal
    		// By definition, the radial distribution function is
    		// pdf(r + 1/2 delr) = histogramPair(r) / nideal(r)
    
    		double constant = 4 / 3.0 * Math.PI * s.e.reduceDensity;
    		double rLower;
    		double rUpper;
    		double nIdeal;
    		for (int i = 0; i < s.e.histogramResolution; i++) {
    			rLower = i * delta;
    			rUpper = rLower + delta;
    			nIdeal = constant * ((Math.pow(rUpper, 3) - Math.pow(rLower, 3)));
    			pairDistFunct[i] = histogramPair[i] / (nIdeal * s.particles.size());
    
    		}
    		// The appropriate distance for a particular element of our
    		// pair distribution function histogram is at the center of the
    		// interval (r + delr), i.e. at rlower + delr/2.0
    		List<R2Point> pairFunction = new ArrayList<R2Point>(s.e.histogramResolution);
    		for (int i = 0; i < s.e.histogramResolution; i++) {
    			pairFunction.add(new R2Point(i * delta + delta / 2, histogramPair[i]));
    		}
    		Graph pairDistributionFunction = new Graph(pairFunction);
    
    		return pairDistributionFunction;
  }

  /**
   * Parameter sigma of the Lennard-Jones potential in reduced units
   */
  public static final double sigmaReduced = 1.0;

  public pds.common.Graph drawForce() {
    // Length of the supercell (in units of sigma)
    		double boxSize = s.e.numberCells * s.e.latcon / s.e.sigma;
    		double maxDistance = 0;
    		double minDistance = 100;
    		for (Particle p : s.particles) {
    			if (p.numberNeighbors > 0) {
    				for (Integer jneig : p.vicinity) {
    					// Find the relative position between the two atoms
    					double relDisX = p.r.x - s.particles.get(jneig).r.x;
    					double relDisY = p.r.y - s.particles.get(jneig).r.y;
    					double relDisZ = p.r.z - s.particles.get(jneig).r.z;
    					// Consider the periodic boundary conditions. The calculation of minimum
    					// image distances is simplified by the use of reduced units:
    					// the length of the sigma parameter is taken to define the
    					// fundamental unit of length in the simulation.
    					relDisX = relDisX - boxSize * (Math.round(relDisX / boxSize));
    					relDisY = relDisY - boxSize * (Math.round(relDisY / boxSize));
    					relDisZ = relDisZ - boxSize * (Math.round(relDisZ / boxSize));
    					// Distance module
    					double modRelDistance = Math
    							.sqrt(Math.pow(relDisX, 2) + Math.pow(relDisY, 2) + Math.pow(relDisZ, 2));
    					if (modRelDistance > maxDistance) {
    						maxDistance = modRelDistance;
    					}
    					if (modRelDistance < minDistance) {
    						minDistance = modRelDistance;
    					}
    				}
    			}
    		}
    		List<Double> distance = new ArrayList<Double>();
    		List<Double> force = new ArrayList<Double>();
    		List<R2Point> distanceForce = new ArrayList<R2Point>();
    		double dcutOffPotential = 24 / s.e.cutoffRadius
    				* (2 * Math.pow(1 / s.e.cutoffRadius, 12) - Math.pow(1 / s.e.cutoffRadius, 6));
    		for (int i = 0; i < s.e.histogramResolution; i++) {
    			distance.add(minDistance + i * (maxDistance - minDistance) / 100);
    			if (distance.get(i) < s.e.cutoffRadius) {
    				force.add((24 / Math.pow(distance.get(i), 1)
    						* (2 * Math.pow(1 / distance.get(i), 12) - Math.pow(1 / distance.get(i), 6))
    						- dcutOffPotential));
    			} else {
    				force.add(0.0);
    			}
    			distanceForce.add(new R2Point(distance.get(i), force.get(i)));
    		}
    		Graph forceTest = new Graph(distanceForce);
    
    		return forceTest;
  }

}
