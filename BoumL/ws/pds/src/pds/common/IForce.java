
package pds.common;

/**
 * This method takes from the class Force maxd and mind from the current state.
 * It will divide the interval between mind and maxd in some divisions which are collected in X. It will calculate the number of particles in each division
 * and collect these values on the Y axis.
 */
public interface IForce {
  /**
   * Kinds of images that can be requested to the results formating module
   *  
   */
  public enum TypeOfPotential {
    LENNARDJONESUNSHIFTED,/**
     * This indicates that the forces that affect particles are to be calculated using the unshifted Lennard-Jones kind of potential.
     *  
     */

    LENNARDJONESSHIFTED,/**
     * This indicates that the forces that affect particles are to be calculated using the shifted Lennard-Jones potential to avoid discontinuity in the potential energy.
     */

    LENNARDJONESTOTALSHIFTED,/**
     * This indicates that the forces that affect particles are to be calculated using the total shifted Lennard-Jones potential to avoid discontinuity in the potential energy
     * and the force
     */
;
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
  double calcForce() ;

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
  double evalPotential(double distance) ;

  /**
   * This method returns a graph in the form of a list of points (x,y) that 
   * describe the potential for the range of distance vales observed in the 
   * current experiment.
   * The number of points is taken from the histogram size constant that comes from Input.
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
  Graph drawPotential() ;

  /**
   * This method takes from the class Force maxd and mind from the current state.
   * It will divide the interval between mind and maxd in some divisions which 
   * are collected in X. It will calculate the number of particles in each division
   * and collect these values on the Y axis.
   */
  Graph drawPairDistributionFunction() ;

  Graph drawForce() ;

}
