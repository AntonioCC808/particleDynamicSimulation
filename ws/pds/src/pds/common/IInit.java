
package pds.common;

public interface IInit {
  /**
   * Kinds of images that can be requested to the results formating module
   *  
   */
  public enum TypeOfCellLattice {
    FCC,/**
     * This indicates that the particles in the cell lattice follow a Face Centered Cubic pattern.
     *  
     */

    BCC,/**
     *  This indicates that the particles in the cell lattice follow a Body Centered Cubic pattern
     */
;
  }

  
  /**
   * Create the initial state for the given experiment and returns it with the 
   * initial positions, velocities and neighbours lists for all particles.
   * 
   * @param experiment is the experiment object with the data needed to configure and 
   * run the simulation
   * @return a fresh new State ready to start the simulation
   */
  State doInit(InputData experiment) ;

  /**
   * Calculates a Graph showing the statistical distribution of velocities for 
   * all the particles of the system at the given state.
   * 
   * @param s a reference to the State object for which the distribution will be calculated
   * @return a reference to the Graph object that holds the velocities distribution function 
   */
  Graph velocityDistribution(State s) ;

  /**
   * Returns a Graph showing the statistical distribution of velocities for 
   * all the particles of the system as they were assigned at the initial state.
   * 
   * @return a reference to the Graph object that holds the velocities distribution function 
   */
  Graph getInitialVelocityDistribution() ;

}
