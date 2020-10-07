
package pds.common;

public interface INeighbours {
  /**
   * This method verifies the need to recalculate the vicinity of the particles
   * and do it if necessary. 
   * It accumulates the maxDisplacements input argument and re-calculates the 
   * Vicinities if the accumulated displacement exceeds the maximum tolerable.
   * If the vicinities need to be recalculated it calculates all the vicinities lists 
   * and reset the accumulated displacement to 0.0.
   * 
   * @param maxDisplacement the value of the largest displacement suffered in the 
   * ongoing step by any of the particles in the system.
   * @return true if the vicinities have been recalculated false otherwise
   */
  boolean calcVecinity(double maxDisplacement) ;

}
