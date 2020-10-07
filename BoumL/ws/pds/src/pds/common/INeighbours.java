
package pds.common;

public interface INeighbours {
  /**
   * This methode verifies the need to recalculate the vecinity of the particles
   * and do it if necessary. 
   * It accumulates the maxDisplacements input argument and re-calculates the 
   * vecinities if the accumulated displacement exceeds the maximum tolerable.
   * If the vecinities need to be recalculated it calculates all the vecinities lists 
   * and reset the accumulated displacement to 0.0.
   * 
   * @param maxDisplacement the value of the largest displacement suffered in the 
   * ongoing step by any of the particles in the system.
   * @return true if the vecinities have been recalculated false otherwise
   */
  boolean calcVecinity(double maxDisplacement) ;

}
