
package pds.neighboursModule;

import pds.common.*;
import java.util.*;
/**
 * This class calculates the current neighbours lists for all particles.
 */
public class Neighbours implements pds.common.INeighbours {
  /**
   * The state for which the vecinities of all particles are calculated.
   * This value is initialized by the constructor.
   */
  protected pds.common.State s;

  /**
   * This static method returns a Neighbours object that implements the interface INeighbours
   * 
   * @param s the state for whose particles the vicinity will be calculated when needed
   * @return an object of a class that implements the interface INeighbours
   */
  public static pds.common.INeighbours createNeighbours(pds.common.State s)
  {
    	 
    		  return new Neighbours(s);
    
  }

  /**
   * This constructor receives the state for which the Neighbours object
   * will be created.
   * 
   * @param s the state for whose particles the vecinity will be calculated when needed
   */
  public Neighbours(pds.common.State s) {
          this.s = s;
  }

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
  public boolean calcVecinity(double maxDisplacement) {
 //Some variables that will be needed
	  //Length of the supercell (in units of sigma)
	  double boxSize=s.e.numberCells*s.e.latcon/s.e.sigma;
	  if (maxDisplacement>(s.e.cutoffList-s.e.cutoffRadius)) {
		  
		  for (int iatom=0;iatom<s.particles.size()-1;iatom++) {//iatom Counter for loops on atoms
			  //Find the relative position between the two atoms
			  s.particles.get(iatom).vicinity.clear();
			  
			  for (int jneig=iatom+1;jneig<s.particles.size();jneig++) { //Counter for loops on possible neighbors
				  //Relative distance between two atoms
				  double revDisX=s.particles.get(iatom).r.x-s.particles.get(jneig).r.x;
				  double revDisY=s.particles.get(iatom).r.y-s.particles.get(jneig).r.y;
				  double revDisZ=s.particles.get(iatom).r.z-s.particles.get(jneig).r.z;
				  // Consider the periodic boundary conditions. The calculation of minimum 
				  // image distances is simplified by the use of reduced units: 
				  // the length of the sigma parameter is taken to define the
				  //fundamental unit of length in the simulation. 
				  revDisX=revDisX-boxSize*(Math.round(revDisX/boxSize) );
				  revDisY=revDisY-boxSize*(Math.round(revDisY/boxSize));
				  revDisZ=revDisZ-boxSize*(Math.round(revDisZ/boxSize));
				  //Distance module
				  double modRelDistance=Math.sqrt(Math.pow(revDisX, 2)+Math.pow(revDisY, 2)+Math.pow(revDisZ, 2));
				  //Compare the square of the distance with the square of the outer radius of the Verlet list.
				  if (modRelDistance<s.e.cutoffList) {
					  s.particles.get(iatom).vicinity.add(jneig);
				  }	  
				  s.particles.get(iatom).numberNeighbors=s.particles.get(iatom).vicinity.size();
			  }
			  s.particles.get(iatom).d=0;
		  }
		  s.particles.get(s.particles.size()-1).d=0;
		  return true;
	  } 
	  return false;
  }

}
