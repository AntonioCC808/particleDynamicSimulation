
package pds.common;

import java.util.*;
public class Particle implements Cloneable {
	 
	/**
	   * The position of the particle
	   *  
	   */
	  public R3Point r;

	  /**
	   * The velocity of the particle
	   *  
	   */
	  public R3Point v;

	  /**
	   * The acceleration of the particle
	   *  
	   */
	  public R3Point a;

	  /**
	   * The maximum possible displacement of the particle since the last
	   * time the vicinity has been calculated. This attribute is updated by
	   * the velocity calculation module at each step. It is reset to 0.0 
	   * by the neighbors module when the vicinity is recalculated.
	   */
	  public double d;
	
	  
	  /**
	   * Indexes of the particles that are in the vicinity of this particle 
	   */
	  public Set<Integer> vicinity;
	  
	  /**
	   * Number of neighbors for this particle 
	   */
	  public int numberNeighbors;
	  
	  
 

  /**
   * This constructor creates the particle assigning only the 
   * position of the particle using the coordinates given 
   * as arguments x, y and z. The velocity and the acceleration
   * are also created but initialized with 0.0
   */
  public Particle(double x, double y, double z) {
	  this.r = new R3Point(x,y,z);
	  this.v = new R3Point(0.0,0.0,0.0);
	  this.a = new R3Point(0.0,0.0,0.0);
	  this.vicinity=new HashSet<Integer> ();
	  this.numberNeighbors=0;


}

  /**
   * Returns a copy of the current object
   * 
   * @return a copy of the current object
   */
  public Object clone() {
      Particle obj=null;
      try{
          obj=(Particle) super.clone();
          obj.a = (R3Point) obj.a.clone();
          obj.v = (R3Point) obj.v.clone();
          obj.r = (R3Point) obj.r.clone();
      }catch(CloneNotSupportedException ex){
          System.out.println(" It was not possible to clone...");
      }
      return obj;
}

 

}
