
package pds.common;

public class R3Point extends R2Point implements Cloneable {
  /**
   * Coordinate in the z axis for the three dimensional value
   *  
   */
  public double z;

  /**
   * Calculates the distance between two points in R3 space (module of the difference)
   * 
   * @param the point to which the distance will be calculated
   * @return the module of the distance between this and the point given as argument, 
   * it is returned as a double value.
   *  
   */
  public double distanceTo(R3Point r3Point) {
    	  // Calculation of the relative position between the two atoms
    		double relavPosX = r3Point.x - x;
    		double relavPosY = r3Point.y - y;
    		double relavPosZ = r3Point.z - z;
    		// Module of the the relative distance
    		double distance = Math.sqrt(Math.pow(relavPosX, 2) + Math.pow(relavPosY, 2) + Math.pow(relavPosZ, 2));
    		// The distance between two points in space will be the subtracting of their
    		// modules
    		return distance;
  }

  public R3Point(double x, double y, double z) {
    	  super(x,y);
    	  this.z = z;
  }

  public R3Point() {
    	  super();
  }

  /**
   * Returns a copy of the current object
   * 
   * @return a copy of the current object
   */
  public Object clone() {
            Object obj=null;
            try{
                obj=super.clone();
            }catch(Exception ex){  //Use CloneNotSupportedException if possible
                System.out.println(" It was not possible to clone...");
            }
            return obj;
  }

  /**
   * Calculates the Kinetic Energy of the particle with Vx, Vy and Vz
   */
  public double kinetic() {
    double kineticEnergy = 1 / 2.0 * (Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2));
    		return kineticEnergy;
  }

}
