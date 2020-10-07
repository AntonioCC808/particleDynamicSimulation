
package pds.common;

public class EnergyLog implements Cloneable {
  /**
   * Returns a copy of the current object
   * 
   * @return a copy of the current object
   */
  public Object clone() {
            Object obj=null;
            try{
                obj=super.clone();
            }catch(Exception ex){ //Use CloneNotSupportedException if possible
                System.out.println(" It was not possible to clone...");
            }
            return obj;
  }

  /**
   * *
   *  * The potential  energy of the system at the current state, it is hold in  a double primitive value
   *  
   */
  public double u;

  /**
   * *
   *  * The kinetic energy of the system at the current state, it is hold in  a double primitive value
   *  
   */
  public double k;

  /**
   * The total  energy of the system at the current state (sum of the kinetic and potential energy), it is hold in  a double primitive value
   */
  public double total;

}
