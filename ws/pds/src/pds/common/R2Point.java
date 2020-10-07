
package pds.common;

public class R2Point implements Cloneable {
	/**
	 * The x coordinate in an cartesian space
	 */
	public double x;

	/**
	 * The y coordinate in a cartesian space
	 */
	public double y;

	/**
	 * Constructor of the class that assign values to x and y
	 */
	public R2Point(double x, double y) {
		this.x = x;
		this.y = y;
	}

	/**
	 * Returns a copy of the current object
	 * 
	 * @return a copy of the current object
	 */
	public Object clone() {
		Object obj = null;
		try {
			obj = super.clone();
		} catch (Exception ex) { // Use CloneNotSupportedException if possible
			System.out.println(" It was not possible to clone...");
		}
		return obj;
	}

}
