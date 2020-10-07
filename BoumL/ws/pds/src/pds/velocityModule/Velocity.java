
package pds.velocityModule;

import pds.common.*;
public class Velocity implements pds.common.IVelocity {
  /**
   * The state for which velocities and positions are calculated.
   * This value is initialized by the constructor.
   */
  protected pds.common.State s;

  /**
   * The object that will be used to do the calculus of accelerations,
   * It is initialized by the constructor.
   */
  protected pds.common.IForce f;

  /**
   * This static method returns a Velocity object that implements the interface IVelocity
   * 
   * @param s the state for whose particles the velocities and positions will be 
   * calculated when requested
   * @param f the force calculation object that will assign accelerations
   * @return an object of a class that implements the interface IVelocity
   */
  public static pds.common.IVelocity createVelocity(pds.common.State s, pds.common.IForce f)
  {
    
    		  return new Velocity(s, f);
  }

  /**
   * This constructor receives the state for which the Velocity object
   * will be created and the force object that calculates the acceleration
   * of all particles.
   * 
   * @param s the state for whose particles the velocities and positions 
   * will be calculated when requested
   * @param f the force calculation object that will assign accelerations
   */
  public Velocity(pds.common.State s, pds.common.IForce f) {
    	this.s = s;
    	this.f = f;
  }

  /**
   * This operation calculates the positions and velocities for the given 
   * state and force objects. It returns the largest of the displacements of the 
   * particles of the system since the last calculation of the neighbours lists.
   * It must update the kinetic energy of the current state.
   * 
   * @return the largest of the displacements of the particles since the last time the vicinity was calculated. 
   */
  public double calcVelocity() {
    double kineticEnergy = 0;
    
    		// The first part of the algorithm is a Taylor series which advances positions
    		// form t to dt+t
    		// and velocities from t to dt/2+t. After this the force routine is called
    		for (int iatom = 0; iatom < s.particles.size(); iatom++) {
    			// Update the positions
    			// First we allocate the new positions in a local variable so that we can work
    			// with the previous
    			// positions to calculate the displacement for each particle
    
    			// X coordinate
    			double newX = 0;
    			newX = s.particles.get(iatom).r.x + s.e.timeStep * s.particles.get(iatom).v.x
    					+ 0.5 * Math.pow(s.e.timeStep, 2) * s.particles.get(iatom).a.x;
    
    			// Y coordinate
    			double newY = 0;
    			newY = s.particles.get(iatom).r.y + s.e.timeStep * s.particles.get(iatom).v.y
    					+ 0.5 * Math.pow(s.e.timeStep, 2) * s.particles.get(iatom).a.y;
    
    			// Z coordinate
    			double newZ = 0;
    			newZ = s.particles.get(iatom).r.z + s.e.timeStep * s.particles.get(iatom).v.z
    					+ 0.5 * Math.pow(s.e.timeStep, 2) * s.particles.get(iatom).a.z;
    
    			// Creation of an object R3Point to calculate the distance between the previous
    			// and the current
    			// position of the atom iatom
    			R3Point newPositions = new R3Point(newX, newY, newZ);
    			// Calculation of the displacement d
    			s.particles.get(iatom).d = s.particles.get(iatom).d + s.particles.get(iatom).r.distanceTo(newPositions);
    
    			// Now, we can allocate the new positions in the positions of the atom iatom
    			// because the previous positions will not use any more.
    			s.particles.get(iatom).r.x = newX;
    			s.particles.get(iatom).r.y = newY;
    			s.particles.get(iatom).r.z = newZ;
    
    			// Compute the velocities at half step
    			// Vx component
    			s.particles.get(iatom).v.x = s.particles.get(iatom).v.x + s.e.timeStep / 2.0 * s.particles.get(iatom).a.x;
    			// Vy component
    			s.particles.get(iatom).v.y = s.particles.get(iatom).v.y + s.e.timeStep / 2.0 * s.particles.get(iatom).a.y;
    			// Vz component
    			s.particles.get(iatom).v.z = s.particles.get(iatom).v.z + s.e.timeStep / 2.0 * s.particles.get(iatom).a.z;
    		}
    
    		// Compute the forces for the atomic positions at (t + \delta t)
    		s.energy.u = f.calcForce();
    
    		// The second part of the algorithm advances velocities from t+dt/2 to t+dt
    		for (int iatom = 0; iatom < s.particles.size(); iatom++) {
    			// Compute the velocities at half step
    
    			// Vx component
    			s.particles.get(iatom).v.x = s.particles.get(iatom).v.x + s.e.timeStep / 2.0 * s.particles.get(iatom).a.x;
    
    			// Vy component
    			s.particles.get(iatom).v.y = s.particles.get(iatom).v.y + s.e.timeStep / 2.0 * s.particles.get(iatom).a.y;
    
    			// Vz component
    			s.particles.get(iatom).v.z = s.particles.get(iatom).v.z + s.e.timeStep / 2.0 * s.particles.get(iatom).a.z;
    			kineticEnergy = kineticEnergy + s.particles.get(iatom).v.kinetic();
    		}
    		s.energy.k = kineticEnergy;
    		// TotalEnergy
    		s.energy.total = s.energy.k + s.energy.u;
    		s.nStep++;
    
    		// Now it will be calculated the largest of the displacements of the particles
    		// since the last time the vicinity was calculated to check later if it is
    		// necessary to recalculate
    		// the list of neighbors.
    
    		/////////////////////// DUDA////////////////////
    		// Maximum displacement of all the atoms
    		double maxDisplacement = 0;
    		// Find the maximum displacement between all the atoms in the supercell
    		for (int iatom = 0; iatom < s.particles.size(); iatom++) {
    			if (maxDisplacement < s.particles.get(iatom).d) {
    				maxDisplacement = s.particles.get(iatom).d;
    			}
    		}
    		return maxDisplacement;
  }

}
