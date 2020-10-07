
package pds.common;

import java.util.*;
/**
 * This class collects the data needed to go through the steps of the simulation.
 * It holds the complete lattice in the form of a list of particles.
 */
public class State {
  public List<Particle> particles;

  /**
   * *
   *  * The configuration and simulation data that defines the experiment being executed
   *  
   */
  public InputData e;

  /**
   * *
   *  * number of the step in the simulation corresponding to the state object
   *  * step 0 is the initialization state
   *  
   */
  public int nStep;

  /**
   * Default constructor
   */
  public State() {
    super();
  }

  /**
   * This constructor is to be used only at the initialization of the system.
   * The initialization step will be set to zero.
   * @param e is the experiment data object used to create this initial state
   * @param particles is the list of all particles of the system
   * @energy an object of the class OutLog with the initial potential and kinetic energies
   */
  public State(InputData e, List<Particle> particles, EnergyLog energy) {
    this.nStep     = 0;
    this.e         = e;
    this.particles = particles;
    this.energy    = energy;
  }

  /**
   * This operation returns a copy of the state with the relevant information for 
   * the creation of the animations.
   */
  public State copy(List<Particle> scope) {
    	  State s = new State();
    	  s.nStep = this.nStep;
    	  s.e = this.e;
    	  s.particles = new ArrayList<Particle>(scope.size());
    	  for (Particle p: scope) s.particles.add((Particle) p.clone());
    	  
    	  return s;
  }

  /**
   * The object that holds the values of the energy at the current state
   */
  public EnergyLog energy;

}
