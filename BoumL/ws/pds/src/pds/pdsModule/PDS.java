
package pds.pdsModule;

import pds.common.*;
import pms.initModule.*;
import pds.outputModule.*;
import pds.forceModule.*;
import pds.velocityModule.*;
import pds.neighboursModule.*;
import java.util.*;
public class PDS {
  /**
   * This is the current step of the simulation
   */
  private int step;

  /**
   * This is the handler for the concurrent thread in which the
   * simulation will run if requested.
   */
  private Thread simulationThread;

  /**
   * Default constructor
   */
  public PDS() {
  }

  /**
   * This operation runs a simulation of the given experiment up to completion and produces 
   * the output results. If it is running in a concurrent thread it can be
   * paused and continued.
   * 
   * @param experiment is the experiment object with the data needed to configure and 
   * run the simulation
   * @return returns the ResultsFile object with the names of all the requested 
   * resulting images and data files generated. If any of those files could not 
   * be generated the corresponding string would be a null reference.
   */
  private pds.common.ResultsFiles simulation(pds.common.InputData experiment) {
            
            	 
    		// prepare the list to take the log of the system (energies)
    		experiment.energyLog = new ArrayList<EnergyLog>(experiment.nSteps);
    		// prepare the list to take the movie samples
    		experiment.movie = new ArrayList<State>(experiment.nSteps);
    
    		/* Initial setup for the simulation */
    
    		// Create the initialization object
    		IInit init = pms.initModule.Init.createInit();
    		// Invoke doInit() to get the initial state
    		State s = init.doInit(experiment);
    		experiment.InitialVelocityDistribution = init.getInitialVelocityDistribution();
    
    		// Create the neighbors list management object
    		INeighbours neighbours = pds.neighboursModule.Neighbours.createNeighbours(s);
    
    		// Set the initial neighbors lists
    		neighbours.calcVecinity(5);
    
    		State dataFile = s.copy(s.particles);
    		// Create the force calculation object
    		IForce force = pds.forceModule.Force.createForce(s);
    
    		// Create the Velocity calculation object
    		IVelocity velocity = pds.velocityModule.Velocity.createVelocity(s, force);
    
    		// calculate the initial accelerations
    		// and take the potential energy it in the system log object
    		s.energy.u = force.calcForce();
    
    		// Set p the initial total energy of the system.
    		s.energy.total = s.energy.u + s.energy.k;
    
    		// Graphs to plot later the potential energy and the force for the initial state
    		experiment.potentialGraph = force.drawPotential();
    		experiment.forceGraph = force.drawForce();
    
    		// store the first state log in the experiment log list
    		experiment.energyLog.add((EnergyLog) s.energy.clone());
    
    		// create the list of particles in the visualization scope
    		List<Particle> sampleScope = new ArrayList<Particle>();
    
    
    		// take the particles of interest from the complete list into the sampleScope
    		for (Particle p : s.particles) {
    			sampleScope.add(p);
    		}
    
    		/* Do the simulation */
    		simulationLoop: for (step = 0; step < experiment.nSteps - 1; step++) {
    
    			// We simulate a time step for the current state.
    			// and update the vicinity if needed
    			neighbours.calcVecinity(velocity.calcVelocity());
    
    			// store the current state log in the experiment log list
    			experiment.energyLog.add((EnergyLog) s.energy.clone());
    
    			// get a copy of the state with the sample of particles
    			// that fit in the visualization scope
    			experiment.movie.add(s.copy(sampleScope));
    
    			
    		}
    		experiment.PairDistributionFunctionGraph = force.drawPairDistributionFunction();
    		// Create the output production object
    		IOutput out = new pds.outputModule.Output(experiment,
    				new ResultsFiles(experiment.outputFolder, experiment.nickname), dataFile);
    		return out.results();
  }

  public static void main(String [] args) throws pds.common.InputData.ExperimentNotCreated
  {
    		String nomFich = "Input\\InputFile";
    		if (args.length != 0 && args[0] != null)
    			nomFich = args[0];
    		InputData e = new InputData(nomFich);
    
    		PDS run=new PDS();
    
    		ResultsFiles r = run.simulation(e);
  }

}