
package pds.common;

import java.util.*;
import java.io.*;
/**
 * Container of all the information needed to perform the simulation. It reads the data 
 * from the text file passed as parameter.
 */
public class InputData {
  /**
   * Number of steps of the simulation
   */
  public int nSteps;

  /**
   * Time step used for the simulation.
   */
  public double timeStep;

  /**
   * Simulation identifier.
   */
  public String nickname;

  /**
   * Parameter with length dimension of the biparametric potential. Sigma for the Lennard-Jones potential, which is the distance for which the potentialbetween particles is zero .
   */
  public double sigma;

  /**
   * Energy dimension parameter of the biparametric potential. Epsilon for the
   * Lennard-Jones potential, which is the depth of the potential.
   */
  public double epsilon;

  /**
   * Reduced density of the bulk of atoms.
   */
  public double reduceDensity;

  /**
   * Mass of the simulated atoms.
   */
  public double mass;

  /**
   * Atomic number of the simulated atom.
   */
  public int atomicNumber;

  /**
   * Initial temperature.
   */
  public double temperature;

  /**
   * Number of cells in each direction (x,y,z) of the supercell.
   */
  public int numberCells;

  /**
   * Cutoff radius considered for the Lennard-Jones potential (in sigma units).
   */
  public double cutoffRadius;

  /**
   * Cutoff distance to compute the Verlet neighbour list
   */
  public double cutoffList;

  /**
   * 
   * Type of potential used in the simulation.
   * 
   */
  public IForce.TypeOfPotential potential;

  /**
   * 
   * Type of the initial lattice.
   * 
   */
  public IInit.TypeOfCellLattice initialLattice;

  /**
   * Full path of the folder where all files will be generated.
   */
  public final String outputFolder;

  /**
   * Number of segments in the X axis used for the various histograms 
   * plotted as results of the simulation. This number is also used as the 
   * number of values of distance used to generate the image of the potential. 
   */
  public int histogramResolution;

  /**
   * List of data collected from each step of the simulation. It has the
   * values of the kinetic and potential energies at each step of the simulation.
   */
  public List<EnergyLog> energyLog;

  /**
   * List of states that will be used for creating the animation. Only 
   * a fraction of the steps is stored in this list, and only a fraction of the lattice is 
   * stored in the set of particles
   */
  public List<State> movie;

  /**
   * 
   * This reference serve to provide the InitialVelocityDistribution graph to the output module.
   * 
   */
  public Graph InitialVelocityDistribution;

  /**
   * Constructor of the object. Reads the file whose full path is given as input 
   * argument, initializes the attributes and writes the log file.
   * 
   * @param simulationData a string with the full path of the file that gives values to all the
   * attributes of the object.
   * @throws ExperimentNotCreated 
   */
  public InputData(String simulationData) throws InputData.ExperimentNotCreated {
	// Variables which will be compared in the switch with the data read from the
		// input file
		final String ReducedTimeStep = "ReducedTimeStep";
		final String NumberSteps = "NumberSteps";
		final String NickName = "NickName";
		final String Sigma = "Sigma";
		final String Epsilon = "Epsilon";
		final String ReducedDensity = "ReducedDensity";
		final String AtomicNumber = "AtomicNumber";
		final String Temperature = "Temperature";
		final String NumberCells = "NumberCells";
		final String CutOffRadius = "CutOffRadius";
		final String CutOffList = "CutOffList";
		final String TypeOfPotential = "TypeOfPotential";
		final String TypeOfLattice = "TypeOfLattice";
		final String Mass = "Mass";
		final String HistogramResolution = "HistogramResolution";
		String s;
		// creates the experiment with the data in the text file
		// Use "." as decimal separator
		Locale.setDefault(Locale.ENGLISH);
		Scanner in = null;
		try {
			// Open the file
			in = new Scanner(new FileReader(simulationData));
			in.hasNext();
			while (in.hasNext()) {
				s = in.next();
				switch (s) {

				case ReducedTimeStep:
					timeStep = in.nextDouble();
					break;
				case NumberSteps:
					nSteps = in.nextInt();
					break;
				case NickName:
					nickname = in.next();
					break;
				case Sigma:
					sigma = in.nextDouble();
					break;
				case Epsilon:
					epsilon = in.nextDouble();
					break;
				case ReducedDensity:
					reduceDensity = in.nextDouble();
					break;
				case AtomicNumber:
					atomicNumber = in.nextInt();
					break;
				case Mass:
					mass = in.nextDouble();
					break;
				case Temperature:
					temperature = in.nextDouble();
					break;
				case NumberCells:
					numberCells = in.nextInt();
					break;
				case CutOffRadius:
					cutoffRadius = in.nextDouble();
					break;
				case CutOffList:
					cutoffList = in.nextDouble();
					break;
				case HistogramResolution:
					histogramResolution = in.nextInt();
					break;
				case TypeOfPotential:
					String typep=in.next();
					if (typep.equals("LENNARDJONES-TOTALSHIFTED")) {
						potential = IForce.TypeOfPotential.LENNARDJONESTOTALSHIFTED;
					}else if(typep.equals("LENNARDJONES-SHIFTED")){
						potential = IForce.TypeOfPotential.LENNARDJONESSHIFTED;
					}else if(typep.equals("LENNARDJONES-UNSHIFTED")) {
						potential = IForce.TypeOfPotential.LENNARDJONESUNSHIFTED;
					}
					break;
				case TypeOfLattice:
					String typel=in.next();
					if (typel.equals("FCC")) {
						initialLattice = IInit.TypeOfCellLattice.FCC;
					}else if(typel.equals("BCC")){
						initialLattice = IInit.TypeOfCellLattice.BCC;
					}
					break;
				default:
				}

			}
		} catch (FileNotFoundException e) {
			System.out.println("Error opening file: " + simulationData);
		} finally {
			if (in != null) {
				// close the file
				in.close();
			}
		}
		outputFolder = "Output\\";
		histogramResolution = 100;
		latcon = sigma * Math.pow(4 / reduceDensity, 1 / 3.0);
  }

  @SuppressWarnings("serial")
  public class ExperimentNotCreated extends Exception {
  }

  /**
   * 
   * This reference serve to provide the potential graph to the output module.
   * 
   */
  public Graph potentialGraph;

  /**
   * 
   * This reference serve to provide the PairDistributionFunction graph to the output module.
   * 
   */
  public Graph PairDistributionFunctionGraph;

  public Graph forceGraph;

  /**
   * Lattice Constant
   */
  public double latcon;

}
