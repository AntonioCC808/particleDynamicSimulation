
package pds.outputModule;

import java.util.*;
import pds.common.*;
import java.io.*;
import pds.common.*;


@SuppressWarnings("unused")
public class Output implements pds.common.IOutput {

	/**
	 * 
	 * The experiment object, with (almost) all information to generate output, with
	 * particles positions and with all constants necessary to convert units.
	 * 
	 */
	protected pds.common.InputData experiment;

	/**
	 * 
	 * The results files object, whose names are taken, edited (with the extension),
	 * and then returned.
	 * 
	 */
	protected pds.common.ResultsFiles resultsFiles;

	/**
	 * 
	 * The results files object, whose names are taken, edited (with the extension),
	 * and then returned.
	 * 
	 */
	protected pds.common.State s;

	/**
	 * This converts the internal units to Armstrong (for the output in XCryDen).
	 */
	private double lengthConverterA;

	/**
	 * The factor to convert the velocity units from internal to the ones of the
	 * movies.
	 */
	private double velocityConverterA;

	/**
	 * The factor to convert the force (or accelerations) units from internal to the
	 * ones of the movies.
	 */
	private double forceConverterA;

	/**
	 * This is the number of steps the movie will have.
	 */
	private int nSteps;

	/**
	 * This static method returns an OutputModule object that implements the
	 * interface IOutput. It does nothing else, as experiment data is not given yet.
	 * 
	 * @return an object that implements the interface IOutput
	 */
	public static pds.common.IOutput createOutput(pds.common.InputData experiment, pds.common.ResultsFiles resultsFiles,
			pds.common.State s) {
		return new Output(experiment, resultsFiles, s);
	}

	/**
	 * This constructor receives the inputData and the name of the results files for
	 * which the Output object will be created.
	 * 
	 * @param experiment Is the experiment object with the data needed to retrieve
	 *                   the simulation results (without the extensions).
	 */
	public Output(pds.common.InputData experiment, pds.common.ResultsFiles resultsFiles, pds.common.State s) {
		this.experiment = experiment;
		this.resultsFiles = resultsFiles;
		this.s = s;
		this.nSteps = experiment.nSteps;
		this.lengthConverterA = experiment.sigma / 1.0E-10;// Transform sigma units to Angstrom
		this.forceConverterA = experiment.epsilon / experiment.sigma * experiment.mass ;// Transform
																														// reduced
																														// units
																														// to
																														// Hartree/Angstrom
		this.velocityConverterA = 1.38E-23 / experiment.mass;
	}

	/**
	 * This makes the results for the experiment, making the text files. The names
	 * of all the files are written in the resultFiles object (Overwriting it).
	 * 
	 * @return returns The ResultsFile object with the names of all the requested
	 *         resulting images and data files generated. If any of those files
	 *         could not be generated the corresponding string would be a null
	 *         reference. This includes the extensions of each file. Also, this is
	 *         the same as the resultFiles given object.
	 */
	public ResultsFiles results() {
		// This function creates the files from the graphs that have been stored along
		// the simulation
		positionsTextFile(s);
		velocitiesTextFile(s);
		velocityTextFile(experiment.InitialVelocityDistribution);
		energyTextFile(experiment.energyLog);
		potentialAndForceTextFile(experiment.potentialGraph, experiment.forceGraph);
		neighbourTextFile(s);
		pairDistributionTextFile(experiment.PairDistributionFunctionGraph);
		positionMovie(experiment.movie);
		return resultsFiles;
	}

	/**
	 * Method that makes a .txt file with the velocity (in its units) and with the
	 * Frequency of each one. This .txt file contains 2 columns, one with the
	 * velocities and another with the frequencies.
	 * 
	 * @param velocityData the Graph object with the velocities and frequencies
	 */
	public void positionsTextFile(State s) {

		PrintWriter out = null;
		try {
			// Open the file
			out = new PrintWriter(new FileWriter(resultsFiles.positionsData + ".txt"));
			// Write the data in the file with this format
			out.println("# Positions of the atoms in each direction in reduced units at the initial state for  "
					+ experiment.nickname);
			out.printf(" %-25s %20s %20s %20s %n", "Atom", "r(x)", "r(y)", "r(z)");

			for (int i = 0; i < s.particles.size(); i++) {
				out.printf(" %-25s %20.6f %20.6f %20.6f %n", i + 1, s.particles.get(i).r.x, s.particles.get(i).r.y,
						s.particles.get(i).r.z);
			}
		} catch (IOException e) {
			System.out.println("Error al abrir " + resultsFiles.positionsData);
			resultsFiles.positionsData = null;
		} finally {
			if (out != null)
				// Close the file
				out.close();

		}
	}

	/**
	 * Method that makes a .txt file with the velocity (in its units) and with the
	 * Frequency of each one. This .txt file contains 2 columns, one with the
	 * velocities and another with the frequencies.
	 * 
	 * @param velocityData the Graph object with the velocities and frequencies
	 */
	public void velocitiesTextFile(State s) {

		PrintWriter out = null;
		try {
			// Open the file
			out = new PrintWriter(new FileWriter(resultsFiles.velocitiesData + ".txt"));
			// Write the data in the file with this format
			out.println("# Velocities of the atoms in each direction in reduced units at the initial state for  "
					+ experiment.nickname);
			out.println(
					"# The variance of the normal (Gaussian) distribution is = " + Math.sqrt(experiment.temperature));
			out.printf(" %-25s %20s %20s %20s %n", "Atom", "V(x)", "V(y)", "V(z)");

			for (int i = 0; i < s.particles.size(); i++) {
				out.printf(" %-25s %20.5f %20.5f %20.5f %n", i + 1, s.particles.get(i).v.x, s.particles.get(i).v.y,
						s.particles.get(i).v.z);
			}
		} catch (IOException e) {
			System.out.println("Error al abrir " + resultsFiles.velocitiesData);
			resultsFiles.velocitiesData = null;
		} finally {
			if (out != null)
				// Close the file
				out.close();

		}
	}

	/**
	 * Method that makes a .txt file with the velocity (in its units) and with the
	 * Frequency of each one. This .txt file contains 2 columns, one with the
	 * velocities and another with the frequencies.
	 * 
	 * @param velocityData the Graph object with the velocities and frequencies
	 */
	public void velocityTextFile(Graph velocityData) {

		PrintWriter out = null;
		try {
			// Open the file
			out = new PrintWriter(new FileWriter(resultsFiles.velocitiesHistogramData + ".txt"));
			// Write the data in the file with this format
			out.println("# Histogram for the intial velocities of " + experiment.nickname);
			out.println("# The velocities are initially chosen randomly from a Gaussian distribution                ");
			out.println("# The mean value of the normal (Gaussian) distribution is zero");
			out.println(
					"# The variance of the normal (Gaussian) distribution is = " + Math.sqrt(experiment.temperature));
			out.printf("# %-25s %20s %n", "Velocity(reduced units)", "Frequency");

			for (int i = 0; i < velocityData.points.size(); i++) {
				out.printf(" %-25f %20.0f %n", velocityData.points.get(i).x, velocityData.points.get(i).y);
			}
		} catch (IOException e) {
			System.out.println("Error al abrir " + resultsFiles.velocitiesHistogramData);
			resultsFiles.velocitiesHistogramData = null;
		} finally {
			if (out != null)
				// Close the file
				out.close();

		}
	}

	/**
	 * Method that makes a .txt file with the neighbour list of each atom. This .txt
	 * file contains the positions of the neighbors for each atom.
	 * 
	 * @param s the initial state in which the neighbour list and the positions of
	 *          the atoms are
	 */
	public void neighbourTextFile(pds.common.State s) {

		PrintWriter out = null;
		try {
			// Open the file
			out = new PrintWriter(new FileWriter(resultsFiles.neighbourData + ".txt"));
			// Write the data in the file with this format
			out.println("# Neighbour List for the initial state of " + experiment.nickname);
			out.println("# The value of the cutOff List is rl= " + experiment.cutoffList);
			for (int i = 0; i < s.particles.size(); i++) {
				out.println("# List of neighbors of atom " + (1 + i));
				out.println("# Position of atom: " + (1 + i) + ": " + s.particles.get(i).r.x + "     "
						+ s.particles.get(i).r.y + "     " + s.particles.get(i).r.z);
				out.println("# Number of Neighbors: " + s.particles.get(i).numberNeighbors);
				out.printf("# %-25s %45s %n", "Index", "Position of neighbors (reduced units)");
				for (Integer jneig : s.particles.get(i).vicinity) {
					out.printf(" %-25d %15f %15f %15f %n", jneig + 1, s.particles.get(jneig).r.x,
							s.particles.get(jneig).r.y, s.particles.get(jneig).r.z);
				}
				out.println();
			}
		} catch (IOException e) {
			System.out.println("Error al abrir " + resultsFiles.neighbourData);
			resultsFiles.neighbourData = null;
		} finally {
			if (out != null)
				// Close the file
				out.close();

		}
	}

	/**
	 * Method that makes a .txt file with the force between 2 particles 'r' (in its
	 * units), the value of the pair distribution function in that point and the
	 * value of the potential in the same 'r'. This .txt file contains 3 columns,
	 * one with the values of 'r', another with the values of the pair distribution
	 * function 'G(r)', and another with the values of the potential 'V(r)'.
	 * 
	 * @param potentialData a Graph with the distances between the atoms and the
	 *                      corresponding potential energy
	 * @param forceData     a Graph with the distances between the atoms and the
	 *                      corresponding force
	 */
	public void potentialAndForceTextFile(Graph potentialData, Graph forceData) {
		PrintWriter out = null;
		try {
			// Open the file
			out = new PrintWriter(new FileWriter(resultsFiles.potentialAndForce + ".txt"));
			// Write the data in the file with this format
			out.println(
					"# Data of the potential and magnitude of the force versus distance for " + experiment.nickname);
			out.println(
					"# The distances are values between the maximum and minimum distance between the particles of the first state");
			out.println("# The potential and the force are shifted to avoid discontinuity at the cutOff Radius rc= "
					+ experiment.cutoffRadius);
			out.println("# Epsilon and Sigma in reudced units are 1");
			out.printf("# %-25s %25s %17s %n", "Distance(sigma units)", "Potential Energy", "Force");

			for (int i = 0; i < potentialData.points.size(); i++) {
				out.printf("%-25f %20f %23f %n", potentialData.points.get(i).x, potentialData.points.get(i).y,
						forceData.points.get(i).y);
			}
		} catch (IOException e) {
			System.out.println("Error al abrir " + resultsFiles.potentialAndForce);
			resultsFiles.potentialAndForce = null;
		} finally {
			if (out != null)
				// Close the file
				out.close();

		}

	}

	/**
	 * Method that makes a .txt file with the distance between 2 particles 'r' (in
	 * its units), the value of the pair distribution function in that point and the
	 * value of the potential in the same 'r'. This .txt file contains 3 columns,
	 * one with the values of 'r', another with the values of the pair distribution
	 * function 'G(r)', and another with the values of the potential 'V(r)'.
	 * 
	 * @param pairDistributionData A 3 dimensional array. The first coordinate can
	 *                             be 0 (for the x points, the distances), 1 (for
	 *                             the y points, the values of 'G(r)'), or 2 (for
	 *                             the y points, the values of 'V(r)'). The second
	 *                             index refeers to each point.
	 */
	public void pairDistributionTextFile(Graph pairDistributionData) {
		PrintWriter out = null;
		try {
			// Open the file
			out = new PrintWriter(new FileWriter(resultsFiles.pairData + ".txt"));
			// Write the data in the file with this format
			out.println("# Pair distribution Function " + experiment.nickname);
			out.printf("# %-25s %20s %n", "Distance(reduced units)", "G(r)");

			for (int i = 0; i < pairDistributionData.points.size(); i++) {
				out.printf(" %-25f %20.0f %n", pairDistributionData.points.get(i).x,
						pairDistributionData.points.get(i).y);
			}
		} catch (IOException e) {
			System.out.println("Error al abrir " + resultsFiles.pairData);
			resultsFiles.pairData = null;
		} finally {
			if (out != null)
				// Close the file
				out.close();

		}

	}

	/**
	 * Method that makes a .txt file with the steps and with the potential, kinetic
	 * and total energy per atom (in its units). This .txt file contains 4 columns,
	 * one with the steps, another with the values of the potential energies,
	 * another with the values of the kinetic energies and the last one with the
	 * values of the total energy.
	 * 
	 * @param energyLog A List with the potential, kinetic and total energies for
	 *                  each step
	 */
	public void energyTextFile(List<EnergyLog> energyLog) {
		PrintWriter out = null;
		try {
			// Open the file
			out = new PrintWriter(new FileWriter(resultsFiles.energyData + ".txt"));
			// Write the data in the file with this format
			out.println("# Potential, Kinetic and Total energy per atom (divided by the total number of atoms) of the "
					+ experiment.nickname);
			out.println("# from the first step to the " + experiment.nSteps);
			out.println(
					"# The potential energy is calculated with the Lennard-Jones potential shiftted to avoid discontinuties");
			out.println("# The kinetic energy is the sum of 1/2v**2 of each particle");
			out.println("# The total energy is the sum of the kinetic and potential energy ");
			out.printf("# %-25s %25s %25s %25s %n", "Step", "Potential Energy per atom", "Kinetic Energy per atom",
					"Total Energy per atom");

			for (int i = 0; i < energyLog.size(); i++) {
				out.printf(" %-25d %20f %20f %30f %n", i + 1, energyLog.get(i).u / s.particles.size(),
						energyLog.get(i).k / s.particles.size(), energyLog.get(i).total / s.particles.size());
			}
		} catch (IOException e) {
			System.out.println("Error al abrir " + resultsFiles.energyData);
			resultsFiles.energyData = null;
		} finally {
			if (out != null)
				// Close the file
				out.close();

		}

	}

	/**
	 * This creates a velocity movie, with the particle's positions and the forces
	 * (proportional to accelerations) represented as arrows. It exports a text file
	 * which, latter on, can be read by xCryDen program (externally). The name is
	 * written in files, from ResultsFiles.
	 */
	public void positionMovie(List<State> movie) {
		double A=experiment.numberCells*lengthConverterA;
		PrintWriter out = null;
		try {
			// Open the file
			out = new PrintWriter(new FileWriter(resultsFiles.positionMovieData + ".axsf"));
			// Write the data in the file with this format
			out.println("ANIMSTEPS "+(movie.size()));
			out.println("CRYSTAL");
			out.println("PRIMVEC");
			out.printf("%-10f %-10f %-10f %n", experiment.latcon*experiment.numberCells/1E-10, 0.0, 0.0 );
			out.printf("%-10f %-10f %-10f %n", 0.0, experiment.latcon*experiment.numberCells/1E-10, 0.0 );
			out.printf("%-10f %-10f %-10f %n",  0.0, 0.0, experiment.latcon*experiment.numberCells/1E-10 );
			for (int i = 0; i < movie.size(); i++) {
				out.println("PRIMCOORD " + (1 + i));
				out.println(movie.get(0).particles.size()+" "+ 1);
				for (int j = 0; j < movie.get(i).particles.size(); j++) {
					out.printf(" %-25d %20.4f %20.4f %20.4f %n", experiment.atomicNumber, movie.get(i).particles.get(j).r.x*lengthConverterA,movie.get(i).particles.get(j).r.y*lengthConverterA,movie.get(i).particles.get(j).r.z*lengthConverterA);
				}
			}
		} catch (IOException e) {
			System.out.println("Error al abrir " + resultsFiles.positionMovieData);
			resultsFiles.positionMovieData = null;
		} finally {
			if (out != null)
				// Close the file
				out.close();

		}

	}

	

}
