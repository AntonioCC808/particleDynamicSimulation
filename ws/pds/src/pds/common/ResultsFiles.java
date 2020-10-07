
package pds.common;

public class ResultsFiles {
	
	/**
	 * Name of the text file with the initial velocities distribution histogram data
	 */
	public String positionsData;
	
	/**
	 * Name of the text file with the initial velocities on each direction
	 */
	public String velocitiesData;
	
	/**
	 * Name of the text file with the initial velocities distribution histogram data 
	 */
	public String velocitiesHistogramData;


	/**
	 * Name of the text file with the potential and the pair distribution function
	 * listed together
	 */
	public String potentialAndForce;

	/**
	 * Name of the text file with the potential, kinetic and total energies listed
	 * for all the steps of the simulation.
	 */
	public String energyData;
	
	/**
	 * Name of the text file which contains the neigbour list for the first state,
	 *  with the indexes of the neighbors for each particle and the total number of neigbours
	 */
	public String neighbourData;
	

	/**
	 * Name of the text file which contains the pair distribution function data list for the last state,
	 *  with the frequency for each distance distance
	 */
	public String pairData;

	/**
	 * This contains the name of the text file to be readen by xCryDen to generate
	 * the position and force movie.
	 */
	public String positionMovieData;


	

	/**
	 * This constructor makes the ResultsFiles, and saves it's names using the given
	 * prefix and suffix. The extensions (.txt, .jpg ...) are not included yet.
	 * 
	 * @param prefix This is the name of the directory (absolute), and includes the
	 *               last "\" symbol for the folder, as an string.
	 * @param suffix The name to be added to each file (experiment name). Do not
	 *               include any special character.
	 */
	public ResultsFiles(String prefix, String suffix) {
		/*
		 * The names of the files are composed of the prefix, a name specifying what it
		 * is (that is included here).
		 */
		positionsData = prefix + "positionsData_" + suffix;
		velocitiesData = prefix + "velocitiesData_" + suffix;
		velocitiesHistogramData = prefix + "velocitiesHistogramData_" + suffix;
		neighbourData = prefix + "neighbourData_" + suffix;
		potentialAndForce = prefix + "potentialAndForceData_" + suffix;
		energyData = prefix + "energyData_" + suffix;
		pairData = prefix + "pairData_" + suffix;
		positionMovieData = prefix + "positionMovieData_" + suffix;
	}

}
