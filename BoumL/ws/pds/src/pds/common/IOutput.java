
package pds.common;

/**
 * The interface for the output class, with the methods to represent each of the results.
 */
public interface IOutput {
  /**
   * Returns the results of a simulation done for the given experiment. The 
   * operation produces the output results files and returns the names of those files. 
   * 
   * @param experiment is the experiment object with the data needed to retrive
   * the simulation results
   * @return returns the ResultsFile object with the names of all the requested 
   * resulting images and data files generated. If any of those files could not 
   * be generated the corresponding string would be a null reference.
   */
  ResultsFiles results();

}
