#ifndef TRAJECTORYFILEWRITER_H
#define TRAJECTORYFILEWRITER_H

#include "MDParameters.h"
#include <string>
#include <vector>

/*!
 * This class is responsible for writing out the atom coordinates during a MD simulation.
 */
class TrajectoryFileWriter {
  public:
    TrajectoryFileWriter(const MDParameters& parameters,
                         std::string finalCoordFilename,
                         std::string trajFilename);
    void writeBeforeRun();
    void writeAfterRun();
    void writeFinalCoordinates(const std::vector<double>& positions, const std::vector<double>& velocities);
    void writeOutTrajectoryStep(const std::vector<double>& positions);

  private:
    void writeFinalCoordinatesInBinaryForm(const std::vector<double>& positions, const std::vector<double>& velocities);
    void writeFinalCoordinatesInAsciiForm(const std::vector<double>& positions, const std::vector<double>& velocities);
    void writeOutTrajectoryStepInBinaryForm(const std::vector<double>& positions);
	// Ascii form = PDB File
    void writeOutTrajectoryStepInAsciiForm(const std::vector<double>& positions);

    const MDParameters& par;
    const std::string finalCoordinatesFilename;
    const std::string trajectoryCoordinatesFilename;
	int model = 0; // n-th model for pdb
};

#endif // TRAJECTORYFILEWRITER_H
