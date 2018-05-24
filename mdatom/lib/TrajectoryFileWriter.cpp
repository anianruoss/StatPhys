#include "TrajectoryFileWriter.h"
#include "BinaryIO.h"
#include "CoordinatesAndVelocitiesInitializer.h" // For MAXTITLE value
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <utility>

using namespace std;

TrajectoryFileWriter::TrajectoryFileWriter(const MDParameters &parameters,
                                           std::string finalCoordFilename,
                                           std::string trajFilename)
  : par(parameters),
    finalCoordinatesFilename(std::move(finalCoordFilename)),
    trajectoryCoordinatesFilename(std::move(trajFilename)) {
}

/* xyz spec:
 * per step:
 * 1 line n_atoms
 * 1 line comment
 * per atom (1 per line):
 * Symbol, x, y, z: separated by whitespace
 */

void TrajectoryFileWriter::writeBeforeRun() {
    ofstream fout1; // trajectory output
    if (par.trajectoryOutput) {
        fout1.open(trajectoryCoordinatesFilename, ios::out);
        if (fout1.bad()) {
            throw std::runtime_error("can't open " + trajectoryCoordinatesFilename);
        }

        return;
        fout1 << par.title << endl;
    }
}

void TrajectoryFileWriter::writeFinalCoordinates(const std::vector<double>& positions,
                                                 const std::vector<double>& velocities) {
    if (par.finalXVOutput == FinalCoordinateFileFormat::ascii) {
        writeFinalCoordinatesInAsciiForm(positions, velocities);
    }
    else {
        writeFinalCoordinatesInBinaryForm(positions, velocities);
    }
}

void TrajectoryFileWriter::writeFinalCoordinatesInBinaryForm(const std::vector<double>& positions,
                                                             const std::vector<double>& velocities) {
    ofstream fout2;
    fout2.open(finalCoordinatesFilename, ios::out | ios::binary);
    if (fout2.bad()) {
        throw std::runtime_error("can't open " + finalCoordinatesFilename);
    }
    fout2.write(par.title.c_str(), MAXTITLE);
    BinaryIO::write(fout2, positions);
    BinaryIO::write(fout2, velocities);
}

void TrajectoryFileWriter::writeFinalCoordinatesInAsciiForm(const std::vector<double>& positions,
                                                            const std::vector<double>& velocities) {
    ofstream fout2;
    fout2.open(finalCoordinatesFilename, ios::out);
    if (fout2.bad()) {
        throw std::runtime_error("can't open " + finalCoordinatesFilename);
    }
    fout2 << par.title << "\n" ;
    fout2 << par.numberAtoms << "\n" ;
    for (int j = 0; j < par.numberAtoms; j++) {
        fout2 << setw(6) << j;
        for (int m = 0; m < 3; m++) {
            fout2 << setw(15) << positions[3 * j + m];
        }
        for (int m = 0; m < 3; m++) {
            fout2 << setw(15) << velocities[3 * j + m];
        }
        fout2 << "\n";
    }
}

void TrajectoryFileWriter::writeOutTrajectoryStep(const std::vector<double>& positions) {
    if (par.trajectoryOutput) {
        if (par.trajectoryOutputFormat == TrajectoryFileFormat::binary) {
            writeOutTrajectoryStepInBinaryForm(positions);
        } else if (par.trajectoryOutputFormat == TrajectoryFileFormat::ascii) {
            writeOutTrajectoryStepInAsciiForm(positions);
        }
    }
}

void TrajectoryFileWriter::writeOutTrajectoryStepInBinaryForm(const std::vector<double>& positions) {
    ofstream fileBW;
    fileBW.open(trajectoryCoordinatesFilename, ios::out | ios::app | ios::binary);
    if (fileBW.bad()) {
        throw runtime_error("I/O ERROR: cannot write to file: " + trajectoryCoordinatesFilename);
    }
    BinaryIO::write(fileBW, positions);
}

void TrajectoryFileWriter::writeOutTrajectoryStepInAsciiForm(const std::vector<double>& positions) {
	// This function writes Molfiles as defined here:
	// https://sourceforge.net/p/jmol/code/HEAD/tree/trunk/Jmol-datafiles/mol/ctfile.pdf
	
    const string atom = "Ar";
    ofstream fileFW;
    fileFW.open(trajectoryCoordinatesFilename, ios::out | ios::app);
    if (fileFW.bad()) {
        throw runtime_error("I/O ERROR: cannot write to file: " + trajectoryCoordinatesFilename);
    }
	// Header block
	fileFW << par.title << endl << endl << endl;
	// Bonds are between subsequent atoms
	// Counts line
	int N = par.numberAtoms;
	fileFW << setw(3) << N << setw(3) << N-1 << setw(3) << 0 << setw(3) << 0
		<< setw(3) << 0 << std::endl;
	fileFW << fixed << setprecision(4);
	// Atom block
    for (int i = 0; i < par.numberAtoms; i++){
        for (int c = 0; c < 3; c++){
			fileFW << setw(10) << positions[i*3+c];
        }
		fileFW << setw(3) << atom;
        fileFW << endl;
    }
	// Bond block
	for (int i = 1; i < N; i++){
		fileFW << setw(3) << i << setw(3) << i + 1 << setw(3) << 1 << endl;
	}
	fileFW << "M  END" << endl;
}
