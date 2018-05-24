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

void TrajectoryFileWriter::writeAfterRun() {
    ofstream fileFW;
    fileFW.open(trajectoryCoordinatesFilename, ios::out | ios::app);
    if (fileFW.bad()) {
        throw runtime_error("I/O ERROR: cannot write to file: " + trajectoryCoordinatesFilename);
    }
	// Bonds
	int N = par.numberAtoms;
	for (int step = 0; step < model; step++){
		for (int atom = 1; atom < N; atom++){
			fileFW << "CONECT" << setw(5) << step*N + atom << setw(5) << step*N + atom + 1 << endl;
			fileFW << "CONECT" << setw(5) << step*N + atom + 1 << setw(5) << step*N + atom << endl;
		}
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
	// This function writes PDB Files as defined here:
	// It writes the minimum amount of information necessary
	
    const string atom = "AR";
    ofstream fileFW;
    fileFW.open(trajectoryCoordinatesFilename, ios::out | ios::app);
    if (fileFW.bad()) {
        throw runtime_error("I/O ERROR: cannot write to file: " + trajectoryCoordinatesFilename);
    }
	// Header block
	// Bonds are between subsequent atoms
	// Counts line
	int N = par.numberAtoms;
	fileFW << "MODEL" << setw(5) << ' ' << model+1 << endl;
    for (int i = 0; i < par.numberAtoms; i++){
		fileFW << fixed << setprecision(3);
		fileFW << "HETATM" << setw(5) << model * N + i + 1<< " " << setw(4) << i + 1<< " UNK";
		fileFW << setw(6) << 1 << setw(4) << " ";
        for (int c = 0; c < 3; c++){
			// angström used
			// fileFW << setw(8) << 10 * positions[i*3+c];
			// perhaps not
			fileFW << setw(8) << positions[i*3+c];
        }
		fileFW << setw(6) << "  1.00" << setw(6) << "  0.00" << setw(10) << ' ';
		fileFW << right << setw(2) << atom << "  " << endl;
    }
	fileFW << "ENDMDL" << endl;
	model ++;
}
