#include "MDRun.h"
#include "PeriodicBoundaryConditions.h"
#include "CenterOfMassCalculator.h"
#include "TrajectoryFileWriter.h"
#include <cmath>

MDRun::MDRun(const MDParameters& parameters, MDRunOutput& out, TrajectoryFileWriter& trajectoryFileWriter)
  : par(parameters),
    output(out),
    trajectoryWriter(trajectoryFileWriter),
    forceCalculator(parameters),
    radialDistribution(parameters.numberRadialDistrPoints, parameters.radialDistrCutoffRadius) {
}

void MDRun::run(std::vector<double> &x, std::vector<double> &v) {
    forces.resize(x.size());
    synchronizedPositions.resize(x.size());
    radialDistribution.setZero();

    initializeVariables();
    initializeTemperature(v);

    output.printInitialTemperature(properties[1] / fac);
    output.printIterationStart();

    /* dynamics step */
    double time = par.initialTime;
    for (int nstep = 0; nstep < par.numberMDSteps; nstep++) {
        time += par.timeStep;
        performStep(x, v, nstep, time);
    }

    printAverages(time);
}

void MDRun::initializeVariables() {
    nat3 = 3 * par.numberAtoms;
    const double boltzmannConstant = 8.3144598e-3; // units: K^-1 ps^-2 u nm^2
    fac = nat3 * boltzmannConstant / 2.;
    ekin0 = fac * par.targetTemperature;
    halfTimeStep = par.timeStep / 2;
    dtm = par.timeStep / par.atomicMass;
    vol = par.boxSize[0] * par.boxSize[1] * par.boxSize[2];

    nhpr = 100 * par.propertyPrintingInterval;
    nlsq = par.numberMDSteps / 10;
    if (nlsq < 10) {
        nlsq = 10;
    }
    if (nlsq < par.propertyPrintingInterval) {
        nlsq = par.propertyPrintingInterval;
    }
    if (nhpr > nlsq) {
        nhpr = nlsq;
    }
    for (int i = 0; i < numberProperties; i++) {
        properties[i] = 0.;
        averages[i] = 0.;
        fluctuations[i] = 0.;
    }
}

void MDRun::initializeTemperature(const std::vector<double>& velocities) {
    double kineticEnergy = 0;
    for (int j3 = 0; j3 < nat3; j3++) {
        kineticEnergy += velocities[j3] * velocities[j3];
    }
    kineticEnergy *= (par.atomicMass / 2.);
    properties[1] = kineticEnergy;
    if (par.mdType == SimulationType::constantTemperature) {
        if (kineticEnergy < 1.e-6) {
            ekg = ekin0;
        }
        else {
            ekg = kineticEnergy;
        }
    }
}

void MDRun::performStep(std::vector<double>& positions, std::vector<double>& velocities, int nstep, double time) {
    /* put atoms in central periodic box */
    PeriodicBoundaryConditions::recenterAtoms(par.numberAtoms, positions, par.boxSize);

    /* calculate forces, potential energy, virial
     * and contribution to the radial distribution function
     */
    forceCalculator.calculate(positions, forces);
    radialDistribution.addInstantaneousDistribution(forceCalculator.getInstantaneousRadialDistribution());
    double vir = forceCalculator.getVirial();
    properties[2] = forceCalculator.getPotentialEnergy();
    properties[3] = vir;

    /* determine velocity scaling factor, when coupling to a bath */
    double scal = 1;
    if (par.mdType == SimulationType::constantTemperature) {
        double dtt = par.timeStep / par.temperatureCouplingTime;
        scal = std::sqrt(1 + dtt * (ekin0 / ekg - 1));
    }


    /* perform leap-frog integration step,
     * calculate kinetic energy at time t-dt/2 and at time t,
     * and calculate pressure
     */
    double oldKineticEnergy = 0.;
    double newKineticEnergy = 0.;

	/* Store current positions for enforcing constraints */
	std::vector<double> prev_pos = positions;
    for (int j3 = 0; j3 < nat3; j3++) {
        double oldVelocity = velocities[j3];
        double newVelocity = (oldVelocity + forces[j3] * dtm) * scal;
        oldKineticEnergy += newVelocity * newVelocity;
        newKineticEnergy += (oldVelocity + newVelocity) * (oldVelocity + newVelocity);
        velocities[j3] = newVelocity;
        positions[j3] += newVelocity * par.timeStep;
	}
	/* # Enforce constraints
	 * _comments: Markdown with inline TeX_
	 *
	 * $f_c = \frac{\mu}{2 \Delta t^2} \frac{d^2 - d'^2}{d' \cdot d} \cdot d $
	 * where 
	 *
	 * * $d$ is the constrained bond vector before integration
	 * * $d'$ is the bond vector after integration step (not fully constrained)
	 * * $\mu$ is the reduced mass ($\equiv \frac{m_a}{2}$ here)
	 *
	 * The loop body consists of:
	 *
	 * 1. Compute d, d' from prev_pos and positions
	 * 2. Compute f from d, d'
	 * 3. Compute $\Delta r_i$
	 */

    bool constrained = false;

	int shake_it = 0;
    while (!constrained) {
		std::cout << "\nshake_it :" << shake_it << " ";
        std::vector<double> delta_pos(positions.size(), 0.);
        constrained = true;

        for (int atom = 0; atom < par.numberAtoms-1; ++atom) {
            double d_start[3], d_end[3], f_c[3];
            double ds, de, dotP;
            ds = de = dotP = 0.;

            for (int i = 0; i < 3; ++i) {
                d_start[i] = prev_pos[(atom+1)*3 + i] - prev_pos[atom*3 + i];
                d_end[i] = positions[(atom+1)*3 + i] - positions[atom*3 + i];

                ds += d_start[i] * d_start[i];
                de += d_end[i] * d_end[i];

                dotP += d_start[i] * d_end[i];
            }

            ds = std::sqrt(ds);
            de = std::sqrt(de);

			double error = (std::abs(ds - de) / ds);
            if (error >= shake_rel_tol)
                constrained = false;

			// Note that we are actually iterating over constraints, not atoms
			// Note that the output happens before computing constraint forces
			std::cout << " atom " << atom << " err " << error << '\t';

            for (int i = 0; i < 3; ++i)
                f_c[i] = d_start[i] * mu_d_2_tsq * (ds*ds - de*de) / dotP;

            for (int i = 0; i < 3; ++i) {
                // delta_pos[atom*3 + i] += Dt2_d_m * f_c[i];
                delta_pos[atom*3 + i] -= Dt2_d_m * f_c[i];
                delta_pos[(atom+1)*3 + i] += Dt2_d_m * f_c[i];
                // delta_pos[(atom+1)*3 + i] -= Dt2_d_m * f_c[i];
            }
        }
		shake_it ++;
		if (shake_it > 200){
			std::cerr << "no convergence" << std::endl;
			break;
		}

        for (int j3 = 0; j3 < nat3; ++j3)
            positions[j3] = positions[j3] + delta_pos[j3];
    }

    /*

    int max_shake_it = 20;
	bool constrained = false;  // true iff all constraints are satisifed (i.e. within tolerance)

	for (int shake_it = 0; shake_it < max_shake_it and not constrained; shake_it++){
		// per shake iteration
		std::vector<double> new_positions(positions.size());
		constrained = true;
		for (int atom = 0; atom < par.numberAtoms; atom ++){
			double f_c[3] = {0, 0, 0};
			// Compute d, d'
			for (int neighbor: {atom - 1, atom + 1}){
				 * // Wrap around instead of continue
				 * if (neighbor < 0)
				 * 	neighbor += par.numberAtoms;
				 * else if (neighbor >= par.numberAtoms)
				 * 	neighbor -= par.numberAtoms;
				// Don't wrap around
				if (neighbor < 0 or neighbor >= par.numberAtoms)
					continue;
				// Compute d', d
				double d_a[3], d[3];
				for (int i = 0; i < 3; i++){
					d_a[i] = positions[neighbor * 3 + i] - positions[atom * 3 + i];
					d[i] = prev_pos[neighbor * 3 + i] - prev_pos[atom * 3 + i];
				}

				double diff = 0, d_abs = 0, error;
				for(int i = 0; i < 3; i++){
					diff += std::pow(d[i] - d_a[i], 2);
					d_abs += std::pow(d[i], 2);
				}
				error = std::sqrt(diff)/std::sqrt(d_abs);
				if (error > shake_rel_tol)
					constrained = false;

				// Compute squares of d', d and d * d'
				double d_a_2 = 0, d_2 = 0, d_d_a = 0;
				for (int i = 0; i < 3; i++){
					d_a_2 += std::pow(d_a[i], 2);
					d_2 += std::pow(d[i], 2);
					d_d_a += d[i] * d_a[i];
				}

				// Compute $f_c$ (constraint force)
				double factor = mu_d_2_tsq * (d_2 - d_a_2) / d_d_a;
				for (int i = 0; i < 3; i++)
					f_c[i] += d[i] * factor;
			}
			// Per atom: Compute $\Delta r_i$
			double r_i[3];
			for (int j = 0; j < 3; j++)
				r_i[j] = Dt2_d_m * f_c[j];
			// Update new position of this atom
			for (int j = 0; j < 3; j++)
				new_positions[atom * 3 + j] = prev_pos[atom * 3 + j] + r_i[j];
        }

		positions = new_positions;
        if (shake_it + 1 == max_shake_it)
			std::cerr << "Maximum number of iterations reached!" << std::endl;
	}
    */

    oldKineticEnergy *= (par.atomicMass / 2.);
    newKineticEnergy *= (par.atomicMass / 8.);
    properties[1] = newKineticEnergy;
    properties[0] = properties[1] + properties[2];
    double pres = 2. * (newKineticEnergy - vir) / (vol * 3.);
    properties[4] = pres;
    properties[5] = scal;
    if (par.mdType == SimulationType::constantTemperature) {
        ekg = oldKineticEnergy;
    }

    /* update arrays for averages and fluctuations */
    for (int m = 0; m < numberProperties; m++) {
        averages[m] += properties[m];
        fluctuations[m] += properties[m] * properties[m];
    }

    printOutputForStep(positions, velocities, nstep, time);
}

void MDRun::printOutputForStep(const std::vector<double> &positions, const std::vector<double> &velocities, int nstep, double time) {
    if ((nstep + 1) == (nstep + 1) / par.trajectoryOutputInterval * par.trajectoryOutputInterval) {
        trajectoryWriter.writeOutTrajectoryStep(positions);
    }

    if (nstep == (nstep + 1) / nhpr * nhpr) {
        output.printPropertiesHeader();
    }

    if ((nstep + 1) == (nstep + 1) / par.propertyPrintingInterval * par.propertyPrintingInterval || nstep == 0) {
        output.printProperties(nstep, time, properties);
    }

    /* calculate and print center of mass motion
     * once in nlsq steps, at time t-dt/2
     * The positions must be back-calculated for t-dt/2, because of the time shift between x and v (leap-frog)
     */
    if ((nstep + 1) == (nstep + 1) / nlsq * nlsq) {
        for (int j3 = 0; j3 < nat3; j3++) {
            synchronizedPositions[j3] = positions[j3] - velocities[j3] * halfTimeStep;
        }
        CenterOfMassCalculator cm;
        cm.update(par.numberAtoms, synchronizedPositions, velocities, par.atomicMass);
        cm.printResults(output);
    }
}

void MDRun::printAverages(double time) {
    double tspan = par.numberMDSteps;
    for (int m = 0; m < numberProperties; m++) {
        averages[m] = averages[m] / tspan;
        fluctuations[m] = std::sqrt(std::abs(fluctuations[m] / tspan - averages[m] * averages[m]));
    }
    output.printAverages(par.numberMDSteps, time, averages);
    output.printRMSFluctuations(par.numberMDSteps, time, fluctuations);
    output.printAverageAndRMSTemperature(averages[1] / fac, fluctuations[1] / fac);
}

const AveragedRadialDistribution& MDRun::getRadialDistribution() const {
    return radialDistribution;
}
