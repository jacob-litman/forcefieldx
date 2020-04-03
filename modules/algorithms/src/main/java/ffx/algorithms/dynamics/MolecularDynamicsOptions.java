//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.dynamics;

import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.dynamics.integrators.*;
import ffx.algorithms.dynamics.thermostats.*;
import ffx.algorithms.mc.MonteCarloListener;
import ffx.numerics.Constraint;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.MolecularAssembly;
import ffx.potential.cli.WriteoutOptions;
import ffx.utilities.Constants;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

/**
 * Options for molecular dynamics, intended to make many MD flags
 * immutable without requiring excessively long constructors.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class MolecularDynamicsOptions {
    private static final Logger logger = Logger.getLogger(MolecularDynamicsOptions.class.getName());

    final MolecularAssembly[] assemblies;
    final CompositeConfiguration properties;
    final IntegratorEnum iEnum;
    final Integrator integrator;
    final ThermostatEnum tEnum;
    final Thermostat thermostat;
    final Potential potential;
    /**
     * Timestep in picoseconds.
     */
    public final double dt;
    final int numSnapshotsToKeep;

    // Should only ever be accessed by the created MD object!
    final double[] x;
    final double[] v;
    final double[] a;
    final double[] aPrevious;
    final double[] mass;
    
    MonteCarloListener mcListener = null;
    Barostat barostat;
    public MolecularDynamics.VerbosityLevel vLevel = MolecularDynamics.VerbosityLevel.DEFAULT_VERBOSITY;
    File dynFile;
    public boolean initVelocities;
    boolean automaticWriteouts;
    int dynSleepTime;
    final double restartInterval;
    final double snapshotInterval;
    final double logInterval;
    boolean saveSnapshotAsPDB;
    
    public MolecularDynamicsOptions(MolecularAssembly[] assemblies, Potential potential,
                                    DynamicsOptions dynOpts, WriteoutOptions writeout,
                                    CompositeConfiguration properties) {
        this.assemblies = Arrays.copyOf(assemblies, assemblies.length);
        this.properties = properties;
        this.potential = potential;
        this.barostat = potential instanceof Barostat ? (Barostat) potential : null;
        IntegratorEnum reqI = dynOpts.integrator;
        ThermostatEnum reqT = dynOpts.thermostat;
        double targetTemperature = dynOpts.getTemp();
        dt = dynOpts.getDt() * Constants.FSEC_TO_PSEC;
        numSnapshotsToKeep = properties.getInteger("dynamicsSnapshotMemory", 0);
        dynSleepTime = properties.getInt("dynamics-sleep-nanos", MolecularDynamics.DEFAULT_DYNAMICS_SLEEP_TIME);
        List<Constraint> constraints = potential.getConstraints();

        int numberOfVariables = potential.getNumberOfVariables();
        x = new double[numberOfVariables];
        potential.getCoordinates(x);
        v = new double[numberOfVariables];
        a = new double[numberOfVariables];
        aPrevious = new double[numberOfVariables];
        mass = potential.getMass();

        restartInterval = dynOpts.getCheckpoint();
        snapshotInterval = dynOpts.getSnapshotInterval();
        logInterval = dynOpts.getReport();

        saveSnapshotAsPDB = writeout.getFileType().equalsIgnoreCase("PDB");

        MolecularAssembly molecularAssembly = assemblies[0];

        // TODO: Better way of doing this.
        boolean oMMLogging = potential instanceof ForceFieldEnergyOpenMM;

        switch (reqI) {
            case RESPA:
                Respa respa = new Respa(numberOfVariables, x, v, a, aPrevious, mass);
                int in = molecularAssembly.getProperties().getInt("respa-dt", 4);
                if (in < 2) {
                    in = 2;
                }
                if (!oMMLogging) {
                    respa.setInnerTimeSteps(in);
                }
                logger.log(Level.FINE, format(" Created a RESPA integrator with %d inner time steps.", in));
                integrator = respa;
                break;
            case STOCHASTIC:
                double friction = properties.getDouble("friction", 91.0);
                logger.log(Level.INFO, format(" Friction set at %.3f collisions/picosecond", friction));

                Stochastic stochastic = new Stochastic(friction, numberOfVariables, x, v, a, mass);
                if (properties.containsKey("randomseed")) {
                    stochastic.setRandomSeed(properties.getInt("randomseed", 0));
                }
                integrator = stochastic;
                // The stochastic dynamics integration procedure will thermostat
                // the system. The ADIABTIC thermostat just serves to report the
                // temperature and initialize velocities if necessary.
                reqT = ThermostatEnum.ADIABATIC;
                break;
            case BEEMAN:
                integrator = new BetterBeeman(numberOfVariables, x, v, a, aPrevious, mass);
                break;
            case VERLET:
            case VELOCITYVERLET:
            default:
                integrator = new VelocityVerlet(numberOfVariables, x, v, a, mass);
        }
        iEnum = reqI;
        integrator.addConstraints(constraints);

        switch (reqT) {
            case BERENDSEN:
                double tau = properties.getDouble("tau-temperature", 0.2);
                thermostat = new Berendsen(numberOfVariables, x, v, mass, potential.getVariableTypes(), targetTemperature, tau, constraints);
                break;
            case BUSSI:
                tau = properties.getDouble("tau-temperature", 0.2);
                thermostat = new Bussi(numberOfVariables, x, v, mass, potential.getVariableTypes(), targetTemperature, tau, constraints);
                if (properties.containsKey("randomseed")) {
                    thermostat.setRandomSeed(properties.getInt("randomseed", 0));
                }
                break;
            case ADIABATIC:
            default:
                thermostat = new Adiabatic(numberOfVariables, x, v, mass, potential.getVariableTypes(), constraints);
                break;
        }
        tEnum = reqT;

        if (properties.containsKey("randomseed")) {
            thermostat.setRandomSeed(properties.getInt("randomseed", 0));
        }

        // For Stochastic dynamics, center of mass motion will not be removed.
        if (integrator instanceof Stochastic) {
            thermostat.setRemoveCenterOfMassMotion(false);
        }
    }
}
