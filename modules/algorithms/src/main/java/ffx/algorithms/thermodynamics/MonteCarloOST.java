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
package ffx.algorithms.thermodynamics;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.lang.System.nanoTime;

import ffx.algorithms.dynamics.DynamicAlgorithm;
import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.numerics.math.RunningStatistics;
import ffx.utilities.MDListener;
import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.algorithms.cli.DynamicsOptions;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.algorithms.mc.BoltzmannMC;
import ffx.algorithms.mc.LambdaMove;
import ffx.algorithms.mc.MDMove;
import ffx.algorithms.thermodynamics.OrthogonalSpaceTempering.Histogram;
import ffx.numerics.Potential;

import javax.annotation.Nullable;

import static ffx.utilities.Constants.NS2SEC;

/**
 * Sample a thermodynamic path using the OST method, with the time-dependent
 * bias built up using Metropolis Monte Carlo steps.
 * <p>
 * The algorithm generates coordinate (X) MC moves using molecular dynamics at a
 * fixed lambda value (i.e. using OpenMM), followed by MC lambda moves.
 * <p>
 * 1.) At a fixed Lambda, run a defined length MD trajectory to "move"
 * coordinates and dU/dL on an approximate potential U* (i.e. no OST Bias).
 * <p>
 * 2.) Accept / Reject the MD move with probability exp[-Beta(dU - dU*)] where
 * dU is the change in AMOEBA + Bias energy and dU* is the change in AMOEBA +
 * Kinetic energy from the MD.
 * <p>
 * 3.) Randomly change the value of Lambda.
 * <p>
 * 4.) Accept / Reject the Lambda move using the AMOEBA + OST Bias energy.
 * <p>
 * 5.) Add to the time dependent 2D bias using the current values of Lambda and
 * dU/dL.
 *
 * @author Michael J. Schnieders
 * @author Hernan Beranbe
 * @author Mallory R. Tollefson
 * @author Jacob Litman
 * @since 1.0
 */
public class MonteCarloOST extends BoltzmannMC implements DynamicAlgorithm {

    /**
     * Logger object to print out information for this class.
     */
    private static final Logger logger = Logger.getLogger(MonteCarloOST.class.getName());
    /**
     * Energy conservation during MD moves should generally be within ~0.1
     * kcal/mol. A change in total energy of 1.0 kcal/mol or more is of
     * significant concern that the time step is too large, or lambda moves are
     * too aggressive.
     */
    private static final double ENERGY_CONSERVATION_TOLERANCE = 10.0;
    /**
     * Potential object used to retrieve the coordinates for the system.
     */
    private final Potential potential;
    /**
     * OST object used to retrieve OST energy throughout the simulation.
     */
    private final OrthogonalSpaceTempering orthogonalSpaceTempering;
    /**
     * MDMove object for completing MC-OST molecular dynamics moves.
     */
    private final MDMove mdMove;
    /**
     * Molecular Dynamics object encapsulated in the MDMove.
     */
    private final MolecularDynamics molDyn;
    /**
     * Number of steps to take per MC-OST round.
     */
    private final long stepsPerMove;
    /**
     * Lambda move object for completing MC-OST lambda moves.
     */
    private final LambdaMove lambdaMove;
    /**
     * Double that keeps track of our lambda value.
     */
    private double lambda = 1.0;
    /**
     * Boolean that tells algorithm that we are in the equilibration phase of MC-OST.
     */
    private boolean equilibration = false;

    /**
     * Controls the effect of verbose by logging at FINE vs. INFO.
     */
    private final Level verboseLoggingLevel;
    /**
     * Deposit a bias once every N MC cycles. Defaults to 1.
     */
    private final int biasDepositionFrequency;
    /**
     * True if MC-OST is the top-level algorithm
     */
    private final boolean topLevelAlgorithm;

    /**
     * Current logging level to use.
     */
    private Level currentLevel = Level.INFO;
    private final Level regularLevel;
    private final Level verboseLevel;
    
    private final double[] proposedCoordinates;

    private final List<MCActions> sinceLastMetropolis = new ArrayList<>(3);
    private final Deque<Runnable> reverts = new ArrayDeque<>(2);

    private double initFFE;
    private double initKE;
    private double initBias;
    // Initial total PE is initFFE + initBias.
    private double currentFFE;
    private double currentKE;
    private double currentBias;
    
    private RunningStatistics acceptXYZ;
    private RunningStatistics acceptLambda;
    private RunningStatistics acceptCombined;

    private static final String equilString = "Equilibration";
    private static final String prodString = "Production";
    private static final String equilAlgoString = "";
    private static final String algoFormatString = " %s %s MOST sampling with lambda and HMC moves";

    /**
     * <p>
     * Constructor for MonteCarloOST.</p>
     *
     * @param opts Options to apply to this MC-OST.
     */
    public MonteCarloOST(MonteCarloOstOptions opts, boolean topLevel) {
        this.potential = opts.potentialEnergy;
        this.orthogonalSpaceTempering = opts.orthogonalSpaceTempering;
        this.verboseLoggingLevel = opts.verbose ? Level.INFO : Level.FINE;
        this.molDyn = opts.moldyn;
        this.stepsPerMove = opts.cycleLength;
        this.topLevelAlgorithm = topLevel;
        this.regularLevel = opts.regularLevel;
        this.verboseLevel = opts.verboseLevel;

        // Create the MC MD and Lambda moves.
        mdMove = new MDMove(potential, opts.moldyn, stepsPerMove, opts.temperature);
        boolean discreteLambda = opts.properties.getBoolean("discrete-lambda", false);

        if (discreteLambda) {
            double expected = orthogonalSpaceTempering.getHistogram().dL;
            if (opts.moveSize != expected) {
                logger.warning(String.format(" Passed in discrete move size %.5f, " +
                        "but expected %.5f by inspecting histogram!", opts.moveSize, expected));
                opts.moveSize = expected;
            }
        }

        if (opts.properties.containsKey("randomseed")) {
            int randomSeed = opts.properties.getInt("randomseed", 0);
            logger.info(format(" Setting random seed for lambdaMove to %d ", randomSeed));
            lambdaMove = new LambdaMove(orthogonalSpaceTempering, discreteLambda, opts.moveSize, randomSeed);
            setRandomSeed(randomSeed);
        } else {
            lambdaMove = new LambdaMove(orthogonalSpaceTempering, discreteLambda, opts.moveSize, ThreadLocalRandom.current().nextLong());
        }

        // Changing the value of lambda will be handled by this class, as well as adding the time dependent bias.
        orthogonalSpaceTempering.setPropagateLambda(false);
        biasDepositionFrequency = opts.biasDepositionFrequency;
        if (biasDepositionFrequency < 1) {
            throw new IllegalArgumentException("The property mc-ost-biasf must be a positive integer, found " + biasDepositionFrequency + " !");
        } else if (biasDepositionFrequency > 1) {
            logger.info(format(" MC-OST will deposit a bias only once per %d MC cycles (mc-ost-biasf).", biasDepositionFrequency));
        }
        
        proposedCoordinates = new double[potential.getNumberOfVariables()];
    }

    public static class MonteCarloOstOptions {
        // Constructor variables.
        final Potential potentialEnergy;
        final OrthogonalSpaceTempering orthogonalSpaceTempering;
        final DynamicsOptions dynamics;
        final CompositeConfiguration properties;
        final MolecularDynamics moldyn;
        final double temperature;
        final long cycleLength;

        // Property variables.
        final int biasDepositionFrequency;

        // Settable by exterior methods.
        public boolean verbose = false;
        public double moveSize;
        Level regularLevel = Level.INFO;
        Level verboseLevel = Level.INFO;

        public MonteCarloOstOptions(Potential p, OrthogonalSpaceTempering ost,
                                    DynamicsOptions dyn, CompositeConfiguration properties,
                                    MolecularDynamics md, long cycleLength) {
            this.potentialEnergy = p;
            this.orthogonalSpaceTempering = ost;
            this.dynamics = dyn;
            this.properties = properties;
            this.moldyn = md;
            this.cycleLength = cycleLength;

            this.temperature = dynamics.getTemp();

            biasDepositionFrequency = properties.getInt("mc-ost-biasf", 1);

            // TODO: Set thermostat to ADIABATIC upstream of here, before the MD is created.
            ThermostatEnum tstat = dynamics.thermostat;
            if (!tstat.equals(ThermostatEnum.ADIABATIC)) {
                throw new IllegalArgumentException(format("MC-OST requires the ADIABATIC thermostat, found %s.", tstat));
            }

            IntegratorEnum integ = dynamics.integrator;
            if (!integ.knownReversible || !integ.knownDeterministic) {
                throw new IllegalArgumentException(format("MC-OST requires " +
                        "a reversible deterministic integrator (e.g. VERLET, RESPA), found %s!", integ));
            }
        }
    }

    public MolecularDynamics getMD() {
        return mdMove.getMD();
    }

    /**
     * Sets the value of the boolean equilibration variables to true or false to
     * either allow an equilibration step or skip it.
     *
     * @param equilibration a boolean.
     */
    public void setEquilibration(boolean equilibration) {
        this.equilibration = equilibration;
    }

    /**
     * Calls on the OST method set lambda to update lambda to the current value
     * in this class
     *
     * @param lambda a double.
     */
    public void setLambda(double lambda) {
        this.lambda = lambda;
        orthogonalSpaceTempering.setLambda(lambda);
    }

    /**
     * Returns the current value of lambda
     *
     * @return lambda
     */
    public double getLambda() {
        return lambda;
    }

    private void resetLevel() {
        this.currentLevel = regularLevel;
    }

    private void singleAction(MCActions act) {
        assert act != MCActions.COINFLIP;
        switch (act) {
            case LAMBDA: {
                singleStepLambda();
            }
            case HMC: {
                singleStepMD();
            }
            case METROPOLIZE: {
                metropolize();
            }
        }
    }

    /**
     * The goal is to sample lambda and coordinates (X) separately to converge
     * the ensemble average dU/dL for every state (lambda) along the
     * thermodynamic path.
     * <p>
     * 1.) At a fixed lambda, run a defined length MD trajectory to "move"
     * coordinates and dU/dL.
     * <p>
     * 2.) Accept / Reject the MD move using the total Hamiltonian (Kinetic
     * energy + OST energy).
     * <p>
     * 3.) Randomly change the value of Lambda.
     * <p>
     * 4.) Accept / Reject the Lambda move using the OST energy.
     * <p>
     * 5.) Add to the bias.
     */
    public void sample(long steps, MCAlgorithm algorithm) {
        // Validate the starting value of lambda.
        lambda = orthogonalSpaceTempering.getLambda();
        lambda = lambdaMove.validateLambda(lambda);
        orthogonalSpaceTempering.setLambda(lambda);

        int n = potential.getNumberOfVariables();
        double[] gradient = new double[n];
        double[] currentCoordinates = new double[n];
        long numMoves = steps / stepsPerMove;

        acceptXYZ = new RunningStatistics();
        acceptLambda = new RunningStatistics();
        acceptCombined = new RunningStatistics();

        // Initialize the current coordinates.
        potential.getCoordinates(currentCoordinates);

        // Update time dependent bias.
        Histogram histogram = orthogonalSpaceTempering.getHistogram();

        // Compute the current OST potential energy.
        double currentOSTEnergy = orthogonalSpaceTempering.energyAndGradient(currentCoordinates, gradient);

        // Collect the current dU/dL, Force Field Energy and Bias Energy.
        double currentdUdL = orthogonalSpaceTempering.getForceFielddEdL();
        initFFE = orthogonalSpaceTempering.getForceFieldEnergy();
        initBias = orthogonalSpaceTempering.getBiasEnergy();
        initKE = molDyn.getKineticEnergy();
        currentFFE = initFFE;
        currentBias = initBias;
        currentKE = initBias;

        MCActions[] acts = algorithm.actions;
        int nActs = acts.length;

        // Initialize MC move instances.
        for (int imove = 0; imove < numMoves; imove++) {
            long totalMoveTime = -nanoTime();

            String epS = equilibration ? equilString : prodString;
            String algS = equilibration ? equilAlgoString : algorithm.algoString;
            logger.info(String.format(algoFormatString, epS, algS));

            for (int i = 0; i < nActs; i++) {
                if (acts[i] == MCActions.COINFLIP) {
                    sinceLastMetropolis.add(MCActions.COINFLIP);
                    assert i < (nActs - 2);
                    boolean coinflip = random.nextBoolean();
                    int first;
                    int second;
                    if (coinflip) {
                        first = i + 1;
                        second = i + 2;
                    } else {
                        first = i + 2;
                        second = i + 1;
                    }
                    singleAction(acts[first]);
                    singleAction(acts[second]);

                    // Manually increment i.
                    i += 2;
                } else {
                    singleAction(acts[i]);
                }
            }

            resetLevel();

            // Run MD in an approximate potential U* (U star) that does not include the OST bias.
            /*long mdMoveTime = -nanoTime();
            mdMove.move();
            mdMoveTime += nanoTime();
            logger.log(verboseLoggingLevel, format("  Total time for MD move: %6.3f", mdMoveTime * NS2SEC));

            // Get the starting and final kinetic energy for the MD move.
            double currentKineticEnergy = mdMove.getInitialKinetic();
            double proposedKineticEnergy = mdMove.getKineticEnergy();

            // Get the new coordinates.
            potential.getCoordinates(proposedCoordinates);

            // Compute the Total OST Energy as the sum of the Force Field Energy and Bias Energy.
            long proposedOSTEnergyTime = -nanoTime();

            double proposedOSTEnergy;
            try {
                proposedOSTEnergy = orthogonalSpaceTempering.energyAndGradient(proposedCoordinates, gradient);
            } catch (EnergyException e) {
                mdMove.revertMove();
                logger.log(Level.INFO, " Unstable MD Move skipped.");
                continue;
            }
            proposedOSTEnergyTime += nanoTime();

            logger.fine(format("  Time to complete MD OST energy method call %6.3f", proposedOSTEnergyTime * NS2SEC));

            // Retrieve the proposed dU/dL, Force Field Energy and Bias Energy.
            double proposeddUdL = orthogonalSpaceTempering.getForceFielddEdL();
            double proposedForceFieldEnergy = orthogonalSpaceTempering.getForceFieldEnergy();
            double proposedBiasEnergy = orthogonalSpaceTempering.getBiasEnergy();

            // The Metropolis criteria is based on the sum of the OST Energy and Kinetic Energy.
            double currentTotalEnergy = currentOSTEnergy + currentKineticEnergy;
            double proposedTotalEnergy = proposedOSTEnergy + proposedKineticEnergy;

            logger.log(verboseLoggingLevel, format("\n  %8s %12s %12s %12s %12s", "", "Kinetic", "Potential", "Bias", "Total"));
            logger.log(verboseLoggingLevel, format("  Current  %12.4f %12.4f %12.4f %12.4f",
                    currentKineticEnergy, currentForceFieldEnergy, currentBiasEnergy, currentTotalEnergy));
            logger.log(verboseLoggingLevel, format("  Proposed %12.4f %12.4f %12.4f %12.4f",
                    proposedKineticEnergy, proposedForceFieldEnergy, proposedBiasEnergy, proposedTotalEnergy));
            logger.log(verboseLoggingLevel, format("  Delta    %12.4f %12.4f %12.4f %12.4f",
                    proposedKineticEnergy - currentKineticEnergy,
                    proposedForceFieldEnergy - currentForceFieldEnergy,
                    proposedBiasEnergy - currentBiasEnergy,
                    proposedTotalEnergy - currentTotalEnergy));

            double energyChange = mdMove.getEnergyChange();
            if (abs(energyChange) > ENERGY_CONSERVATION_TOLERANCE) {
                mdMove.revertMove();
                logger.warning(" MC Move skipped due to lack of MD energy conservation");
                continue;
            }

            if (orthogonalSpaceTempering.insideHardWallConstraint(orthogonalSpaceTempering.getLambda(), proposeddUdL) &&
                    evaluateMove(currentTotalEnergy, proposedTotalEnergy)) {
                // Accept MD move.
                acceptMD++;
                double percent = (acceptMD * 100.0) / (imove + 1);
                logger.info(format(" Accept [FL=%8.3f, E=%12.4f]  -> [FL=%8.3f,E=%12.4f] (%5.1f%%)",
                        currentdUdL, currentOSTEnergy, proposeddUdL, proposedOSTEnergy, percent));
                currentOSTEnergy = proposedOSTEnergy;
                currentdUdL = proposeddUdL;
                currentForceFieldEnergy = proposedForceFieldEnergy;
                currentBiasEnergy = proposedBiasEnergy;
                arraycopy(proposedCoordinates, 0, currentCoordinates, 0, n);
            } else {
                double percent = (acceptMD * 100.0) / (imove + 1);
                logger.info(format(" Reject [FL=%8.3f, E=%12.4f]  -> [FL=%8.3f,E=%12.4f] (%5.1f%%)",
                        currentdUdL, currentOSTEnergy, proposeddUdL, proposedOSTEnergy, percent));
                mdMove.revertMove();
            }
            mdMoveAndEvalTime += nanoTime();

            logger.fine(format("\n  Total time to run and evaluate MD move: %6.3f", mdMoveAndEvalTime * NS2SEC));

            // During equilibration, do not change Lambda or contribute to the OST bias.
            if (!equilibration) {
                // Update Lambda.
                logger.info(" MC Lambda Step");

                long lambdaMoveTime = -nanoTime();
                double currentLambda = orthogonalSpaceTempering.getLambda();
                lambdaMove.move();
                double proposedLambda = orthogonalSpaceTempering.getLambda();

                // Compute the Total OST Energy as the sum of the Force Field Energy and Bias Energy.
                long proposedOSTEnergyTime2 = -nanoTime();
                proposedOSTEnergy = orthogonalSpaceTempering.energyAndGradient(currentCoordinates, gradient);
                proposedOSTEnergyTime2 += nanoTime();

                logger.log(verboseLoggingLevel, format("  Time to complete Lambda OST energy method call %6.3f ", proposedOSTEnergyTime2 * NS2SEC));

                // Retrieve the proposed dU/dL, Force Field Energy and Bias Energy.
                proposedForceFieldEnergy = orthogonalSpaceTempering.getForceFieldEnergy();
                proposeddUdL = orthogonalSpaceTempering.getForceFielddEdL();

                logger.log(verboseLoggingLevel, format("\n  Current  OST     %12.3f at L=%5.3f.", currentOSTEnergy, currentLambda));
                logger.log(verboseLoggingLevel, format("  Proposed OST     %12.3f at L=%5.3f.", proposedOSTEnergy, proposedLambda));
                logger.log(verboseLoggingLevel, format("  MC Energy change: %12.3f (kcal/mol).", proposedOSTEnergy - currentOSTEnergy));

                if (orthogonalSpaceTempering.insideHardWallConstraint(proposedLambda, proposeddUdL) &&
                        evaluateMove(currentOSTEnergy, proposedOSTEnergy)) {
                    acceptLambda++;
                    double percent = (acceptLambda * 100.0) / (imove + 1);
                    logger.info(format("  Accept [ L=%8.3f,E=%12.4f]   -> [ L=%8.3f,E=%12.4f] (%5.1f%%)",
                            currentLambda, currentOSTEnergy, proposedLambda, proposedOSTEnergy, percent));
                    currentForceFieldEnergy = proposedForceFieldEnergy;
                    currentdUdL = proposeddUdL;
                    lambda = proposedLambda;
                } else {
                    double percent = (acceptLambda * 100.0) / (imove + 1);
                    logger.info(format("  Reject [ L=%8.3f,E=%12.4f]   -> [ L=%8.3f,E=%12.4f] (%5.1f%%)",
                            currentLambda, currentOSTEnergy, proposedLambda, proposedOSTEnergy, percent));
                    lambdaMove.revertMove();
                    lambda = currentLambda;
                }

                lambdaMoveTime += nanoTime();
                logger.log(verboseLoggingLevel, format("  Lambda move completed in %6.3f", lambdaMoveTime * NS2SEC));

                if (imove % biasDepositionFrequency == 0) {
                    histogram.addBias(currentdUdL, currentCoordinates, null);
                } else {
                    // TODO: Step down to FINE when we know this works.
                    logger.log(Level.INFO, format(" Cycle %d: skipping bias deposition.", imove));
                }

                logger.log(verboseLoggingLevel, format("  Added Bias at [L=%5.3f, FL=%9.3f]", lambda, currentdUdL));

                // Compute the updated OST bias.
                currentBiasEnergy = histogram.computeBiasEnergy(lambda, currentdUdL);

                // Update the current OST Energy to be the sum of the current Force Field Energy and updated OST Bias.
                currentOSTEnergy = currentForceFieldEnergy + currentBiasEnergy;

                // TODO: Call listeners.
                /*if (automaticWriteouts) {
                    long mdMoveNum = imove * stepsPerMove;
                    boolean trySnapshot = lambda >= orthogonalSpaceTempering.getLambdaWriteOut();
                    EnumSet<MolecularDynamics.WriteActions> written = mdMove.writeFilesForStep(mdMoveNum, trySnapshot, true);
                    if (written.contains(MolecularDynamics.WriteActions.RESTART)) {
                        orthogonalSpaceTempering.writeAdditionalRestartInfo(false);
                    }
                }*/
            /*}

            totalMoveTime += nanoTime();
            logger.info(format(" Round complete in %6.3f sec.", totalMoveTime * NS2SEC));*/
        }
    }

    /**
     * Propose a lambda move.
     *
     * @return The proposed lambda.
     */
    private double singleStepLambda() {
        lambdaMove.move();
        double proposedLambda = orthogonalSpaceTempering.getLambda();
        logger.log(verboseLoggingLevel, format(" Proposed lambda: %5.3f.", proposedLambda));
        sinceLastMetropolis.add(MCActions.LAMBDA);
        reverts.push(mdMove::revertMove);
        return proposedLambda;
    }

    /**
     * Run MD in an approximate potential U* (U star) that does not include the OST bias.
     */
    private void singleStepMD() {
        long mdMoveTime = -nanoTime();
        mdMove.move();
        mdMoveTime += nanoTime();
        sinceLastMetropolis.add(MCActions.HMC);
        reverts.push(mdMove::revertMove);
        logger.log(verboseLoggingLevel, format(" Total time for MD move: %6.3f", mdMoveTime * NS2SEC));
    }

    private void resetRevertStack(boolean accepted) {
        if (accepted) {
            while (reverts.size() > 0) {
                reverts.pop().run();
            }
        } else {
            while (reverts.size() > 0) {
                reverts.pop();
            }
        }
    }

    private boolean metropolize() {
        potential.getCoordinates(proposedCoordinates);
        double proposedPE = potential.energy(proposedCoordinates, false);
        double proposedKE = mdMove.getKineticEnergy();
        double proposedH = proposedPE + proposedKE;
        double currentH = currentFFE + currentKE + currentBias;

        // Retrieve the proposed dU/dL, Force Field Energy and Bias Energy.
        double proposeddUdL = orthogonalSpaceTempering.getForceFielddEdL();
        double proposedLambda = orthogonalSpaceTempering.getLambda();
        double proposedFFE = orthogonalSpaceTempering.getForceFieldEnergy();
        double proposedBias = orthogonalSpaceTempering.getBiasEnergy();

        double deltaKE = proposedKE - currentKE;
        double deltaFFE = proposedFFE - currentFFE;
        double deltaBias = proposedBias - currentBias;
        double deltaH = proposedH - currentH;

        double deltaHStar = deltaKE + deltaFFE;
        
        boolean accepted;
        boolean hardWallReject = false;
        if (!orthogonalSpaceTempering.insideHardWallConstraint(proposedLambda, proposeddUdL)) {
            accepted = false;
            hardWallReject = true;
        } else if (Math.abs(deltaHStar) > ENERGY_CONSERVATION_TOLERANCE) {
            accepted = false;
        } else if (!evaluateMove(currentH, proposedH)) {
            accepted = false;
        } else {
            accepted = true;
        }

        double acceptVal = accepted ? 1.0 : 0.0;

        if (sinceLastMetropolis.contains(MCActions.COINFLIP)) {
            acceptCombined.addValue(acceptVal);
        } else if (sinceLastMetropolis.size() == 1) {
            if (sinceLastMetropolis.get(0) == MCActions.HMC) {
                acceptXYZ.addValue(acceptVal);
            } else if (sinceLastMetropolis.get(0) == MCActions.LAMBDA) {
                acceptLambda.addValue(acceptVal);
            }
        }

        boolean doLog = logger.isLoggable(currentLevel);

        if (doLog) {
            StringBuilder sb = new StringBuilder("\n   Kinetic    Potential         Bias        Total\n");
            sb.append(String.format("  Current  %12.4f %12.4f %12.4f %12.4f\n", currentKE, currentFFE, currentBias, currentH));
            sb.append(String.format("  Proposed %12.4f %12.4f %12.4f %12.4f\n", proposedKE, proposedFFE, proposedBias, proposedH));
            sb.append(String.format("  Delta    %12.4f %12.4f %12.4f %12.4f", deltaKE, deltaFFE, deltaBias, deltaH));

            logger.log(currentLevel, sb.toString());
            //resetLevel();
        }

        resetRevertStack(accepted);
        
        /*

            logger.log(verboseLoggingLevel, format("\n  %8s %12s %12s %12s %12s", "", "Kinetic", "Potential", "Bias", "Total"));
            logger.log(verboseLoggingLevel, format("  Current  %12.4f %12.4f %12.4f %12.4f",
                    currentKineticEnergy, currentForceFieldEnergy, currentBiasEnergy, currentTotalEnergy));
            logger.log(verboseLoggingLevel, format("  Proposed %12.4f %12.4f %12.4f %12.4f",
                    proposedKineticEnergy, proposedForceFieldEnergy, proposedBiasEnergy, proposedTotalEnergy));
            logger.log(verboseLoggingLevel, format("  Delta    %12.4f %12.4f %12.4f %12.4f",
                    proposedKineticEnergy - currentKineticEnergy,
                    proposedForceFieldEnergy - currentForceFieldEnergy,
                    proposedBiasEnergy - currentBiasEnergy,
                    proposedTotalEnergy - currentTotalEnergy));

            double energyChange = mdMove.getEnergyChange();
            if (abs(energyChange) > ENERGY_CONSERVATION_TOLERANCE) {
                mdMove.revertMove();
                logger.warning(" MC Move skipped due to lack of MD energy conservation");
                continue;
            }

            if (orthogonalSpaceTempering.insideHardWallConstraint(orthogonalSpaceTempering.getLambda(), proposeddUdL) &&
                    evaluateMove(currentTotalEnergy, proposedTotalEnergy)) {
                // Accept MD move.
                acceptMD++;
                double percent = (acceptMD * 100.0) / (imove + 1);
                logger.info(format(" Accept [FL=%8.3f, E=%12.4f]  -> [FL=%8.3f,E=%12.4f] (%5.1f%%)",
                        currentdUdL, currentOSTEnergy, proposeddUdL, proposedOSTEnergy, percent));
                currentOSTEnergy = proposedOSTEnergy;
                currentdUdL = proposeddUdL;
                currentForceFieldEnergy = proposedForceFieldEnergy;
                currentBiasEnergy = proposedBiasEnergy;
                arraycopy(proposedCoordinates, 0, currentCoordinates, 0, n);
            } else {
                double percent = (acceptMD * 100.0) / (imove + 1);
                logger.info(format(" Reject [FL=%8.3f, E=%12.4f]  -> [FL=%8.3f,E=%12.4f] (%5.1f%%)",
                        currentdUdL, currentOSTEnergy, proposeddUdL, proposedOSTEnergy, percent));
                mdMove.revertMove();
            }
            mdMoveAndEvalTime += nanoTime();

            logger.fine(format("\n  Total time to run and evaluate MD move: %6.3f", mdMoveAndEvalTime * NS2SEC));
         */
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected double currentEnergy() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void storeState() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void revertStep() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getX() {
        return molDyn.getX();
    }

    @Override
    public double[] copyX(@Nullable double[] x) {
        return molDyn.copyX(x);
    }

    @Override
    public double[] getV() {
        return molDyn.getV();
    }

    @Override
    public double[] copyV(@Nullable double[] v) {
        return molDyn.copyV(v);
    }

    @Override
    public double[] getA() {
        return molDyn.getA();
    }

    @Override
    public double[] copyA(@Nullable double[] a) {
        return molDyn.copyA(a);
    }

    @Override
    public double[] getPriorAcceleration() {
        return molDyn.getPriorAcceleration();
    }

    @Override
    public double[] copyPriorAcceleration(@Nullable double[] aPrev) {
        return molDyn.copyPriorAcceleration(aPrev);
    }

    @Override
    public boolean topLevelAlgorithm() {
        return false;
    }

    @Override
    public void distributeListeners(long steps, long exteriorSteps) {

    }

    @Override
    public Optional<Crystal> getCrystal() {
        Crystal c = potential instanceof CrystalPotential ? ((CrystalPotential) potential).getCrystal() : null;
        return Optional.ofNullable(c);
    }

    @Override
    public boolean actionForStep(long step) {
        return false;
    }

    private enum MCActions {
        HMC, LAMBDA, COINFLIP, METROPOLIZE;
    }

    public enum MCAlgorithm {
        ONE_STEP("1-step coinflip ", MCActions.COINFLIP, MCActions.HMC, MCActions.LAMBDA, MCActions.METROPOLIZE),
        TWO_STEP("2-step staggered ", MCActions.LAMBDA, MCActions.HMC, MCActions.METROPOLIZE);

        private final String algoString;
        private final MCActions[] actions;

        MCAlgorithm(String alS, MCActions... actions) {
            this.algoString = alS;
            this.actions = actions;
        }
    }
}
