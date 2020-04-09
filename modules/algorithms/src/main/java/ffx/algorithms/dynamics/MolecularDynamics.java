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

import java.io.File;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.function.DoubleConsumer;
import java.util.function.DoubleSupplier;
import java.util.function.Supplier;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Stream;

import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;

import ffx.algorithms.listeners.CompoundListener;
import ffx.algorithms.listeners.DynListener;
import ffx.algorithms.listeners.IntervalListener;
import ffx.algorithms.listeners.TrajectoryListener;
import ffx.potential.parsers.SystemFilter;
import ffx.utilities.MDListener;
import ffx.crystal.CrystalPotential;
import ffx.utilities.Constants;
import org.apache.commons.collections4.queue.CircularFifoQueue;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import edu.rit.pj.Comm;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import ffx.algorithms.dynamics.integrators.Integrator;
import ffx.algorithms.dynamics.integrators.Respa;
import ffx.algorithms.dynamics.integrators.Stochastic;
import ffx.algorithms.dynamics.thermostats.Adiabatic;
import ffx.algorithms.dynamics.thermostats.Thermostat;
import ffx.algorithms.mc.MonteCarloListener;
import ffx.crystal.Crystal;
import ffx.numerics.Constraint;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parsers.DYNFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.potential.utils.EnergyException;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;

import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.NS2SEC;

/**
 * Run NVE, NVT, or NPT molecular dynamics.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class MolecularDynamics extends AbstractDynAlg implements Runnable, Terminatable  {

    private static final Logger logger = Logger.getLogger(MolecularDynamics.class.getName());
    /**
     * MolecularAssembly to run dynamics on.
     */
    protected final MolecularAssembly molecularAssembly;
    /**
     * List of MolecularAssembly instances.
     */
    protected final AssemblyInfo[] assemblies;
    /**
     * Propogate dynamics on this potential surface.
     */
    protected final Potential potential;
    /**
     * Monte Carlo listener.
     */
    private MonteCarloListener monteCarloListener;
    /**
     * An Algorithm Listener to send updates to the GUI.
     */
    protected final AlgorithmListener algorithmListener;
    /**
     * Thermostat instance.
     */
    protected final Thermostat thermostat;
    /**
     * Integrator instance.
     */
    protected final Integrator integrator;
    /**
     * Flag to indicate use of constant pressure.
     */
    protected final boolean constantPressure;
    /**
     * Any Barostat that may be in use.
     */
    protected final Barostat barostat;
    /**
     * Flag to indicate MD should be terminated.
     */
    private boolean terminate = false;
    /**
     * Number of MD steps to take.
     */
    private long nSteps = 1000;
    /**
     * State of the dynamics.
     */
    private DynamicsState dynamicsState;
    /**
     * Total simulation time.
     */
    protected double totalSimTime = 0.0;
    /**
     * Indicates how verbose MD should be.
     */
    private final VerbosityLevel verbosityLevel;
    /**
     * Log some of the more frequent messages at this level. Always at or below basicLogging.
     */
    private final Level intermediateLogging;
    /**
     * Log basic information at this level. Always at or above intermediateLogging.
     */
    protected final Level basicLogging;
    /**
     * Dynamics restart file.
     */
    private final File restartFile;
    /**
     * Flag to indicate loading of restart file.
     */
    private boolean loadRestart = false;
    /**
     * Filter to parse the dynamics restart file.
     */
    private final DYNFilter dynFilter;
    /**
     * Flag to indicate dynamics has been initialized.
     */
    private boolean initialized;
    /**
     * Flag to indicate a run has finished.
     */
    protected boolean done;
    /**
     * Flag to indicate velocities should be initialized.
     */
    protected boolean initVelocities = true;
    /**
     * Number of dynamics variables.
     */
    private final int numberOfVariables;
    /**
     * Coordinates.
     */
    protected final double[] x;
    /**
     * Velocities.
     */
    protected final double[] v;
    /**
     * Accelerations.
     */
    protected final double[] a;
    /**
     * Previous accelerations.
     */
    protected final double[] aPrevious;
    /**
     * The gradient.
     */
    protected final double[] gradient;
    /**
     * Mass for each degree of freedom.
     */
    protected final double[] mass;
    /**
     * Time step (picoseconds).
     */
    protected final double dt;
    /**
     * TIme step in fs
     */
    private final double timeStepFS;
    /**
     * Whether MD handles writing restart/trajectory files itself (true), or will be
     * commanded by another class (false) to do it. The latter is true for MC-OST, for example.
     */
    protected final boolean automaticWriteouts;
    public static final double DEFAULT_RESTART_INTERVAL = 1.0;
    public static final double DEFAULT_LOG_INTERVAL = 0.25;
    public static final double DEFAULT_TRAJECTORY_INTERVAL = 10.0;
    private final List<MDListener> allListeners = new ArrayList<>();
    protected final List<MDListener> interiorListeners = new ArrayList<>();
    private final List<MDListener> exteriorListeners = new ArrayList<>();
    /**
     * Time between writing out restart/checkpoint files in picoseconds.
     */
    protected final double restartInterval;
    /**
     * Timesteps between writing out restart/checkpoint files. Set by the init method.
     */
    protected final int restartFrequency;
    /**
     * Time between appending to the trajectory file in picoseconds.
     */
    private final double snapshotInterval;
    /**
     * Timesteps between adding a frame to the trajectory file. Set by the init method.
     */
    protected final int snapshotFrequency;
    /**
     * Time between logging information to the screen in picoseconds.
     */
    private final double logInterval;
    /**
     * TImesteps between logging information to the screen. Set by the init method.
     */
    protected final int logFrequency;
    /**
     * Target temperature. ToDo: use the Thermostat instance.
     */
    protected double targetTemperature;
    /**
     * Current temperature.
     */
    double currentTemperature;
    /**
     * Temperature of the system as of the last init call (and start of the last dynamics run).
     */
    protected double initialTemp;
    /**
     * Current kinetic energy.
     */
    double currentKineticEnergy;
    /**
     * Kinetic energy of the system as of the last init call (and start of the last dynamics run).
     */
    protected double initialKinetic;
    /**
     * Current potential energy.
     */
    double currentPotentialEnergy;
    /**
     * Potential energy of the system as of the last init call (and start of the last dynamics run).
     */
    protected double initialPotential;
    /**
     * Current total energy.
     */
    double currentTotalEnergy;
    /**
     * Total energy of the system as of the last init call (and start of the last dynamics run).
     */
    protected double initialTotal;
    /**
     * Save snapshots in PDB format.
     */
    private final boolean saveSnapshotAsPDB;
    /**
     * Keep some old coordinate snapshots around.
     */
    private final int numSnapshotsToKeep;
    /**
     * Circular FIFO queues will simply discard old elements.
     */
    private final CircularFifoQueue<CoordinateSnapshot> lastSnapshots;
    /**
     * MC notification flag.
     */
    private MonteCarloNotification mcNotification = MonteCarloNotification.NEVER;
    /**
     * By default, wait 25000 ns (25 us) in between polling the dynamics thread.
     */
    public static final int DEFAULT_DYNAMICS_SLEEP_TIME = 25000;
    /**
     * Wait this many nanoseconds in between polling the dynamics thread.
     */
    private final int dynSleepTime;
    /**
     * List of (extra) suppliers of extra lines to add to snapshot headers.
     */
    private final List<Supplier<Stream<String>>> snapshotLineProviders = new ArrayList<>();
    /**
     * Logs MD-related information.
     */
    protected final ThermoLogListener tll;
    /**
     * The layer above MD, or the MD itself if it's the top layer.
     */
    private DynamicAlgorithm aboveLayer;

    /**
     * <p>
     * dynamicsFactory.</p>
     *
     * @param opts     Container object for various options.
     * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
     * @param engine   a {@link MolecularDynamics.DynamicsEngine} object.
     * @return a {@link MolecularDynamics} object.
     */
    public static MolecularDynamics dynamicsFactory(MolecularDynamicsOptions opts, AlgorithmListener listener,
                                                    DynamicsEngine engine) {
        engine = engine == null ? defaultEngine(opts) : engine;
        switch (engine) {
            case OPENMM:
                Potential potentialEnergy = opts.potential;
                // TODO: Replace this with calls to the leaves of a proper tree structure.
                // Unfortunately, neither Java, nor Apache Commons, nor Guava has an arbitrary tree implementing Collection.
                // Nor does javax.swing have a quick "get me the leaves" method that I was able to find.
                boolean ommLeaves = potentialEnergy.getUnderlyingPotentials().stream().
                        anyMatch((Potential p) -> p instanceof ForceFieldEnergyOpenMM);
                ommLeaves = ommLeaves || potentialEnergy instanceof ForceFieldEnergyOpenMM;
                if (ommLeaves) {
                    // TODO: Fix MD-OMM.
                    /*return new MolecularDynamicsOpenMM(assembly,
                            potentialEnergy, properties, listener, requestedThermostat, requestedIntegrator);*/
                    return null;
                } else {
                    throw new IllegalArgumentException(format(" Requested OpenMM engine %s, but at least one leaf of the potential %s is not an OpenMM force field!", engine, potentialEnergy));
                }
            case FFX:
            default:
                return new MolecularDynamics(opts, listener);
        }
    }

    private static DynamicsEngine defaultEngine(MolecularDynamicsOptions opts) {
        MolecularAssembly molecularAssembly = opts.assemblies[0];
        Potential potentialEnergy = opts.potential;
        CompositeConfiguration properties = molecularAssembly.getProperties();
        String mdEngine = properties.getString("MD-engine");
        if (mdEngine != null) {
            if (mdEngine.equalsIgnoreCase("OMM")) {
                logger.info(" Creating OpenMM Dynamics Object");
                return DynamicsEngine.OPENMM;
            } else {
                logger.info(" Creating FFX Dynamics Object");
                return DynamicsEngine.FFX;
            }
        } else {
            // TODO: Replace this with a better check.
            boolean ommLeaves = potentialEnergy.getUnderlyingPotentials().stream().
                    anyMatch((Potential p) -> p instanceof ForceFieldEnergyOpenMM);
            ommLeaves = ommLeaves || potentialEnergy instanceof ForceFieldEnergyOpenMM;
            if (ommLeaves) {
                return DynamicsEngine.OPENMM;
            } else {
                return DynamicsEngine.FFX;
            }
        }
    }

    /**
     * <p>
     * Constructor for MolecularDynamics.</p>
     *
     * @param opts     Container object for various options.
     * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
     */
    protected MolecularDynamics(MolecularDynamicsOptions opts, AlgorithmListener listener) {
        this.assemblies = Arrays.stream(opts.assemblies).
                map(AssemblyInfo::new).
                toArray(AssemblyInfo[]::new);
        this.molecularAssembly = opts.assemblies[0];
        this.algorithmListener = listener;
        this.potential = opts.potential;
        this.barostat = opts.barostat;
        constantPressure = barostat != null;

        dynSleepTime = opts.dynSleepTime;
        
        numberOfVariables = potential.getNumberOfVariables();
        mass = opts.mass;
        x = opts.x;
        v = opts.v;
        a = opts.a;
        aPrevious = opts.aPrevious;
        gradient = new double[numberOfVariables];

        integrator = opts.integrator;
        thermostat = opts.thermostat;

        targetTemperature = thermostat.getTargetTemperature();

        numSnapshotsToKeep = opts.numSnapshotsToKeep;
        // Cannot construct a CircularFifoQueue of zero length.
        lastSnapshots = new CircularFifoQueue<>(Math.max(numSnapshotsToKeep, 1));

        verbosityLevel = opts.vLevel;
        basicLogging = verbosityLevel.basicLevel;
        intermediateLogging = verbosityLevel.intmdLevel;
        this.restartFile = opts.dynFile;
        this.dynFilter = new DYNFilter(molecularAssembly.getName());
        dt = opts.dt;
        timeStepFS = dt * Constants.PSEC_TO_FSEC;

        restartInterval = opts.restartInterval;
        restartFrequency = (int) Math.max(1.0, restartInterval / dt);
        logInterval = opts.logInterval;
        logFrequency = (int) Math.max(1.0, logInterval / dt);
        snapshotInterval = opts.snapshotInterval;
        snapshotFrequency = (int) Math.max(1.0, snapshotInterval / dt);

        saveSnapshotAsPDB = opts.saveSnapshotAsPDB;
        automaticWriteouts = opts.automaticWriteouts;

        if (potential instanceof LambdaInterface) {
            snapshotLineProviders.add(() -> {
                String lamString = String.format("Lambda: %.8f", ((LambdaInterface) potential).getLambda());
                return Stream.of(lamString);
            });
        }

        Comm world = Comm.world();
        if (world != null && world.size() > 1) {
            snapshotLineProviders.add(() -> {
                String rankString = String.format("Rank: %d", world.rank());
                return Stream.of(rankString);
            });
        }

        generateMDListeners();
        tll = new ThermoLogListener(logFrequency, 0);
        allListeners.add(tll);
        aboveLayer = this;

        done = true;
    }

    /**
     * Gives MD a reference to the algorithm it "belongs" to.
     * 
     * @param layer Algorithm directly above MD.
     */
    public void setAboveLayer(DynamicAlgorithm layer) {
        if (initialized) {
            logger.warning(" Attempting to set a layer above an initialized MD!");
        }
        aboveLayer = layer;
    }

    public void addLineProvider(Supplier<Stream<String>> lineProvider) {
        snapshotLineProviders.add(lineProvider);
    }

    protected List<Runnable> restartActions() {
        return Collections.emptyList();
    }

    protected List<Runnable> snapshotActions() {
        return Collections.emptyList();
    }

    private String[] getExtraLines() {
        return snapshotLineProviders.stream().
                flatMap(Supplier::get).
                toArray(String[]::new);
    }

    private void generateMDListeners() {
        DynListener dynL = new DynListener(restartFile, dynFilter, restartActions(),
                restartFrequency, 0, this, basicLogging);
        allListeners.add(dynL);

        List<MDListener> snapshotListeners = new ArrayList<>(assemblies.length);
        for (AssemblyInfo ai : assemblies) {
            File outFi = saveSnapshotAsPDB ? ai.pdbFile : ai.archiveFile;
            SystemFilter filter = saveSnapshotAsPDB ? ai.pdbFilter : ai.xyzFilter;
            TrajectoryListener trjL = new TrajectoryListener(outFi, filter, this::getExtraLines,
                    snapshotFrequency, 0, Collections.emptyList(), basicLogging);
            // The trajectory listener's file will be updated through the AssemblyInfo.
            ai.setTrajectoryListener(trjL);
            snapshotListeners.add(trjL);
        }
        CompoundListener snapListener = new CompoundListener(snapshotFrequency, 0, snapshotActions(), snapshotListeners);
        allListeners.add(snapListener);

        allListeners.addAll(potential.getAdditionalListeners(snapshotFrequency, restartFrequency, logFrequency, 0));
    }

    // TODO: Consider either re-implementing reinit method or just building a new MD object each time.

    /**
     * Not meaningful for FFX dynamics (no need to obtain velocities/accelerations from
     * a different program, especially one running on a GPU). Is a no-op.
     *
     * @param obtainVA Not meaningful for this implementation.
     */
    public void setObtainVelAcc(boolean obtainVA) {
        // Not meaningful for FFX dynamics.
    }

    /**
     * <p>
     * Getter for the field <code>thermostat</code>.</p>
     *
     * @return a {@link ffx.algorithms.dynamics.thermostats.Thermostat} object.
     */
    public Thermostat getThermostat() {
        return thermostat;
    }

    /**
     * <p>
     * Setter for the field <code>monteCarloListener</code>.</p>
     *
     * @param listener a {@link MonteCarloListener} object.
     * @param when     a {@link MolecularDynamics.MonteCarloNotification} object.
     */
    public void setMonteCarloListener(MonteCarloListener listener, MonteCarloNotification when) {
        if (initialized) {
            logger.warning(" Attempting to set the Monte Carlo listener for an already-initialized MD!");
        }
        monteCarloListener = listener;
        mcNotification = when;
    }

    /**
     * <p>
     * init</p>
     *
     * @param nSteps             Number of MD steps
     * @param temperature        Temperature in Kelvins.
     * @param initVelocities     Initialize new velocities from a Maxwell-Boltzmann distribution.
     */
    public void init(final long nSteps, final double temperature, final boolean initVelocities) {
        // Return if already running.
        if (!done) {
            logger.warning(" Programming error - attempt to modify parameters of a running MolecularDynamics instance.");
            return;
        }

        if (topLevelAlgorithm()) {
            distributeListeners(nSteps, 0);
        }

        if (integrator instanceof Stochastic) {
            if (constantPressure) {
                logger.log(basicLogging, "\n Stochastic dynamics in the NPT ensemble");
            } else {
                logger.log(basicLogging, "\n Stochastic dynamics in the NVT ensemble");
            }
        } else if (!(thermostat instanceof Adiabatic)) {
            if (constantPressure) {
                logger.log(basicLogging, "\n Molecular dynamics in the NPT ensemble");
            } else {
                logger.log(basicLogging, "\n Molecular dynamics in the NVT ensemble");
            }
        } else {
            if (constantPressure) {
                logger.severe("\n NPT Molecular dynamics requires a thermostat");
            } else {
                logger.log(basicLogging, "\n Molecular dynamics in the NVE ensemble");
            }
        }

        this.nSteps = nSteps;
        totalSimTime = 0.0;

        loadRestart = restartFile.exists() && !initialized;

        this.targetTemperature = temperature;
        this.initVelocities = initVelocities;
        done = false;

        if (loadRestart) {
            logger.info("  Continuing from " + restartFile.getAbsolutePath());
        }

        if (!verbosityLevel.isQuiet) {
            logger.info(format("  Number of steps:     %8d", nSteps));
            logger.info(format("  Time step:           %8.3f (fsec)", timeStepFS));
            logger.info(format("  Print interval:      %8.3f (psec)", logInterval));
            logger.info(format("  Save interval:       %8.3f (psec)", snapshotInterval));
            if (assemblies.length > 1) {
                for (int i = 0; i < assemblies.length; i++) {
                    AssemblyInfo ai = assemblies[i];
                    logger.info(format("  Archive file %3d: %s", (i + 1), ai.archiveFile.getName()));
                }
            } else {
                logger.info(format("  Archive file:     %s", assemblies[0].archiveFile.getName()));
            }
            logger.info(format("  Restart file:     %s", restartFile.getName()));
        }
    }

    /**
     * Blocking molecular dynamics. When this method returns, the MD run is
     * done.
     *
     * @param nSteps             Number of MD steps
     * @param temperature        Temperature in Kelvins.
     * @param initVelocities     Initialize new velocities from a Maxwell-Boltzmann distribution.
     */
    public void dynamic(final long nSteps, final double temperature, final boolean initVelocities) {
        // Return if already running;
        // Could happen if two threads call dynamic on the same MolecularDynamics instance.
        if (!done) {
            logger.warning(" Programming error - a thread invoked dynamic when it was already running.");
            return;
        }

        init(nSteps, temperature, initVelocities);

        Thread dynamicThread = new Thread(this);
        dynamicThread.start();
        synchronized (this) {
            try {
                while (dynamicThread.isAlive()) {
                    wait(0, dynSleepTime);
                }
            } catch (InterruptedException e) {
                String message = " Molecular dynamics interrupted.";
                logger.log(Level.WARNING, message, e);
            }
        }
        if (!verbosityLevel.isQuiet) {
            logger.info(" Done with an MD round.");
        }
    }

    /**
     * Performs basic pre-MD operations such as loading the restart file.
     */
    void preRunOps() throws IllegalStateException {
        done = false;
        terminate = false;

        // Set the target temperature.
        thermostat.setTargetTemperature(targetTemperature);
        boolean quiet = verbosityLevel.isQuiet();
        thermostat.setQuiet(quiet);
        if (integrator instanceof Stochastic) {
            Stochastic stochastic = (Stochastic) integrator;
            stochastic.setTemperature(targetTemperature);
        }

        // Set the step size.
        integrator.setTimeStep(dt);

        if (!initialized) {
            // Initialize from a restart file.
            if (loadRestart) {
                Crystal crystal = molecularAssembly.getCrystal();
                if (!dynFilter.readDYN(restartFile, crystal, x, v, a, aPrevious)) {
                    String message = " Could not load the restart file - dynamics terminated.";
                    logger.log(Level.WARNING, message);
                    done = true;
                    throw new IllegalStateException(message);
                } else {
                    molecularAssembly.setCrystal(crystal);
                }
            } else {
                // Initialize using current atomic coordinates.
                potential.getCoordinates(x);
                // Initialize atomic velocities from a Maxwell-Boltzmann distribution or set to 0.
                if (initVelocities) {
                    thermostat.maxwell(targetTemperature);
                } else {
                    fill(v, 0.0);
                }
            }
        } else {
            // If MD has already been run (ie. Annealing or RepEx), then initialize velocities if requested.
            if (initVelocities) {
                thermostat.maxwell(targetTemperature);
            }
        }
    }

    /**
     * Initializes energy fields, esp. potential energy.
     */
    private void initializeEnergies() {
        // Compute the current potential energy.
        try {
            currentPotentialEnergy = potential.energyAndGradient(x, gradient);
        } catch (EnergyException ex) {
            writeStoredSnapshots();
            throw ex;
        }
        initialPotential = currentPotentialEnergy;

        // Initialize current and previous accelerations.
        if (!loadRestart || initialized || integrator instanceof Respa) {
            // For the Respa integrator, initial accelerations are from the slowly varying forces.
            if (integrator instanceof Respa) {
                potential.setEnergyTermState(Potential.STATE.SLOW);
                potential.energyAndGradient(x, gradient);
            }

            for (int i = 0; i < numberOfVariables; i++) {
                a[i] = -KCAL_TO_GRAM_ANG2_PER_PS2 * gradient[i] / mass[i];
            }

            if (aPrevious != null) {
                arraycopy(a, 0, aPrevious, 0, numberOfVariables);
            }
        }

        // Compute the current kinetic energy.
        thermostat.computeKineticEnergy();
        currentKineticEnergy = thermostat.getKineticEnergy();
        initialKinetic = currentKineticEnergy;
        currentTemperature = thermostat.getCurrentTemperature();
        initialTemp = currentTemperature;
        currentTotalEnergy = currentKineticEnergy + currentPotentialEnergy;
        initialTotal = currentTotalEnergy;
    }

    /**
     * Pre-run operations (mostly logging) that require knowledge of system energy.
     */
    void postInitEnergies() {
        initialized = true;

        logger.log(basicLogging, format("\n  %8s %12s %12s %12s %8s %8s", "Time", "Kinetic", "Potential", "Total", "Temp", "CPU"));
        logger.log(basicLogging, format("  %8s %12s %12s %12s %8s %8s", "psec", "kcal/mol", "kcal/mol", "kcal/mol", "K", "sec"));
        logger.log(basicLogging, format("  %8s %12.4f %12.4f %12.4f %8.2f",
                "", currentKineticEnergy, currentPotentialEnergy, currentTotalEnergy, currentTemperature));

        // Store the initialized state.
        storeState();
    }

    /**
     * Post-run cleanup operations.
     */
    void postRun() {
        // Add the potential energy of the slow degrees of freedom.
        if (integrator instanceof Respa) {
            potential.setEnergyTermState(Potential.STATE.BOTH);
        }

        // Log normal completion.
        if (!terminate) {
            logger.log(basicLogging, format(" Completed %8d time steps\n", nSteps));
        }

        // Reset the done and terminate flags.
        done = true;
        terminate = false;

        if (monteCarloListener != null && mcNotification == MonteCarloNotification.AFTER_DYNAMICS) {
            monteCarloListener.mcUpdate(thermostat.getCurrentTemperature());
        }
    }

    protected class ThermoLogListener extends IntervalListener {
        private long time;

        public ThermoLogListener(long frequency, long priorSteps) {
            super(frequency, priorSteps);
            time = -System.nanoTime();
        }

        protected void setTime() {
            time = -System.nanoTime();
        }

        @Override
        public boolean actionForStep(long step) {
            if (checkWrite(step)) {
                time += System.nanoTime();
                logger.log(intermediateLogging, format(" %7.3e %12.4f %12.4f %12.4f %8.2f %8.3f",
                        totalSimTime, currentKineticEnergy, currentPotentialEnergy,
                        currentTotalEnergy, currentTemperature, time * NS2SEC));
                for (AssemblyInfo assembly : assemblies) {
                    // Probably unwise to parallelize this, so that it doesn't
                    // hit the GUI with parallel updates.
                    algorithmListener.algorithmUpdate(assembly.getAssembly());
                }
                setTime();
                return true;
            }
            return false;
        }
    }

    /**
     * Main loop of the run method.
     */
    private void mainLoop() {
        tll.setTime();
        double defaultDeltaPEThresh = 1.0E6;
        boolean removeCOMM = thermostat.getRemoveCenterOfMassMotion();

        // Instead of having if-else blocks inside the main loop, use functional programming!
        DoubleConsumer mcUpdate;
        if (monteCarloListener != null && mcNotification == MonteCarloNotification.EACH_STEP) {
            mcUpdate = monteCarloListener::mcUpdate;
        } else {
            mcUpdate = (double val) -> {};
        }
        DoubleSupplier respaE;
        if (integrator instanceof Respa) {
            respaE = ((Respa) integrator)::getHalfStepEnergy;
        } else {
            respaE = () -> 0.0;
        }

        // Integrate Newton's equations of motion for the requested number of steps,
        // unless early termination is requested.
        for (long step = 1; step <= nSteps; step++) {
            if (step > 1) {
                List<Constraint> constraints = potential.getConstraints();
                // TODO: Replace magic numbers with named constants.
                long constraintFails = constraints.stream().
                        filter((Constraint c) -> !c.constraintSatisfied(x, v, 1E-7, 1E-7)).
                        count();
                if (constraintFails > 0) {
                    logger.info(format(" %d constraint failures in step %d", constraintFails, step));
                }
            }
            /* Notify MonteCarlo handlers such as PhMD or rotamer drivers. */
            mcUpdate.accept(thermostat.getCurrentTemperature());

            // Do the half-step thermostat operation.
            thermostat.halfStep(dt);

            // Do the half-step integration operation.
            integrator.preForce(potential);

            // Compute the potential energy and gradients.
            double priorPE = currentPotentialEnergy;
            try {
                currentPotentialEnergy = potential.energyAndGradient(x, gradient);
            } catch (EnergyException ex) {
                writeStoredSnapshots();
                throw ex;
            }

            // Add the potential energy of the slow degrees of freedom.
            currentPotentialEnergy += respaE.getAsDouble();

            detectAtypicalEnergy(priorPE, defaultDeltaPEThresh);

            // Do the full-step integration operation.
            integrator.postForce(gradient);

            // Compute the full-step kinetic energy.
            thermostat.computeKineticEnergy();

            // Do the full-step thermostat operation.
            thermostat.fullStep(dt);

            // Recompute the kinetic energy after the full-step thermostat operation.
            thermostat.computeKineticEnergy();

            // Remove center of mass motion every ~100 steps.
            // TODO: Consider making this another MDListener.
            int removeCOMMotionFrequency = 100;
            if (removeCOMM && step % removeCOMMotionFrequency == 0) {
                thermostat.centerOfMassMotion(true, false);
            }

            // Collect current kinetic energy, temperature, and total energy.
            currentKineticEnergy = thermostat.getKineticEnergy();
            currentTemperature = thermostat.getCurrentTemperature();
            currentTotalEnergy = currentKineticEnergy + currentPotentialEnergy;

            // Update atomic velocity, acceleration and previous acceleration.
            potential.setVelocity(v);
            potential.setAcceleration(a);
            potential.setPreviousAcceleration(aPrevious);

            // Log the current state every printFrequency steps.
            totalSimTime += dt;
            final long currStep = step; // For streaming purposes.
            interiorListeners.forEach((MDListener mdL) -> mdL.actionForStep(currStep));

            // Check for a termination request.
            if (terminate) {
                logger.info(format("\n Terminating after %8d time steps\n", step));
                break;
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void run() {
        try {
            preRunOps();
        } catch (IllegalStateException ise) {
            return;
        }
        initializeEnergies();
        postInitEnergies();
        mainLoop();
        postRun();
    }

    /**
     * Detects grossly atypical potential energy values that are likely
     * incorrect, and writes snapshots to disc. "Grossly atypical" is defined
     * as: greater than 1.0E100 kcal/mol, less than -1.0E100 kcal/mol,
     * non-finite (NaN/infinite), or exceeding specified delta from the prior
     * potential energy.
     * <p>
     * After prior snapshots have been written to disc, the queue they are
     * stored on is now empty, preventing printing of duplicate snapshots.
     *
     * @param priorPE     Potential energy prior to this step.
     * @param delPEThresh If potential energy changes by this much, trigger a
     */
    void detectAtypicalEnergy(double priorPE, double delPEThresh) {

        // If not keeping snapshots, disable functionality.
        if (numSnapshotsToKeep < 1) {
            return;
        }

        double deltaPE = currentPotentialEnergy - priorPE;

        CoordinateSnapshot currState = new CoordinateSnapshot();
        currState.storeState();
        lastSnapshots.add(currState);

        double maxPEThresh = 1.0E100; // 1.0E100 kcal/mol is well into the territory of the absurd.
        double absPE = Math.abs(currentPotentialEnergy);

        if (absPE > maxPEThresh || !Double.isFinite(currentPotentialEnergy) || Math.abs(deltaPE) > delPEThresh) {
            logger.info(format(" Unusual potential energy %12.5g detected, writing snapshots.", currentPotentialEnergy));
            writeStoredSnapshots();
            currState.revertState(); // May be unnecessary, thanks to the current state always being last on the queue.
            if (absPE > 1.0E100 || !Double.isFinite(currentPotentialEnergy)) {
                logger.severe(format(" Dynamics exiting with atypical potential energy of %12.5g", currentPotentialEnergy));
            }
        }
    }

    /**
     * Performs the inner loop of writing snapshots to disk; used by both
     * detectAtypicalEnergy and a try-catch in dynamics.
     */
    public void writeStoredSnapshots() {
        int numSnaps = lastSnapshots.size();

        File origFile = molecularAssembly.getFile();
        String timeString = LocalDateTime.now().format(DateTimeFormatter.
                ofPattern("HH_mm_ss"));
        PotentialsFunctions potentialsFunctions = new PotentialsUtils();

        String filename = format("%s-%s-SNAP.pdb",
                FilenameUtils.removeExtension(molecularAssembly.getFile().getName()),
                timeString);

        for (int is = 0; is < numSnaps; is++) {
            CoordinateSnapshot oldState = lastSnapshots.poll();
            if (oldState != null) {
                oldState.revertState();
            }
            potentialsFunctions.saveAsPDB(molecularAssembly, new File(potentialsFunctions.versionFile(filename)));
        }
        molecularAssembly.setFile(origFile);
    }

    public MolecularAssembly[] getAssemblies() {
        return Arrays.stream(assemblies).
                map(AssemblyInfo::getAssembly).
                toArray(MolecularAssembly[]::new);
    }

    public void setTrajectoryFiles(File[] outputFiles) {
        logger.info(" Setting files " + Arrays.toString(outputFiles));
        int nFi = assemblies.length;
        assert outputFiles.length == nFi;
        for (int i = 0; i < nFi; i++) {
            AssemblyInfo ai = assemblies[i];
            ai.setArchiveFile(outputFiles[i]);
        }
    }

    /**
     * Get the total system energy (kinetic plus potential).
     *
     * @return total energy.
     */
    public double getTotalEnergy() {
        return currentTotalEnergy;
    }

    /**
     * Get the system kinetic energy.
     *
     * @return kinetic energy.
     */
    public double getKineticEnergy() {
        return currentKineticEnergy;
    }

    /**
     * Get the system potential energy.
     *
     * @return potential energy.
     */
    public double getPotentialEnergy() {
        return currentPotentialEnergy;
    }

    /**
     * Get the current temperature of the system
     *
     * @return currentTemperature
     */
    public double getTemperature() {
        return currentTemperature;
    }

    /**
     * Gets the kinetic energy at the start of the last dynamics run.
     *
     * @return Kinetic energy at the start of the run.
     */
    public double getInitialKineticEnergy() {
        return initialKinetic;
    }

    /**
     * Gets the temperature at the start of the last dynamics run.
     *
     * @return temperature at the start of the run.
     */
    public double getInitialTemperature() {
        return initialTemp;
    }

    /**
     * Gets the potential energy at the start of the last dynamics run.
     *
     * @return potential energy at the start of the run.
     */
    public double getInitialPotentialEnergy() {
        return initialPotential;
    }

    /**
     * Gets the total energy at the start of the last dynamics run.
     *
     * @return total energy at the start of the run.
     */
    public double getInitialTotalEnergy() {
        return initialTotal;
    }

    /**
     * Returns the associated dynamics file.
     *
     * @return Dynamics restart File.
     */
    public File getDynFile() {
        return restartFile;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, " Exception terminating dynamics.\n", e);
                }
            }
        }
    }

    private double[] copyArray(double[] array, double[] source) {
        if (array == null || array.length < numberOfVariables) {
            return Arrays.copyOf(source, numberOfVariables);
        } else {
            System.arraycopy(source, 0, array, 0, numberOfVariables);
            return array;
        }
    }

    @Override
    public double[] getX() {
        return x;
    }

    @Override
    public double[] copyX(double[] x) {
        return copyArray(x, this.x);
    }

    @Override
    public double[] getV() {
        return v;
    }

    @Override
    public double[] copyV(double[] v) {
        return copyArray(v, this.v);
    }

    @Override
    public double[] getA() {
        return a;
    }

    @Override
    public double[] copyA(double[] a) {
        return copyArray(a, this.a);
    }

    @Override
    public double[] getPriorAcceleration() {
        return aPrevious;
    }

    @Override
    public double[] copyPriorAcceleration(double[] aPrev) {
        return copyArray(aPrev, this.aPrevious);
    }

    @Override
    public boolean topLevelAlgorithm() {
        return automaticWriteouts;
    }

    @Override
    public Optional<Crystal> getCrystal() {
        return (potential instanceof CrystalPotential) ? Optional.of(((CrystalPotential) potential).getCrystal()) : Optional.empty();
    }

    @Override
    public boolean actionForStep(long step) {
        // TODO: Act only on "outside" listeners.
        boolean someAction = false;
        for (MDListener listener : allListeners) {
            if (listener.actionForStep(step)) {
                someAction = true;
            }
        }
        return someAction;
    }

    /**
     * A simple container class to hold all the infrastructure associated with a
     * MolecularAssembly for MolecularDynamics; assembly, properties, archive
     * and PDB files, PDB and XYZ filters. Direct access to package-private
     * members breaks encapsulation a bit, but the private inner class shouldn't
     * be used externally anyways.
     */
    protected static class AssemblyInfo {

        private final MolecularAssembly assembly;
        //private final CompositeConfiguration compositeConfiguration;
        File archiveFile = null;
        File pdbFile;
        PDBFilter pdbFilter;
        XYZFilter xyzFilter = null;
        TrajectoryListener trjL = null;

        AssemblyInfo(MolecularAssembly assembly) {
            this.assembly = assembly;
            pdbFile = this.assembly.getFile();
            //compositeConfiguration = this.assembly.getProperties();
            pdbFilter = new PDBFilter(this.assembly.getFile(), this.assembly,
                    this.assembly.getForceField(), this.assembly.getProperties());
        }

        public MolecularAssembly getAssembly() {
            return assembly;
        }

        void setArchiveFile(File file) {
            archiveFile = file;
            trjL.alterFile(file);
        }

        void setTrajectoryListener(TrajectoryListener trjL) {
            this.trjL = trjL;
        }
    }

    /**
     * <p>
     * storeState.</p>
     */
    public void storeState() {
        if (dynamicsState == null) {
            dynamicsState = new DynamicsState();
        }
        dynamicsState.storeState();
    }

    /**
     * <p>
     * revertState.</p>
     *
     * @throws java.lang.Exception if any.
     */
    public void revertState() throws Exception {
        if (dynamicsState == null) {
            throw new Exception();
        }
        dynamicsState.revertState();
    }

    /**
     * More limited version of a DynamicsState, storing only coordinates.
     * TODO: Make DynamicsState more flexible and let it store any combination of variables.
     */
    protected class CoordinateSnapshot {
        final double[] xBak;

        CoordinateSnapshot() {
            xBak = new double[numberOfVariables];
        }

        void storeState() {
            arraycopy(x, 0, xBak, 0, numberOfVariables);
        }

        void revertState() {
            arraycopy(xBak, 0, x, 0, numberOfVariables);
            Atom[] atoms = molecularAssembly.getActiveAtomArray();
            for (int i = 0; i < atoms.length; i++) {
                int i3 = 3 * i;
                double[] newXYZ = new double[3];
                arraycopy(xBak, i3, newXYZ, 0, 3);
                atoms[i].setXYZ(newXYZ);
            }
        }
    }

    protected class DynamicsState {

        double[] xBak, vBak, aBak;
        double[] aPreviousBak, massBak, gradBak;
        double currentKineticEnergyBak, currentPotentialEnergyBak, currentTotalEnergyBak;
        double currentTemperatureBak;

        DynamicsState() {
            xBak = new double[numberOfVariables];
            vBak = new double[numberOfVariables];
            aBak = new double[numberOfVariables];
            aPreviousBak = new double[numberOfVariables];
            massBak = new double[numberOfVariables];
            gradBak = new double[numberOfVariables];
        }

        public void storeState() {
            currentKineticEnergyBak = currentKineticEnergy;
            currentPotentialEnergyBak = currentPotentialEnergy;
            currentTotalEnergyBak = currentTotalEnergy;
            currentTemperatureBak = currentTemperature;
            arraycopy(x, 0, xBak, 0, numberOfVariables);
            arraycopy(v, 0, vBak, 0, numberOfVariables);
            arraycopy(a, 0, aBak, 0, numberOfVariables);
            arraycopy(aPrevious, 0, aPreviousBak, 0, numberOfVariables);
            arraycopy(mass, 0, massBak, 0, numberOfVariables);
            arraycopy(gradient, 0, gradBak, 0, numberOfVariables);
        }

        public void describe(String title) {
            StringBuilder sb = new StringBuilder();
            sb.append(title);
            sb.append("\nx: ");
            Arrays.stream(x).forEach(val -> sb.append(format("%.2g, ", val)));
            sb.append("\nv: ");
            Arrays.stream(v).forEach(val -> sb.append(format("%.2g, ", val)));
            sb.append("\na: ");
            Arrays.stream(a).forEach(val -> sb.append(format("%.2g, ", val)));
            sb.append("\naP: ");
            Arrays.stream(aPrevious).forEach(val -> sb.append(format("%.2g, ", val)));
            sb.append("\nm: ");
            Arrays.stream(mass).forEach(val -> sb.append(format("%.2g, ", val)));
            sb.append("\ng: ");
            Arrays.stream(gradient).forEach(val -> sb.append(format("%.2g, ", val)));
            sb.append(format("\nK,U,E,T: %g %g %g %g\n",
                    currentKineticEnergy, currentPotentialEnergy,
                    currentTotalEnergy, currentTemperature));
            logger.info(sb.toString());
        }

        public void revertState() {
            currentKineticEnergy = currentKineticEnergyBak;
            currentPotentialEnergy = currentPotentialEnergyBak;
            currentTotalEnergy = currentTotalEnergyBak;
            currentTemperature = currentTemperatureBak;
            arraycopy(xBak, 0, x, 0, numberOfVariables);
            arraycopy(vBak, 0, v, 0, numberOfVariables);
            arraycopy(aBak, 0, a, 0, numberOfVariables);
            arraycopy(aPreviousBak, 0, aPrevious, 0, numberOfVariables);
            arraycopy(massBak, 0, mass, 0, numberOfVariables);
            arraycopy(gradBak, 0, gradient, 0, numberOfVariables);

            Atom[] atoms = molecularAssembly.getActiveAtomArray();
            if (atoms.length * 3 == numberOfVariables) {
                double[] vec = new double[3];
                int index = 0;
                for (Atom atom : atoms) {
                    atom.moveTo(x[index], x[index + 1], x[index + 2]);
                    atom.setXYZGradient(gradient[index], gradient[index + 1], gradient[index + 2]);
                    vec[0] = v[index];
                    vec[1] = v[index + 1];
                    vec[2] = v[index + 2];
                    atom.setVelocity(vec);
                    vec[0] = a[index];
                    vec[1] = a[index + 1];
                    vec[2] = a[index + 2];
                    atom.setAcceleration(vec);
                    vec[0] = aPrevious[index];
                    vec[1] = aPrevious[index + 1];
                    vec[2] = aPrevious[index + 2];
                    atom.setPreviousAcceleration(vec);
                    index += 3;
                }
            }
        }
    }

    /**
     * No-op; FFX does not need to occasionally return information from FFX.
     *
     * @param intervalSteps Ignored.
     */
    public void setIntervalSteps(int intervalSteps) {
        if (initialized) {
            logger.warning(" Attempting to set interval steps for an already-initialized MD!");
        }
    }

    /**
     * <p>
     * getTimeStep.</p>
     *
     * @return Timestep in picoseconds.
     */
    public double getTimeStep() {
        return dt;
    }

    /**
     * <p>
     * getIntervalSteps.</p>
     *
     * @return Always 1 for this implementation.
     */
    public int getIntervalSteps() {
        return 1;
    }

    public enum VerbosityLevel {
        VERBOSE(false, Level.INFO, Level.INFO),
        QUIET(true, Level.INFO, Level.FINE),
        SILENT(true, Level.FINE, Level.FINE);

        public static final VerbosityLevel DEFAULT_VERBOSITY = VERBOSE;
        private boolean isQuiet;
        private final Level basicLevel;
        private final Level intmdLevel;

        VerbosityLevel(boolean isQuiet, Level basic, Level intmd) {
            this.isQuiet = isQuiet;
            this.basicLevel = basic;
            this.intmdLevel = intmd;
        }

        public boolean isQuiet() {
            return isQuiet;
        }
    }

    /**
     * Monte Carlo notification enumeration.
     */
    public enum MonteCarloNotification {
        NEVER, EACH_STEP, AFTER_DYNAMICS
    }

    /**
     * Enumerates available molecular dynamics engines; presently limited to the
     * FFX reference engine and the OpenMM engine.
     * <p>
     * Distinct from the force field energy Platform, as the FFX engine can use
     * OpenMM energies, but not vice-versa.
     */
    public enum DynamicsEngine {
        FFX(true, true), OPENMM(false, true);

        // Set of supported Platforms. The EnumSet paradigm is very efficient, as it
        // is internally stored as a bit field.
        private final EnumSet<ForceFieldEnergy.Platform> platforms = EnumSet.noneOf(ForceFieldEnergy.Platform.class);

        /**
         * Constructs a DynamicsEngine using the two presently known types of
         * Platform.
         *
         * @param ffx    Add support for the FFX reference energy platform.
         * @param openMM Add support for the OpenMM energy platforms.
         */
        DynamicsEngine(boolean ffx, boolean openMM) {
            if (ffx) {
                platforms.add(ForceFieldEnergy.Platform.FFX);
            }
            if (openMM) {
                platforms.add(ForceFieldEnergy.Platform.OMM);
                platforms.add(ForceFieldEnergy.Platform.OMM_REF);
                platforms.add(ForceFieldEnergy.Platform.OMM_CUDA);
                platforms.add(ForceFieldEnergy.Platform.OMM_OPENCL);
                platforms.add(ForceFieldEnergy.Platform.OMM_OPTCPU);
            }
        }

        /**
         * Checks if this energy Platform is supported by this DynamicsEngine
         *
         * @param platform The requested platform.
         * @return If supported
         */
        public boolean supportsPlatform(ForceFieldEnergy.Platform platform) {
            return platforms.contains(platform);
        }

        /**
         * Gets the set of Platforms supported by this DynamicsEngine
         *
         * @return An EnumSet
         */
        public EnumSet<ForceFieldEnergy.Platform> getSupportedPlatforms() {
            return EnumSet.copyOf(platforms);
        }
    }
}
