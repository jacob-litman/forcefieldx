package ffx.algorithms.listeners;

import ffx.potential.parsers.SystemFilter;
import ffx.utilities.FileUtils;

import java.io.File;
import java.util.List;
import java.util.function.DoubleSupplier;
import java.util.function.Supplier;
import java.util.logging.Level;
import java.util.logging.Logger;

public class TrajectoryListener extends MDFileListener {
    private static final Logger logger = Logger.getLogger(TrajectoryListener.class.getName());

    private String fileName;
    private final SystemFilter filter;
    private final Supplier<String[]> extraLines;
    private final Level logLevel;
    private final DoubleSupplier lamProvider;
    private final double minLambda;
    private final double maxLambda;

    public TrajectoryListener(File fi, SystemFilter filter, Supplier<String[]> extraLineSupplier,
                              long frequency, long priorSteps, List<Runnable> priorActions,
                              Level basicLogging) {
        this(fi, filter, extraLineSupplier, frequency, priorSteps, priorActions, basicLogging, () -> 1.0, 0.0, 1.0);
    }

    public TrajectoryListener(File fi, SystemFilter filter, Supplier<String[]> extraLineSupplier,
                              long frequency, long priorSteps, List<Runnable> priorActions,
                              Level basicLogging, DoubleSupplier lamProvider, double minLam, double maxLam) {
        super(fi, frequency, priorSteps, priorActions);
        this.filter = filter;
        fileName = FileUtils.relativePathTo(fi).toString();
        logLevel = basicLogging;
        this.minLambda = minLam;
        this.maxLambda = maxLam;
        this.lamProvider = lamProvider;
        this.extraLines = extraLineSupplier;
    }

    @Override
    public void alterFile(File newOutput) {
        super.alterFile(newOutput);
        fileName = FileUtils.relativePathTo(outputFile).toString();
    }

    @Override
    public boolean actionForStep(long step) {
        double lamVal = lamProvider.getAsDouble();
        // This may induce bunching up of snapshots.
        // TODO: Figure out if we want to append a line w/ step #
        if (lamVal <= maxLambda && lamVal >= minLambda && checkWrite(step)) {
            priorActions.forEach(Runnable::run);
            if (filter.writeFile(outputFile, true, extraLines.get())) {
                logger.log(logLevel, String.format(" Appended snap shot to %s", fileName));
                return true;
            } else {
                logger.warning(String.format(" Failed to append snap shot to %s", fileName));
                return false;
            }
        }
        return false;
    }
}
