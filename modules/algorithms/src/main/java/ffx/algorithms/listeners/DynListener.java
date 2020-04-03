package ffx.algorithms.listeners;

import ffx.algorithms.dynamics.DynamicAlgorithm;
import ffx.potential.parsers.DYNFilter;
import ffx.utilities.FileUtils;

import java.io.File;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

public class DynListener extends MDFileListener {
    private static final Logger logger = Logger.getLogger(DynListener.class.getName());
    private final DYNFilter dynFilter;
    private final DynamicAlgorithm algo;
    private final Level basicLogging;
    private String outFileName;

    public DynListener(File dyn, DYNFilter filter, List<Runnable> priorActions, long frequency,
                       long priorSteps, DynamicAlgorithm algorithm, Level basicLogging) {
        super(dyn, frequency, priorSteps, priorActions);
        dynFilter = filter;
        algo = algorithm;
        this.basicLogging = basicLogging;
        outFileName = FileUtils.relativePathTo(outputFile).toString();
    }

    @Override
    public boolean actionForStep(long step) {
        if (checkWrite(step)) {
            priorActions.forEach(Runnable::run);
            if (dynFilter.writeDYN(outputFile, algo.getCrystal().orElse(null),
                    algo.getX(), algo.getV(), algo.getA(), algo.getPriorAcceleration())) {
                logger.log(basicLogging, " Wrote dynamics restart file to " + outFileName);
                return true;
            } else {
                logger.log(basicLogging, " Writing dynamics restart file to " + outFileName + " failed");
                return false;
            }
        }
        return false;
    }

    @Override
    public void alterFile(File newOutput) {
        super.alterFile(newOutput);
        outFileName = FileUtils.relativePathTo(outputFile).toString();
    }
}
