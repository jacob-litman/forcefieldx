package ffx.algorithms.listeners;

import ffx.algorithms.dynamics.DynamicAlgorithm;
import ffx.potential.parsers.DYNFilter;

import java.io.File;

public class DynListener extends MDFileListener {
    private final DYNFilter dynFilter;
    private final DynamicAlgorithm algo;

    public DynListener(File dyn, DYNFilter filter, long frequency, long priorSteps, DynamicAlgorithm algorithm) {
        super(dyn, frequency, priorSteps);
        dynFilter = filter;
        algo = algorithm;
    }

    @Override
    public boolean actionForStep(long step) {
        if (checkWrite(step)) {
            /*
            potential.writeAdditionalRestartInfo(true);
        String dynName = FileUtils.relativePathTo(restartFile).toString();
        if (dynFilter.writeDYN(restartFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
            logger.log(basicLogging, " Wrote dynamics restart file to " + dynName);
        } else {
            logger.log(basicLogging, " Writing dynamics restart file to " + dynName + " failed");
        }
             */

            dynFilter.writeDYN(outputFile, algo.getCrystal().orElse(null),
                    algo.getX(), algo.getV(), algo.getA(), algo.getPriorAcceleration());
            return true;
        }
        return false;
    }
}
