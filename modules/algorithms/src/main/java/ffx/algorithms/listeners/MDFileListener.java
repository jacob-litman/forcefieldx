package ffx.algorithms.listeners;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public abstract class MDFileListener extends IntervalListener {
    protected File outputFile;
    protected final List<Runnable> priorActions;

    public MDFileListener(File fi, long frequency, long priorSteps, List<Runnable> priorActions) {
        super(frequency, priorSteps);
        outputFile = fi;
        this.priorActions = new ArrayList<>(priorActions);
    }

    public void alterFile(File newOutput) {
        outputFile = newOutput;
    }
}
