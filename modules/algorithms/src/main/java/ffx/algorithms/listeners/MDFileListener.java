package ffx.algorithms.listeners;

import ffx.utilities.MDListener;

import java.io.File;

public abstract class MDFileListener implements MDListener {
    protected File outputFile;
    protected final long frequency;
    private long stepsTaken;
    private long nextWrite;

    public MDFileListener(File fi, long frequency, long priorSteps) {
        assert frequency > 0 && priorSteps >= 0;
        outputFile = fi;
        this.frequency = frequency;
        stepsTaken = priorSteps;
        setNextWrite(priorSteps);
    }

    public void alterFile(File newOutput) {
        outputFile = newOutput;
    }

    private void setNextWrite(long step) {
        long cycles = step / frequency;
        nextWrite = frequency * (cycles + 1);
    }

    protected boolean checkWrite(long step) {
        if (step >= nextWrite) {
            setNextWrite(step);
            return true;
        }
        return false;
    }
}
