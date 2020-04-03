package ffx.algorithms.listeners;

import ffx.utilities.MDListener;

import java.util.ArrayList;
import java.util.List;

public abstract class IntervalListener implements MDListener {
    public final long frequency;
    private long nextWrite;

    public IntervalListener(long frequency, long priorSteps) {
        assert frequency > 0 && priorSteps >= 0;
        this.frequency = frequency;
        setNextWrite(priorSteps);
    }

    /**
     * Informs this IntervalListener that this many steps have
     * been taken, and that it should update its next action
     * time accordingly.
     *
     * @param step Steps that have been taken.
     */
    public void setNextWrite(long step) {
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
