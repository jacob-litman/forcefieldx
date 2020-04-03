package ffx.algorithms.listeners;

import ffx.utilities.MDListener;

import java.util.ArrayList;
import java.util.List;

public class CompoundListener extends IntervalListener {
    private final List<Runnable> priorActions;
    private final List<MDListener> listeners;

    public CompoundListener(long frequency, long priorSteps, List<Runnable> priorActions, List<MDListener> listeners) {
        super(frequency, priorSteps);
        this.priorActions = new ArrayList<>(priorActions);
        this.listeners = new ArrayList<>(listeners);
    }

    @Override
    public boolean actionForStep(final long step) {
        if (checkWrite(step)) {
            priorActions.forEach(Runnable::run);
            boolean allActed = listeners.stream().allMatch((MDListener mdL) -> mdL.actionForStep(step));
            assert allActed;
            return allActed;
        }
        return false;
    }
}
