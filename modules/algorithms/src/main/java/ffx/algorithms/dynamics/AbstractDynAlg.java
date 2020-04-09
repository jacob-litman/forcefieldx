package ffx.algorithms.dynamics;

import ffx.algorithms.listeners.IntervalListener;
import ffx.utilities.MDListener;

import java.util.ArrayList;
import java.util.List;

/**
 * DynamicAlgorithm with default implementation of distributeListeners.
 */
public abstract class AbstractDynAlg implements DynamicAlgorithm {
    private static final int DEFAULT_INIT_LISTENER_CAPACITY = 5; // Even 5 is overkill for the forseeable future.
    protected final List<MDListener> allListeners = new ArrayList<>(DEFAULT_INIT_LISTENER_CAPACITY);
    protected final List<MDListener> interiorListeners = new ArrayList<>(DEFAULT_INIT_LISTENER_CAPACITY);
    protected final List<MDListener> exteriorListeners = new ArrayList<>(DEFAULT_INIT_LISTENER_CAPACITY);

    @Override
    public void distributeListeners(long nSteps, long exteriorSteps) {
        interiorListeners.clear();
        exteriorListeners.clear();
        if (topLevelAlgorithm()) {
            interiorListeners.addAll(allListeners);
            interiorListeners.stream().
                    filter((MDListener mdL) -> mdL instanceof IntervalListener).
                    map((MDListener mdL) -> (IntervalListener) mdL).
                    forEach((IntervalListener il) -> il.setNextWrite(0));
        } else {
            for (MDListener listener : allListeners) {
                if (listener instanceof IntervalListener) {
                    IntervalListener il = (IntervalListener) listener;
                    if (il.frequency < nSteps) {
                        interiorListeners.add(listener);
                        il.setNextWrite(0);
                    } else {
                        exteriorListeners.add(listener);
                        il.setNextWrite(exteriorSteps);
                    }
                } else {
                    interiorListeners.add(listener);
                }
            }
        }
    }
}
