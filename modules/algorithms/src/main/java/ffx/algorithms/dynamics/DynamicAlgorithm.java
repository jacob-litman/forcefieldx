package ffx.algorithms.dynamics;

import ffx.algorithms.listeners.IntervalListener;
import ffx.crystal.HasCrystal;
import ffx.utilities.MDListener;

import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.List;

/**
 * Describes an algorithm possessed of a crystal, coordinates,
 * velocity, acceleration, and prior acceleration (e.g. MD).
 * Also extends MDListener to ping all of its listeners.
 */
public abstract class DynamicAlgorithm implements HasCrystal, MDListener {
    private static final int DEFAULT_INIT_LISTENER_CAPACITY = 5; // Even 5 is overkill for the forseeable future.
    protected final List<MDListener> allListeners;
    protected final List<MDListener> interiorListeners;
    protected final List<MDListener> exteriorListeners;

    public DynamicAlgorithm() {
        this(DEFAULT_INIT_LISTENER_CAPACITY);
    }

    public DynamicAlgorithm(int listenerCapacity) {
        allListeners = new ArrayList<>(listenerCapacity);
        interiorListeners = new ArrayList<>(listenerCapacity);
        exteriorListeners = new ArrayList<>(listenerCapacity);
    }

    /**
     * Direct access to coordinates array.
     *
     * @return The coordinates array used.
     */
    public abstract double[] getX();

    /**
     * Copy of the coordinates array. If the passed array is
     * null or too short, returns a new array.
     *
     * @param x Array to be filled.
     * @return Copy of coordinates array.
     */
    public abstract double[] copyX(@Nullable double[] x);

    /**
     * Direct access to velocities array.
     *
     * @return The velocities array used.
     */
    public abstract double[] getV();

    /**
     * Copy of the velocities array. If the passed array is
     * null or too short, returns a new array.
     *
     * @param v Array to be filled.
     * @return Copy of velocities array.
     */
    public abstract double[] copyV(@Nullable double[] v);

    /**
     * Direct access to accelerations array.
     *
     * @return The accelerations array used.
     */
    public abstract double[] getA();

    /**
     * Copy of the accelerations array. If the passed array is
     * null or too short, returns a new array.
     *
     * @param a Array to be filled.
     * @return Copy of accelerations array.
     */
    public abstract double[] copyA(@Nullable double[] a);

    /**
     * Direct access to prior accelerations array.
     *
     * @return The prior accelerations array used.
     */
    public abstract double[] getPriorAcceleration();

    /**
     * Copy of the prior accelerations array. If the passed array is
     * null or too short, returns a new array.
     *
     * @param aPrev Array to be filled.
     * @return Copy of prior accelerations array.
     */
    public abstract double[] copyPriorAcceleration(@Nullable double[] aPrev);

    /**
     * Checks if this is the top-level algorithm. Should ideally
     * be set at construction time.
     *
     * @return Whether this is the top-level algorithm in use.
     */
    public abstract boolean topLevelAlgorithm();

    protected void distributeListeners(long nSteps, long exteriorSteps) {
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

    /**
     * Act on the "exterior" listeners too slow to fit in this
     * algorithm's cycle length. Interior listeners should be triggered
     * by the algorithm itself.
     *
     * @param step MD step update.
     * @return Whether any "exterior" listeners were triggered.
     */
    @Override
    public abstract boolean actionForStep(long step);
}
