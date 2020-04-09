package ffx.algorithms.dynamics;

import ffx.algorithms.IterativeAlgorithm;
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
public interface DynamicAlgorithm extends IterativeAlgorithm {

    /**
     * Direct access to velocities array.
     *
     * @return The velocities array used.
     */
    double[] getV();

    /**
     * Copy of the velocities array. If the passed array is
     * null or too short, returns a new array.
     *
     * @param v Array to be filled.
     * @return Copy of velocities array.
     */
    double[] copyV(@Nullable double[] v);

    /**
     * Direct access to accelerations array.
     *
     * @return The accelerations array used.
     */
    double[] getA();

    /**
     * Copy of the accelerations array. If the passed array is
     * null or too short, returns a new array.
     *
     * @param a Array to be filled.
     * @return Copy of accelerations array.
     */
    double[] copyA(@Nullable double[] a);

    /**
     * Direct access to prior accelerations array.
     *
     * @return The prior accelerations array used.
     */
    double[] getPriorAcceleration();

    /**
     * Copy of the prior accelerations array. If the passed array is
     * null or too short, returns a new array.
     *
     * @param aPrev Array to be filled.
     * @return Copy of prior accelerations array.
     */
    double[] copyPriorAcceleration(@Nullable double[] aPrev);
    
    void distributeListeners(long nSteps, long exteriorSteps);
}
