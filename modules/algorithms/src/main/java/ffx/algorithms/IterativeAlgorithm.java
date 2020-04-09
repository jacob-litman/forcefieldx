package ffx.algorithms;

import ffx.crystal.HasCrystal;
import ffx.utilities.MDListener;

import javax.annotation.Nullable;

/**
 * Describes an algorithm possessed of a crystal and coordinates.
 * Also extends MDListener to ping all of its listeners.
 */
public interface IterativeAlgorithm extends HasCrystal, MDListener {
    /**
     * Direct access to coordinates array.
     *
     * @return The coordinates array used.
     */
    double[] getX();

    /**
     * Copy of the coordinates array. If the passed array is
     * null or too short, returns a new array.
     *
     * @param x Array to be filled.
     * @return Copy of coordinates array.
     */
    double[] copyX(@Nullable double[] x);

    /**
     * Checks if this is the top-level algorithm. Should ideally
     * be set at construction time.
     *
     * @return Whether this is the top-level algorithm in use.
     */
    boolean topLevelAlgorithm();

    /**
     * Act on the "exterior" listeners too slow to fit in this
     * algorithm's cycle length. Interior listeners should be triggered
     * by the algorithm itself.
     *
     * @param step MD step update.
     * @return Whether any "exterior" listeners were triggered.
     */
    @Override
    boolean actionForStep(long step);
}
