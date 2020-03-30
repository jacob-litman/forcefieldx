package ffx.algorithms.dynamics;

import ffx.crystal.HasCrystal;

/**
 * Describes an algorithm possessed of a crystal, coordinates,
 * velocity, acceleration, and prior acceleration (e.g. MD).
 */
public interface DynamicAlgorithm extends HasCrystal {
    double[] getX();
    double[] copyX(double[] x);

    double[] getV();
    double[] copyV(double[] v);

    double[] getA();
    double[] copyA(double[] a);

    double[] getPriorAcceleration();
    double[] copyPriorAcceleration(double[] aPrev);
}
