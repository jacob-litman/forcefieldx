//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.mc;

import java.util.Random;
import java.util.function.DoubleSupplier;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;

import ffx.potential.bonded.LambdaInterface;

/**
 * Define an MC move to update lambda.
 *
 * @author Mallory R. Tollefson
 */
public class LambdaMove implements MCMove {

    private static final Logger logger = Logger.getLogger(LambdaMove.class.getName());
    public static final double DEFAULT_MOVE_SIZE = 0.01;
    public static final double MAX_CONTINOUOUS_MOVE = 0.1;

    /**
     * Current value of lambda, which always refreshed from the OST instance.
     */
    private double priorLambda;
    /**
     * Apply the Lambda move to a LambdaInterface.
     */
    private final LambdaInterface lamPotential;
    /**
     * Random number generator.
     */
    private final Random random;
    /**
     * Lambda move size:
     * 1) The standard deviation for continuous moves from a Gaussian distribution.
     * 2) The step size for discrete moves.
     */
    private final double moveSize;
    /**
     * If true, do continuous moves. Otherwise, use discrete moves.
     */
    private final boolean isContinuous;
    /**
     * Function to generate dL values.
     */
    private final DoubleSupplier dlGen;

    /**
     * <p>Constructor for LambdaMove.</p>
     *
     * @param lamPotential Lambda interface/potential to operate on.
     * @param continuous   If there is a continuous distribution of lambda values (vs. a discrete ladder).
     * @param moveSize     Either the distance between bins (discrete) or standard deviation of move size (continuous).
     * @param randomSeed   PRNG seed.
     */
    public LambdaMove(LambdaInterface lamPotential, boolean continuous, double moveSize, long randomSeed) {
        this(lamPotential, continuous, moveSize, new Random(randomSeed));
    }

    /**
     * Constructor for lambdaMove.
     *
     * @param lamPotential Lambda interface/potential to operate on.
     * @param continuous   If there is a continuous distribution of lambda values (vs. a discrete ladder).
     * @param moveSize     Either the distance between bins (discrete) or standard deviation of move size (continuous).
     * @param random       Source of randomness.
     */
    public LambdaMove(LambdaInterface lamPotential, boolean continuous, double moveSize, Random random) {
        this.lamPotential = lamPotential;
        this.moveSize = moveSize;
        this.isContinuous = continuous;
        dlGen = isContinuous ? this::continuousMove : this::discreteMove;
        this.random = random;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void move() {
        priorLambda = lamPotential.getLambda();

        // Draw a trial move from the distribution.
        double dL = dlGen.getAsDouble();
        double newLambda = mirror(priorLambda, dL);

        // Update the OST instance.
        lamPotential.setLambda(newLambda);
    }

    /**
     * Applies 0-1 mirroring conditions to lam + dL. Skips any moves where
     * dL is greater than 1 or less than -1, and skips 50% of moves from
     * 0 or 1 (exact).
     *
     * @param lam Initial lambda.
     * @param dL  Change in lambda.
     * @return Correctly mirrored lam + dL
     */
    private double mirror(double lam, double dL) {
        // Telescope to public static method because a public static method
        // may be useful in the future.
        return mirror(random, lam, dL);
    }

    /**
     * Applies 0-1 mirroring conditions to lam + dL. Skips any moves where
     * dL is greater than 1 or less than -1, and skips 50% of moves from
     * 0 or 1 (exact).
     *
     * @param random Source of randomness.
     * @param lam    Initial lambda.
     * @param dL     Change in lambda.
     * @return Correctly mirrored lam + dL
     */
    public static double mirror(Random random, double lam, double dL) {
        if (lam == 0.0 || lam == 1.0) {
            boolean skip = random.nextBoolean();
            if (skip) {
                return lam;
            }
        }
        // Eliminate really weird edge cases.
        if (abs(dL) > 1.0) {
            logger.warning(String.format(" Skipping large lambda move of %.3f not between -1 and +1", dL));
            return lam;
        }
        // Math.abs to mirror negative values.
        double newLam = abs(lam + dL);
        // If greater than 1, mirror via 2.0 - val
        return newLam <= 1.0 ? newLam : 2.0 - newLam;
    }

    /**
     * Validate lambda is in the range [0 .. 1].
     * <p>
     * For discrete moves, set Lambda to the closest valid value [0, dL, 2dL, .. 1].
     *
     * @param lambda Input lambda value.
     * @return Validated lambda value.
     */
    public double validateLambda(double lambda) {
        lambda = max(0.0, min(lambda, 1.0));
        if (isContinuous) {
            return lambda;
        }
        double remainder = lambda % moveSize;
        if (remainder < moveSize / 2.0) {
            return max(0.0, lambda - remainder);
        } else {
            return min(lambda + (moveSize - remainder), 1.0);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void revertMove() {
        lamPotential.setLambda(priorLambda);
    }

    /**
     * If true, do continuous moves. Otherwise, use discrete moves.
     */
    public boolean isContinuous() {
        return isContinuous;
    }

    /**
     * Pulls a delta-lambda from a continuous Gaussian distribution.
     *
     * @return A random Gaussian value with width of moveSize.
     */
    private double continuousMove() {
        // Draw a trial move from the distribution.
        return random.nextGaussian() * moveSize;
    }

    /**
     * Pulls a continuous lambda move of width moveSize.
     *
     * @return +/- moveSize (never 0)
     */
    private double discreteMove() {
        // Make a discrete move.
        double dL = moveSize;
        if (random.nextBoolean()) {
            dL = -moveSize;
        }
        return dL;
    }
}
