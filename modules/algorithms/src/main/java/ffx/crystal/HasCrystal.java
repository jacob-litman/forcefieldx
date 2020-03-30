package ffx.crystal;

import java.util.Optional;

/**
 * Describes something that can have a Crystal.
 */
public interface HasCrystal {
    /**
     * The Crystal involved or null.
     *
     * @return Either a Crystal or null.
     */
    Optional<Crystal> getCrystal();
}
