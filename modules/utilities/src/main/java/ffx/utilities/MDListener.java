package ffx.utilities;

public interface MDListener {
    /**
     * May perform an action associated with an MD step based on
     * the listener's criteria. For example, if snapshots are to
     * be written every 500 steps, it will write every 500 steps.
     * @param step MD step update.
     * @return If an action was performed.
     */
    boolean actionForStep(long step);
}
