//******************************************************************************
//
// File:    SharedInteger.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.SharedInteger
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// ark@cs.rit.edu.
//
// This Java source file is part of the Parallel Java Library ("PJ"). PJ is free
// software; you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// PJ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a module
// which is not derived from or based on this library. If you modify this library,
// you may extend this exception to your version of the library, but you are not
// obligated to do so. If you do not wish to do so, delete this exception
// statement from your version.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.pj.reduction;

import java.io.Serial;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Class SharedInteger provides a reduction variable for a value of type
 * <code>int</code>.
 * <P>
 * Class SharedInteger is multiple thread safe. The methods use lock-free atomic
 * compare-and-set.
 * <P>
 * <I>Note:</I> Class SharedInteger is implemented using class
 * java.util.concurrent.atomic.AtomicInteger.
 *
 * @author Alan Kaminsky
 * @version 06-Jun-2007
 */
public class SharedInteger
        extends Number {

    @Serial
    private static final long serialVersionUID = 1L;

// Hidden data members.
    private AtomicInteger myValue;

// Exported constructors.
    /**
     * Construct a new integer reduction variable with the initial value 0.
     */
    public SharedInteger() {
        myValue = new AtomicInteger();
    }

    /**
     * Construct a new integer reduction variable with the given initial value.
     *
     * @param initialValue Initial value.
     */
    public SharedInteger(int initialValue) {
        myValue = new AtomicInteger(initialValue);
    }

// Exported operations.
    /**
     * Returns this reduction variable's current value.
     *
     * @return Current value.
     */
    public int get() {
        return myValue.get();
    }

    /**
     * Set this reduction variable to the given value.
     *
     * @param value New value.
     */
    public void set(int value) {
        myValue.set(value);
    }

    /**
     * Set this reduction variable to the given value and return the previous
     * value.
     *
     * @param value New value.
     * @return Previous value.
     */
    public int getAndSet(int value) {
        return myValue.getAndSet(value);
    }

    /**
     * Atomically set this reduction variable to the given updated value if the
     * current value equals the expected value.
     *
     * @param expect Expected value.
     * @param update Updated value.
     * @return True if the update happened, false otherwise.
     */
    public boolean compareAndSet(int expect,
            int update) {
        return myValue.compareAndSet(expect, update);
    }

    /**
     * Atomically set this reduction variable to the given updated value if the
     * current value equals the expected value. May fail spuriously.
     *
     * @param expect Expected value.
     * @param update Updated value.
     * @return True if the update happened, false otherwise.
     */
    @SuppressWarnings("deprecation")
    public boolean weakCompareAndSet(int expect,
            int update) {
        return myValue.weakCompareAndSet(expect, update);
    }

    /**
     * Add one to this reduction variable and return the previous value.
     *
     * @return Previous value.
     */
    public int getAndIncrement() {
        return myValue.getAndIncrement();
    }

    /**
     * Subtract one from this reduction variable and return the previous value.
     *
     * @return Previous value.
     */
    public int getAndDecrement() {
        return myValue.getAndDecrement();
    }

    /**
     * Add the given value to this reduction variable and return the previous
     * value.
     *
     * @param value Value to add.
     * @return Previous value.
     */
    public int getAndAdd(int value) {
        return myValue.getAndAdd(value);
    }

    /**
     * Add one to this reduction variable and return the new value.
     *
     * @return New value.
     */
    public int incrementAndGet() {
        return myValue.incrementAndGet();
    }

    /**
     * Subtract one from this reduction variable and return the new value.
     *
     * @return New value.
     */
    public int decrementAndGet() {
        return myValue.decrementAndGet();
    }

    /**
     * Add the given value to this reduction variable and return the new value.
     *
     * @param value Value to add.
     * @return New value.
     */
    public int addAndGet(int value) {
        return myValue.addAndGet(value);
    }

    /**
     * Combine this reduction variable with the given value using the given
     * operation. The result is stored back into this reduction variable and is
     * returned.
     *
     * @param value Value.
     * @param op Binary operation.
     * @return (This variable) <I>op</I> (<code>value</code>).
     */
    public int reduce(int value,
            IntegerOp op) {
        for (;;) {
            int oldvalue = myValue.get();
            int newvalue = op.op(oldvalue, value);
            if (myValue.compareAndSet(oldvalue, newvalue)) {
                return newvalue;
            }
        }
    }

    /**
     * Returns a string version of this reduction variable.
     *
     * @return String version.
     */
    public String toString() {
        return Integer.toString(get());
    }

    /**
     * Returns this reduction variable's current value converted to type
     * <code>int</code>.
     *
     * @return Current value.
     */
    public int intValue() {
        return get();
    }

    /**
     * Returns this reduction variable's current value converted to type
     * <code>long</code>.
     *
     * @return Current value.
     */
    public long longValue() {
        return get();
    }

    /**
     * Returns this reduction variable's current value converted to type
     * <code>float</code>.
     *
     * @return Current value.
     */
    public float floatValue() {
        return (float) get();
    }

    /**
     * Returns this reduction variable's current value converted to type
     * <code>double</code>.
     *
     * @return Current value.
     */
    public double doubleValue() {
        return get();
    }

}
