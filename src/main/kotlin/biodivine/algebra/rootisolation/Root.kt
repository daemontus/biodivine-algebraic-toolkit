package biodivine.algebra.rootisolation

import biodivine.algebra.UPoly
import biodivine.algebra.getNumberOfSignChanges
import biodivine.algebra.ia.*
import biodivine.algebra.transformPolyToInterval

/**
 * Represents one (irrational) root of a polynomial. The point of this class is that it allows infinite error
 * operations on roots, such as exact comparison.
 *
 * It achieves this by increasing the error of root isolation on the fly when needed. That is why we
 * do not expose the actual bounds as part of the class API - because they can change after each operation.
 *
 * This also makes this class thread unsafe!
 */
class Root private constructor(
    private val polynomial: UPoly,
    private var lowerBound: NumQ,
    private var upperBound: NumQ
) : Comparable<Root> {

    companion object {

        /**
         * Create a rational root from a linear polynomial
         */
        fun rational(polynomial: UPoly): Root {
            val root = polynomial[0].negate().divide(polynomial[1])
            return Root(polynomial, root, root)
        }

        /**
         * Create an irrational root
         */
        fun irrational(low: NumQ, high: NumQ, polynomial: UPoly): Root {
            return Root(polynomial, low, high)
        }

    }

    init {
        if (lowerBound > upperBound) error("Invalid root! Interval [$lowerBound, $upperBound] is empty.")
    }

    private val error: NumQ
        get() = upperBound - lowerBound

    fun middleValue(other: Root): NumQ {
        if (this >= other) error("Cannot find middle value in [$this, $other] (empty interval)")
        // at this point, both roots have been refined to be comparable!
        return this.upperBound + (other.lowerBound - this.upperBound) / 2
    }

    private fun refine() {
        val middlePoint = lowerBound + (error / 2)
        val signChangesInLower = polynomial.transformPolyToInterval(lowerBound, middlePoint).getNumberOfSignChanges()
        val signChangesInUpper = polynomial.transformPolyToInterval(middlePoint, upperBound).getNumberOfSignChanges()
        when {
            // sanity check
            signChangesInLower + signChangesInUpper != 1 -> error("This should not happen, but I can't prove it.")
            signChangesInLower == 1 -> upperBound = middlePoint
            signChangesInUpper == 1 -> lowerBound = middlePoint
        }
    }

    override operator fun compareTo(other: Root): Int {
        // first, refine the roots until they are comparable
        do {
            // find intersection of isolating intervals and check if the roots are comparable (also solved equality)
            val intersectionLow = max(this.lowerBound, other.lowerBound)
            val intersectionHigh = min(this.upperBound, other.upperBound)
            val comparable = when {
                // there is no intersection - roots are comparable
                intersectionLow > intersectionHigh -> true
                // there is an intersection and polynomials are different - we need more error
                this.polynomial != other.polynomial -> false
                // there is an intersection, but out interval is fully contained in the others
                intersectionLow == this.lowerBound && intersectionHigh == this.upperBound -> {
                    return 0    // shot-circuit: we know roots are equal, no need to work more
                }
                intersectionLow == other.lowerBound && intersectionHigh == other.upperBound -> {
                    return 1    // the same here, but for the other interval
                }
                // there is a non-trivial intersection which needs to be resolved
                else -> false
            }
            // if they are not comparable, refine the one with higher approximation error
            if (!comparable) {
                if (this.error > other.error) this.refine() else other.refine()
            }
        } while (!comparable)

        // At this point, we know that the roots are not equal and that they are comparable (have disjoint intervals)
        // so we can just test which bound is lower
        return this.upperBound.compareTo(other.lowerBound)
    }

    operator fun compareTo(other: NumQ): Int {
        // This is the same principle as in the compareTo method for other roots, but we know the other value is a
        // rational number, so intersections and logic are simpler.

        if (lowerBound == upperBound) { // if this is also a rational number, just compare them
            return lowerBound.compareTo(other)
        }

        do {
            val comparable = other < lowerBound || upperBound < other
            if (!comparable) this.refine()
        } while (!comparable)

        // here, we know that the number is outside of the interval, so it does not matter which bound we compare it to
        return lowerBound.compareTo(other)
    }

    override fun equals(other: Any?): Boolean {
        if (other !is Root) return false
        if (this.polynomial != other.polynomial) return false   // not the same root!
        // We can just use compare to, because there isn't much performance gain there - there is always a chance
        // we would have to refine the isolating interval
        return this.compareTo(other) == 0
    }

    override fun hashCode(): Int {
        // Unfortunately, we can't really use the isolating interval inside the hash-code, because
        // it can be different for equivalent roots, which would break the fact that two equivalent
        // objects need to have the same hash code. However, we at least know that when the polynomials
        // are different, the roots are different (this also works nicely for rational roots).
        return polynomial.hashCode()
    }

    override fun toString(): String {
        return "[$lowerBound, $upperBound]{$polynomial}"
    }


}