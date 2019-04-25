package biodivine.algebra.synth

import biodivine.algebra.NumQ
import biodivine.algebra.ia.*
import cc.redberry.rings.Rational
import cc.redberry.rings.Rings
import java.math.BigInteger
import java.util.*

/**
 * A box represents one state in our system. Note that while the usage
 * of rationals may seem costly, there aren't going to be that many states
 * in the end (and there will be very little manipulation going on there)
 */
class Box(
    vararg items: Interval
) {

    val data: Array<Interval> = items.map { it }.toTypedArray() // wtf polymorphism

    /**
     * Compute volume of the box. Remember that this grows proportionally
     * with the number of dimensions! So a box [10x10] is 10x smaller than
     * [10x10x10] (and reverse for values <1)
     */
    val volume: NumQ
        get() = data.fold(Rings.Q.one) { a, (low, high) ->
            a * (high - low)
        }

    /**
     * Split this box into 2^dim boxes of half the size.
     */
    fun subdivide(): List<Box> {
        val dimensions = data.size
        val result = ArrayList<Box>()
        for (mask in 0 until (1.shl(dimensions))) {
            // mask determines which part of the subdivided interval should be taken - 1=high, 0=low
            val newBox = Box(*Array(dimensions) { d ->
                val (low, high) = data[d]
                val half = (high - low) / 2
                if (mask.shr(d).and(1) == 1) {  // mask is one, take high interval
                    Interval(high - half, high)
                } else {                              // mask is zero, take low interval
                    Interval(low, low + half)
                }
            })
            result.add(newBox)
        }
        return result
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as Box

        if (!data.contentEquals(other.data)) return false

        return true
    }

    override fun hashCode(): Int {
        return data.contentHashCode()
    }

    override fun toString(): String {
        return "Box(data=${Arrays.toString(data)})"
    }


}