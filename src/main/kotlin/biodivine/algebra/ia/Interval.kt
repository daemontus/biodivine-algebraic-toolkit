package biodivine.algebra.ia

import biodivine.algebra.NumQ
import cc.redberry.rings.Rational
import cc.redberry.rings.Rings
import cc.redberry.rings.bigint.BigDecimal
import cc.redberry.rings.bigint.BigInteger
import kotlin.math.max
import kotlin.math.min


data class Interval(
    val low: NumQ, val high: NumQ
) {

    init {
        if (low > high) error("Empty interval [$low,$high]")
    }

    constructor(low: Int, high: Int) : this(
        NumQ(Rings.Q.ring, BigInteger.valueOf(low)),
        NumQ(Rings.Q.ring, BigInteger.valueOf(high))
    )

    val size: NumQ
        get() = high - low

    val hasZero: Boolean
        get() = low <= Rings.Q.zero && high >= Rings.Q.zero

    val center: NumQ
        get() = low + size/2

    operator fun plus(that: Interval): Interval {
        return Interval(this.low + that.low, this.high + that.high)
    }

    operator fun minus(that: Interval): Interval {
        return Interval(this.low - that.high, this.high - that.low)
    }

    operator fun times(that: Interval): Interval {
        val (x1, x2) = this
        val (y1, y2) = that
        val x1y1 = x1 * y1
        val x1y2 = x1 * y2
        val x2y1 = x2 * y1
        val x2y2 = x2 * y2
        return Interval(
            low = min(min(x1y1, x1y2), min(x2y1,x2y2)),
            high = max(max(x1y1, x1y2), max(x2y1,x2y2))
        )
    }

    operator fun div(that: Interval): Interval {
        if (that.low.signum() != that.high.signum()) error("Zero in division $this / $that")
        val (y1, y2) = that
        return this * Interval(
            NumQ(y2.ring, y2.denominator(), y2.numerator()),
            NumQ(y1.ring, y1.denominator(), y1.numerator())
        )
    }

    operator fun times(that: NumQ): Interval {
        return if (that < Rings.Q.zero) {
            Interval(high * that, low * that)
        } else {
            Interval(low * that, high * that)
        }
    }

    infix fun intersects(that: Interval): Boolean = this.high >= that.low && this.low <= that.high

    infix fun intersect(that: Interval): Interval? {
        val low = if (this.low < that.low) that.low else this.low
        val high = if (this.high < that.high) this.high else that.high
        return Interval(low, high)
    }

    operator fun contains(num: NumQ): Boolean = this.low <= num && this.high >= num

    fun isNumber(): Boolean = this.low == this.high
}

operator fun NumQ.plus(that: NumQ): NumQ = this.add(that)
operator fun NumQ.times(that: NumQ): NumQ = this.multiply(that)
operator fun NumQ.minus(that: NumQ): NumQ = this.subtract(that)
operator fun NumQ.times(that: Interval) = Interval(this * that.low, this * that.high)
operator fun NumQ.unaryMinus(): NumQ = this.negate()

operator fun NumQ.div(that: Int): NumQ = this.divide(BigInteger.valueOf(that))
operator fun NumQ.times(that: Int): NumQ = this.multiply(BigInteger.valueOf(that))
operator fun NumQ.div(that: NumQ) = this.divide(that)

fun min(a: NumQ, b: NumQ): NumQ {
    return if (a < b) a else b
}

fun max(a: NumQ, b: NumQ): NumQ {
    return if (a > b) a else b
}

fun NumQ.decimalString(): String {
    return BigDecimal(this.numerator()).divide(BigDecimal(this.denominator())).toString()
}