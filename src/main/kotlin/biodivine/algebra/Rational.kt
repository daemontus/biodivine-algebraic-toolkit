package biodivine.algebra

/*
import kotlin.math.max

/**
 * Rational number is a fraction of two infinite precision integral numbers.
 *
 * Ideally, the two numbers should always be co-prime, however, in order to avoid expensive
 * division/remainder operations, we can defer this type of normalisation for later.
 *
 * Note that both the numerator and denominator can be null. However, the default value
 * for numerator is 0 whereas the default value for the denominator is 1.
 *
 * In particular, this means that if numerator is null, we can safely set denominator to null
 * as well (0/x = 0/y). Of denominator is null, the number is an integer.
 */
class Rational(
    val sign: Int,
    val numerator: IntArray?,
    val denominator: IntArray?
) {

    companion object {
        val ZERO = Rational(0, null, null)
        val EMPTY_ARRAY = IntArray(0)
    }

    /**
     * a   c   a*d + b*c
     * - + - = ---------
     * b   d      b*d
     */
    operator fun plus(that: Rational): Rational = when {
        this.sign == 0 || this.numerator == null -> that
        that.sign == 0 || that.numerator == null -> this
        this.denominator == null && that.denominator == null -> {
            // integer addition
            val newNumerator = IntArray(max(
                IntArrayMath.cardinality(this.numerator), IntArrayMath.cardinality(that.numerator)
            ) + 1)
            IntArrayMath.addition(this.numerator, that.numerator, newNumerator)
            Rational(this.sign * that.sign, newNumerator, null)
        }
        this.denominator == null && that.denominator != null -> that.plus(this)
        this.denominator != null && that.denominator == null -> {
            // fraction plus integer
            val bc = IntArray(IntArrayMath.cardinality(that.numerator) + IntArrayMath.cardinality(this.denominator))
            IntArrayMath.multiplication(this.denominator, that.numerator, bc)
            val newNumerator = IntArray(max(IntArrayMath.cardinality(this.numerator), bc.size) + 1)
            IntArrayMath.addition(this.numerator, bc, newNumerator)
            Rational(this.sign * that.sign, newNumerator, this.denominator)
        }
        this.denominator != null && that.denominator != null -> {
            val ad = IntArray(IntArrayMath.cardinality(this.numerator) + IntArrayMath.cardinality(that.denominator))
            val bc = IntArray(IntArrayMath.cardinality(that.numerator) + IntArrayMath.cardinality(this.denominator))
            IntArrayMath.multiplication(this.numerator, that.denominator, ad)
            IntArrayMath.multiplication(that.numerator, this.denominator, bc)
            val newNumerator = IntArray(max(bc.size, ad.size) + 1)
            IntArrayMath.addition(ad, bc, newNumerator)
            val newDenominator = IntArray(IntArrayMath.cardinality(this.denominator) + IntArrayMath.cardinality(that.denominator))
            IntArrayMath.multiplication(this.denominator, that.denominator, newDenominator)
            Rational(this.sign, newNumerator, newDenominator)
        }
        else -> error("WTF?")
    }

}
*/