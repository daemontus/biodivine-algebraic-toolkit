package biodivine.algebra

import java.math.BigInteger
import java.text.ParseException

/**
 * Decimal number is a number with finite decimal representation:
 * 14582; 154.876; 3.14e15; 0.14e-2;
 *
 * The number consists of a variable integer prefix and a long exponent such that the actual value is given as:
 * value = prefix * 10^exponent
 *
 * Why not just use BigInteger or BigDecimal?
 *
 * BigInteger does not support decimal numbers, just integers,
 * furthermore, no exponent "compression" is performed, so 3.14e21 is represented as a very long array of ints.
 *
 * BigDecimal on the other hand supports decimal representation and also exponent compression, but only supports
 * integer exponents, which is realistic for real-world values but can quickly become restrictive in scientific
 * setting.
 *
 * (Furthermore, the memory footprint of BigInteger and BigDecimal are much bigger for simpler numbers.)
 */
class Decimal private constructor(

    /**
     * The remaining decimal prefix of this number. It is stored in little-endian order (so element 0 represents
     * bits directly above root).
     */
    private val prefix: IntArray?,

    /**
     * The least significant "root" value of this number. Determines the sign of the number
     * and can be also used to quickly represent small numbers without extra allocation for the prefix array.
     *
     * Furthermore, small modulo operations can be performed without extra pointer lookup, etc.
     */
    private val root: Int,

    /**
     * A decimal exponent. If the exponent overflows, we will throw an Arithmetic exception.
     */
    private val exponent: Long
) {

    companion object {

        val ZERO = Decimal(null, 0, 0)
        val ONE = Decimal(null, 1, 0)
        val TWO = Decimal(null, 2, 0)

        /**
         * Parse given string into a decimal number.
         */
        fun valueOf(input: String): Decimal {
            if (input.isEmpty()) decimalError("empty string")
            val expSplit = input.split('e')
            return when (expSplit.size) {
                1 -> parseValue(input)    // input is just a fractional value with no exponent
                2 -> {
                    val (value, exponent) = expSplit
                    val parsedValue = parseValue(value)
                    // note: toLong will safely reject values which don't fit into long
                    Decimal(parsedValue.prefix, parsedValue.root, exponent.toLong())
                }
                else -> decimalError(input)
            }
        }

        private fun parseValue(value: String): Decimal {
            if (value.isEmpty()) decimalError(value)
            val decimalSplit = value.split('.')
            when (decimalSplit.size) {
                1 -> {
                    // value is just a (possibly long) integer, read sign and parse
                    val sign = when (value[0]) {
                        '-' -> -1; else -> +1
                    }
                    val valueString = when (value[0]) {
                        '-', '+' -> value.drop(1)
                        in '0'..'9' -> value
                        else -> decimalError(value)
                    }
                    val parsed = readIntegral(valueString)
                    Decimal(
                        prefix = parsed.drop(1).takeIf { it.isNotEmpty() }?.toIntArray(),
                        root = parsed.first() * sign,
                        exponent = 0L
                    )
                }
                2 -> {
                    val (decimal, fraction) = decimalSplit
                    // value has an integer and a fraction -> read sign, read a combined number as
                    // integer and apply correct exponent
                    val sign = when (decimal[0]) {
                        '-' -> -1; else -> +1
                    }
                    val valueString = when (decimal[0]) {
                        '-', '+' -> value.drop(1) + fraction
                        in '0'..'9' -> value + fraction
                        else -> decimalError(value)
                    }
                    val parsed = readIntegral(valueString)
                    Decimal(
                        prefix = parsed.drop(1).takeIf { it.isNotEmpty() }?.toIntArray(),
                        root = parsed.first() * sign,
                        exponent = -1L * fraction.length
                    )
                }
            }
            error("")
        }

        /**
         * Read a long integral number into a [Decimal] object.
         *
         * The number can have a leading sign, but otherwise no other characters are allowed (no scientific
         * or decimal notation).
         *
         * Allowed: +123456, -153592, 3, 531, -6
         * Not allowed: 3.14, 1245e-2
         */
        private fun readIntegral(value: String): List<Int> {
            if (value.isEmpty()) decimalError("empty")
            var cursor = 0
            val signum = when (value[0]) {
                '+' -> { cursor += 1; 1 }
                '-' -> { cursor += 1; -1 }
                in '0'..'9' -> 1    // do not increase cursor
                else -> decimalError(value)
            }
            val digitCount = value.length - cursor
            val maxSize = ((4L * digitCount) / 31)  // 4 bits per digit, 31 usable bits per integer
            while (cursor < value.length) {
                val digitChar = value[cursor]
                if (digitChar !in '0'..'9') decimalError(value)
                val digit = digitChar.toInt()
            }
        }

        /* Throw parsing error when Decimal number cannot be read. */
        @Throws(ParseException::class)
        private fun decimalError(value: String): Nothing = throw ParseException("Cannot create decimal representation of \"$value\"")

        /**
         * Represent given integer array as a little-endian integral number. Shift the number by given number
         * of decimal places to the left (increasing the value)
         *
         * value * 10^places
         */
        private fun IntArray.decimalShift(places: Int) {

        }

        private fun IntArray.scalarMultiply(value: Int) {
            var shiftBits = 0
            var multiply = value
            while (multiply % 2 == 0) {
                multiply /= 2; shiftBits += 1
            }
        }

    }


}