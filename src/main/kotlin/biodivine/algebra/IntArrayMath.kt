package biodivine.algebra

import java.math.BigInteger
import java.util.*
import kotlin.math.max

/**
 * IntArrayMath provides unbounded numerical operations for long integer numbers represented as IntArrays.
 *
 * Currently, we assume each integer in such array is non-negative, hence each integer can represent
 * 2^31 unique "symbols". The integers are then ordered from the least significant to the most significant
 * (10467 becomes 76401).
 *
 * Each array can contain any number of trailing zeroes and all operations should be able to support such numbers.
 * This is mostly supported in order to allow pre-allocation and re-use of local arrays which would otherwise had
 * to be discarded.
 *
 * Note that we could utilise all 32 bits of the integer using a little bit-magic, but currently we want
 * to keep it simple, so we want to avoid this approach. The space savings are not that significant, in the
 * future, we should think about moving the operations to native implementation which can be vectorized
 * and unsigned at the same time, which should provide a significant speed-up for large numbers.
 * (but remember that a single native call can take ~200 cycles, so it is only relevant for larger numbers,
 * tens or hundreds of items long, depending on the type of operation)
 *
 * Since we allow trailing zeroes, we use *capacity* to refer to size of the array and *cardinality* to the number
 * of integers which are actually used.
 *
 */
object IntArrayMath {

    private const val radix: Long = 1L.shl(31)
    private val radixI = BigInteger(radix.toString())

    private const val itemWidth: Int = 31
    private const val itemMask: Long = 0b1111111111111111111111111111111
    private const val itemDomain: Long = 1L.shl(31)

    /**
     * c = a + b
     *
     * We perform standard algorithm, adding each item and transferring any overflow using the carry variable.
     *
     * @param a First addition argument
     * @param b Second addition argument
     * @param c Output array, capacity at least max(a, b) + 1
     * @throws [ArithmeticException] If the result exceeds the capacity of the output array.
     */
    @Throws(ArithmeticException::class)
    fun addition(a: IntArray, b: IntArray, c: IntArray) {
        if (a.size < b.size) return addition(b, a, c)
        try {
            // |A| >= |B|
            var carry = 0L
            for (i in 0 until b.size) {             // add common items
                val sum = a[i].toLong() + b[i].toLong() + carry
                c[i] = sum.and(itemMask).toInt()
                carry = sum.shr(itemWidth)
            }
            for (i in b.size until a.size) {        // add extra items in A, preserving carry
                val sum = a[i].toLong() + carry
                c[i] = sum.and(itemMask).toInt()
                carry = sum.shr(itemWidth)
            }
            for (i in a.size until c.size) {        // zero remaining items, plus apply carry
                c[i] = carry.toInt()
                carry = 0
            }
        } catch (e: ArrayIndexOutOfBoundsException) {
            throw ArithmeticException("Addition: output array too small. |A| = ${a.size}, |B| = ${b.size}, |C| = ${c.size}")
        }
    }

    /**
     * c = a - b without underflow (assuming a >= b)
     *
     * Subtraction is very similar to addition, except the carry is handled a bit differently and we have to
     * account for the fact that signed integers underflow to negative numbers.
     *
     * @param a Larger subtraction argument
     * @param b Smaller subtraction argument
     * @param c Output array, capacity at least the same as a
     */
    fun subtraction(a: IntArray, b: IntArray, c: IntArray) {
        try {
            var i = 0
            var carry = false
            while (i < a.size || i < b.size) {
                val base = a.getOr(i, 0).toLong()
                val remove = b.getOr(i, 0).toLong() + if (carry) 1 else 0
                val sub = base - remove
                c[i] = ((sub + radix) % radix).toInt()
                carry = sub < 0 // underflow
                i += 1
            }
            if (carry) {
                error("A is smaller than B in subtraction! A: ${Arrays.toString(a)}, B: ${Arrays.toString(b)}")
            }
            while (i < c.size) {
                c[i] = 0
                i += 1
            }
        } catch (e: ArrayIndexOutOfBoundsException) {
            throw ArithmeticException("Subtraction: output array too small. |A| = ${a.size}, |B| = ${b.size}, |C| = ${c.size}")
        }
    }

    /**
     * c = a * b
     *
     * Currently, we only implement standard "school" algorithm for multiplication which is O(n^2), since
     * it first multiplies the larger number by each item in the smaller one, shifting it by its position.
     * The result is then the sum of these numbers.
     *
     * Capacity of c must be at least capacity of a plus capacity of b.
     */
    fun multiplication(a: IntArray, b: IntArray, c: IntArray) {
        if (a.size < b.size) return multiplication(b, a, c)
        // C starts ar zero
        for (i in c.indices) c[i] = 0
        // A is longer
        for (i in b.indices) {
            if (b[i] != 0) {
                constMulShiftAdd(a, b[i], i, c)
            }
        }
    }

    /**
     * Fused, in-place constant multiplication + shifted addition. This operation is used during basic multiplication
     * and basically amounts to the following:
     *
     * c += ((a * b) << shift)
     *
     * (Note that out is modified, not overwritten)
     *
     * Can throw [ArrayIndexOutOfBoundsException] when [c] is not big enough!
     */
    private fun constMulShiftAdd(a: IntArray, b: Int, shift: Int, c: IntArray) {
        var i = 0   // index into the A array, use (i+shift) for C indices
        var carry: Long = 0
        while (i < a.size) {
            // Note that the additions won't overflow long because even if we assume maximal values,
            // two extra additions will fit (also, we have one extra unused bit)
            // 1111 * 1111 = 1110 0001
            // 1110 0001 + 1111 + 1111 = 1111 1111
            val sum = c[i+shift].toLong() + a[i].toLong() * b.toLong() + carry
            c[i+shift] = sum.and(itemMask).toInt()
            carry = sum.shr(itemWidth)
            i += 1
        }
        // Now we have added most of A, but we still need to apply carry!
        while (carry > 0) {
            val sum = c[i+shift].toLong() + carry
            c[i+shift] = sum.and(itemMask).toInt()
            carry = sum.shr(itemWidth)
            i += 1
        }
        // Now the carry is also applied and we leave the rest of c alone.
    }

    /**
     * Standard equality test.
     */
    fun equals(a: IntArray, b: IntArray): Boolean {
        if (a.size < b.size) return equals(b, a)
        var i = a.size - 1
        while (i >= 0) {
            if (a[i] != b.getOr(i, 0)) {
                return false
            }
            i -= 1
        }
        return true
    }

    /**
     * Standard comparison operation ignoring trailing zeroes. Return -1/0/1 such that:
     * a (op) b <=> compare(a, b) (op) 0
     *  - 0 means the numbers are equal
     *  - 1 means a is strictly larger (a > b <=> compare(a, b) > 0)
     *  - -1 means a is strictly smaller (a < b <=> compare(a, b) < 0)
     */
    fun compare(a: IntArray, b: IntArray): Int {
        var i = max(a.size, b.size) - 1
        while (i >= 0) {
            val x = a.getOr(i, 0)
            val y = b.getOr(i, 0)
            if (x > y) return 1
            if (x < y) return -1
            i -= 1
        }
        return 0
    }

    /**
     * Actual number of used items in this array.
     */
    fun cardinality(x: IntArray): Int {
        var c = x.size
        while (c-1 >= 0 && x[c-1] == 0) c -= 1
        return c
    }

    private fun IntArray.getOr(index: Int, default: Int) = if (index < size) get(index) else default

    /**
     * Convert our int array number representation to standard big integer.
     */
    fun toBigInteger(a: IntArray): BigInteger {
        var result = BigInteger.ZERO
        var radix = BigInteger.ONE
        for (i in a) {
            result += BigInteger(i.toString()) * radix
            radix *= radixI
        }
        return result
    }

}