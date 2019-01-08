package biodivine.algebra

import org.junit.Test
import java.math.BigInteger
import java.util.*
import kotlin.math.absoluteValue
import kotlin.random.Random
import kotlin.test.assertEquals
import kotlin.test.assertFails
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class IntArrayMathTest {

    // The random number generator is pretty expensive, so we make N numbers and test
    // addition on every combination, which gives us N*N combinations.
    private val n = 100
    private val k = 50
    private val random = Random(42)
    private val numbers = Array(n) { randomNumber(random, k) }

    @Test
    fun test_toBigInteger_01_zero() {
        assertEquals(BigInteger.ZERO, IntArrayMath.toBigInteger(intArrayOf()))
    }

    @Test
    fun test_toBigInteger_02_single() {
        val a = intArrayOf(0b110101)
        val aI = BigInteger("110101", 2)
        assertEquals(aI, IntArrayMath.toBigInteger(a))
    }

    @Test
    fun test_toBigInteger_03_small() {
        val a = intArrayOf(
            0b1101010111011000100100101001010,
            0b1011010100010111001011101001010
        )

        val aI = BigInteger("1011010100010111001011101001010" + "1101010111011000100100101001010", 2)

        assertEquals(aI, IntArrayMath.toBigInteger(a))
    }

    @Test
    fun test_toBigInteger_04_medium() {
        val a = intArrayOf(
            0b1111110010111101011110001011001,
            0b1110100010101010110101010100010,
            0b0000001011001011111001001000001,
            0b1111110110101011111000010010001
        )

        val aI = BigInteger(
            "1111110110101011111000010010001"
                    + "0000001011001011111001001000001"
                    + "1110100010101010110101010100010"
                    + "1111110010111101011110001011001",
            2
        )

        assertEquals(aI, IntArrayMath.toBigInteger(a))
    }

    @Test
    fun test_toBigInteger_05_large() {

        val a = intArrayOf(
            0b1101001011010101010100010101010,
            0b1010101101110100000101100001010,
            0b0000110111111010101010101111100,
            0b0000000000000000000000000000000,
            0b1111111111111111111111111111111,
            0b0000001111101010111101111010001
        )

        val aI = BigInteger(
            "0000001111101010111101111010001"
                    + "1111111111111111111111111111111"
                    + "0000000000000000000000000000000"
                    + "0000110111111010101010101111100"
                    + "1010101101110100000101100001010"
                    + "1101001011010101010100010101010",
            2
        )

        assertEquals(aI, IntArrayMath.toBigInteger(a))
    }

    @Test
    fun test_equals_01_small() {
        val a = intArrayOf(45)
        val b = intArrayOf(45, 0, 0)
        val c = intArrayOf(0, 45, 0)
        val d = intArrayOf(33)
        assertTrue { IntArrayMath.equals(a, a) }
        assertTrue { IntArrayMath.equals(a, b) }
        assertTrue { IntArrayMath.equals(b, a) }
        assertFalse { IntArrayMath.equals(a, c) }
        assertFalse { IntArrayMath.equals(c, a) }
        assertFalse { IntArrayMath.equals(a, d) }
    }

    @Test
    fun test_equals_02_fuzz() {
        for (a in numbers) {
            for (b in numbers) {
                assertEquals(
                    IntArrayMath.toBigInteger(a) == IntArrayMath.toBigInteger(b),
                    IntArrayMath.equals(a, b)
                )
            }
        }
    }

    @Test
    fun test_compare_01_small() {
        val a = intArrayOf(45)
        val b = intArrayOf(5673, 0, 0)
        val c = intArrayOf(4, 0)
        val d = intArrayOf(5673)

        // a < b
        assertTrue { IntArrayMath.compare(a, b) < 0 }
        // b > a
        assertTrue { IntArrayMath.compare(b, a) > 0 }
        // c < a
        assertTrue { IntArrayMath.compare(c, a) < 0 }
        // b >= d
        assertTrue { IntArrayMath.compare(b, d) >= 0 }
        // b > c
        assertTrue { IntArrayMath.compare(b, c) > 0 }
        // b = d
        assertTrue { IntArrayMath.compare(b, d) == 0 }
        // a < d
        assertTrue { IntArrayMath.compare(a, d) < 0 }
    }

    @Test
    fun test_compare_02_fuzz() {
        val n = 100; val k = 50
        val random = Random(42)
        val numbers = Array(n) { randomNumber(random, k) }
        for (a in numbers) {
            for (b in numbers) {
                assertEquals(
                    IntArrayMath.toBigInteger(a).compareTo(IntArrayMath.toBigInteger(b)),
                    IntArrayMath.compare(a, b)
                )
            }
        }
    }

    @Test
    fun test_addition_01_zero() {
        val a = intArrayOf(36)
        val zero = intArrayOf()
        val expected = intArrayOf(36, 0)
        val c = IntArray(2)

        // c = a + zero
        c.reset()
        IntArrayMath.addition(a, zero, c)
        assertEquals(IntArrayMath.toBigInteger(expected), IntArrayMath.toBigInteger(c))

        // c = zero + a
        c.reset()
        IntArrayMath.addition(zero, a, c)
        assertEquals(IntArrayMath.toBigInteger(expected), IntArrayMath.toBigInteger(c))
    }

    @Test
    fun test_addition_01_single() {
        val a = intArrayOf(12345)
        val b = intArrayOf(678)
        val c = IntArray(2)

        val expected = IntArrayMath.toBigInteger(a) + IntArrayMath.toBigInteger(b)

        c.reset()
        IntArrayMath.addition(a, b, c)
        assertEquals(expected, IntArrayMath.toBigInteger(c))

        c.reset()
        IntArrayMath.addition(b, a, c)
        assertEquals(expected, IntArrayMath.toBigInteger(c))
    }

    @Test
    fun test_addition_02_small() {
        val a = intArrayOf(
            0b1111111111111111111111111111111,
            0b1110101101011101110011011010101
        )
        val b = intArrayOf(
            0b1110110101010110101000100101110,
            0b1111100101010101010100100010001
        )
        val c = IntArray(3)

        val expected = IntArrayMath.toBigInteger(a) + IntArrayMath.toBigInteger(b)

        c.reset()
        IntArrayMath.addition(a, b, c)
        assertEquals(expected, IntArrayMath.toBigInteger(c))

        c.reset()
        IntArrayMath.addition(b, a, c)
        assertEquals(expected, IntArrayMath.toBigInteger(c))
    }

    @Test
    fun test_addition_03_large() {
        val a = intArrayOf(
            0b1101010101101010101010101010001,
            0b1111011111101101000000100101001,
            0b0000000000000000000000000000000,
            0b1111111111111111111111111111111,
            0b1110101101011101110011011010101
        )
        val b = intArrayOf(
            0b1110110101010110101000100101110,
            0b0101010101010101011111101110100,
            0b0000000000000000000000000000000,
            0b0010000101011011000001100101010,
            0b1111111111111111111111111111111,
            0b1111100101010101010100100010001
        )
        val c = IntArray(7)

        val expected = IntArrayMath.toBigInteger(a) + IntArrayMath.toBigInteger(b)

        c.reset()
        IntArrayMath.addition(a, b, c)
        assertEquals(expected, IntArrayMath.toBigInteger(c))

        c.reset()
        IntArrayMath.addition(b, a, c)
        assertEquals(expected, IntArrayMath.toBigInteger(c))
    }

    @Test
    fun test_addition_04_fuzz() {
        // The test is still pretty slow, but this is due to BigInteger conversion, the addition itself
        // is only ~10% of the runtime.
        val c = IntArray(k + 1)
        for (a in numbers) {
            for (b in numbers) {
                c.reset()
                IntArrayMath.addition(a, b, c)
                val expected = IntArrayMath.toBigInteger(a) + IntArrayMath.toBigInteger(b)
                assertEquals(expected, IntArrayMath.toBigInteger(c))
            }
        }
    }

    @Test
    fun test_subtraction_01_small() {
        val a = intArrayOf(45, 0)
        val b = intArrayOf(63)
        val zero = intArrayOf(0)
        val c = IntArray(4)

        c.reset()
        IntArrayMath.subtraction(b, a, c)
        assertTrue { IntArrayMath.equals(intArrayOf(18), c) }
        assertFails { IntArrayMath.subtraction(a, b, c) }

        c.reset()
        IntArrayMath.subtraction(a, zero, c)
        assertTrue { IntArrayMath.equals(a, c) }
        assertFails { IntArrayMath.subtraction(zero, a, c) }
    }

    @Test
    fun test_subtraction_02_fuzz() {
        val c = IntArray(k)
        for (a in numbers) {
            for (b in numbers) {
                val aI = IntArrayMath.toBigInteger(a)
                val bI = IntArrayMath.toBigInteger(b)
                if (IntArrayMath.compare(a, b) >= 0) {   // a >= b
                    c.reset()
                    IntArrayMath.subtraction(a, b, c)
                    val cI = IntArrayMath.toBigInteger(c)
                    assertEquals(aI - bI, cI,
                        "Subtraction of ${Arrays.toString(a)} - ${Arrays.toString(b)}"
                    )
                } else {
                    assertFails("Subtraction of $aI - $bI should fail!") {
                        IntArrayMath.subtraction(a, b, c)
                    }
                }
            }
        }
    }

    @Test
    fun test_multiplication_01_small() {
        val a = intArrayOf(35, 0)
        val b = intArrayOf(45)
        val zero = intArrayOf(0)
        val c = IntArray(3)

        c.reset()
        IntArrayMath.multiplication(a, b, c)
        assertTrue { IntArrayMath.equals(intArrayOf(1575), c) }

        c.reset()
        IntArrayMath.multiplication(b, a, c)
        assertTrue { IntArrayMath.equals(intArrayOf(1575), c) }

        c.reset()
        IntArrayMath.multiplication(a, zero, c)
        assertTrue { IntArrayMath.equals(zero, c) }
    }

    @Test
    fun test_multiplication_02_fuzz() {
        val c = IntArray(k + k)
        for (a in numbers) {
            for (b in numbers) {
                val aI = IntArrayMath.toBigInteger(a)
                val bI = IntArrayMath.toBigInteger(b)
                c.reset()
                IntArrayMath.multiplication(a, b, c)
                val cI = IntArrayMath.toBigInteger(c)
                assertEquals(aI * bI, cI)
            }
        }
    }

    private fun randomNumber(random: Random, maxLength: Int): IntArray {
        val length = (random.nextInt() % maxLength).absoluteValue
        return IntArray(length) { random.nextInt().absoluteValue }
    }

    /* Fill array with some random garbage. */
    private fun IntArray.reset() {
        for (i in indices) {
            this[i] = i
        }
    }

}