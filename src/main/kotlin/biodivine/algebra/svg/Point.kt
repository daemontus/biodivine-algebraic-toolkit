package biodivine.algebra.svg

import biodivine.algebra.NumQ
import biodivine.algebra.ia.minus
import biodivine.algebra.ia.plus
import biodivine.algebra.ia.times
import cc.redberry.rings.Rings

/**
 * Point in 2D space. We assume a 2D standard vector space with addition and multiplication by scalar.
 */
data class Point(val x: NumQ, val y: NumQ) : Comparable<Point> {

    companion object {
        val ZERO = Point(Rings.Q.zero, Rings.Q.zero)
    }

    operator fun plus(other: Point) = Point(x + other.x, y + other.y)

    operator fun minus(other: Point) = Point(x - other.x, y - other.y)

    operator fun times(num: NumQ) = Point(x * num, y * num)

    fun flipY(height: NumQ) = Point(x, height - y)

    override fun compareTo(other: Point): Int = if (x == other.x) y.compareTo(other.y) else x.compareTo(other.x)
}