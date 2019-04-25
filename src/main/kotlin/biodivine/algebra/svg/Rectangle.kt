package biodivine.algebra.svg

import biodivine.algebra.NumQ
import biodivine.algebra.ia.decimalString
import biodivine.algebra.ia.plus
import biodivine.algebra.ia.times
import cc.redberry.rings.Rings

data class Rectangle(
    override val center: Point,
    val dimensions: Point,
    override val style: Style
) : SvgPrimitive<Rectangle> {

    override val anchors: Sequence<Point>
        get() = sequenceOf(leftDown, left, leftUp, up, rightUp, right, rightDown, down)

    val secondaryAnchors: Sequence<Point>
        get() = sequenceOf(leftDown to leftUp, leftUp to rightUp, rightUp to rightDown, rightDown to leftDown).flatMap { (a, b) ->
            sequenceOf((a + (b - a) * p3_10), a + (b - a) * p7_10)
        }

    val sizes: Sequence<NumQ>
        get() = sequenceOf(width, height)

    val leftDown: Point
        get() = center + (dimensions * n1_2)

    val rightUp: Point
        get() = center + (dimensions * p1_2)

    val leftUp: Point
        get() = center + (xy(dimensions.x * n1_2, dimensions.y * p1_2))

    val rightDown: Point
        get() = center + (xy(dimensions.x * p1_2, dimensions.y * n1_2))

    val up: Point
        get() = center + xy(zero, height * p1_2)

    val down: Point
        get() = center + xy(zero, height * n1_2)

    val left: Point
        get() = center + xy(width * n1_2, zero)

    val right: Point
        get() = center + xy(width * p1_2, zero)

    val width: NumQ
        get() = dimensions.x

    val height: NumQ
        get() = dimensions.y

    fun innerPoint(dimension: Dimension, fraction: NumQ): Point {
        return if (dimension == Dimension.X) xy(left.x + fraction * width, center.y)
        else xy(center.x, down.y + fraction * height)
    }

    override fun compileSvg(): String {
        val (x, y) = center + (dimensions * n1_2)
        val (width, height) = dimensions
        return """<rect x="${x.decimalString()}" y="${y.decimalString()}" width="${width.decimalString()}" height="${height.decimalString()}" ${style.compileAttributes()}/>"""
    }

    override fun scale(factor: NumQ): Rectangle = copy(center = center * factor, dimensions = dimensions * factor)

    override fun translate(delta: Point): Rectangle = copy(center = center + delta)

    override fun flipY(height: NumQ): Rectangle = copy(center = center.flipY(height))
}