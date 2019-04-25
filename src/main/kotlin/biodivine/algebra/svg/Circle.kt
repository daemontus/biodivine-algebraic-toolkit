package biodivine.algebra.svg

import biodivine.algebra.NumQ
import biodivine.algebra.ia.decimalString
import biodivine.algebra.ia.times
import biodivine.algebra.ia.unaryMinus

data class Circle(
    override val center: Point,
    val radius: NumQ,
    override val style: Style
) : SvgPrimitive<Circle> {

    override val anchors: Sequence<Point>
        get() = sequenceOf(left, up, right, down)

    val up: Point
        get() = center + Point(zero, radius)

    val down: Point
        get() = center + Point(zero, -radius)

    val left: Point
        get() = center + Point(-radius, zero)

    val right: Point
        get() = center + Point(radius, zero)

    fun anchor(d: Dimension, positive: Boolean): Point = when {
        d == Dimension.X && positive -> right
        d == Dimension.X && !positive -> left
        d == Dimension.Y && positive -> up
        else -> down
    }

    override fun compileSvg(): String =
            """<circle cx="${center.x.decimalString()}" cy="${center.y.decimalString()}" r="${radius.decimalString()}" ${style.compileAttributes()} />"""

    override fun scale(factor: NumQ): Circle =
            copy(center = center * factor, radius = radius * factor)

    override fun flipY(height: NumQ): Circle =
            copy(center = center.flipY(height))

    override fun translate(delta: Point): Circle =
            copy(center = center + delta)
}