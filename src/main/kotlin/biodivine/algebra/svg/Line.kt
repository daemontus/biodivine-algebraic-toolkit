package biodivine.algebra.svg

import biodivine.algebra.NumQ
import biodivine.algebra.ia.decimalString

data class Line(
    val begin: Point,
    val end: Point,
    override val style: Style = Style.STROKE
) : SvgPrimitive<Line> {

    override val anchors: Sequence<Point>
        get() = sequenceOf(begin, end)

    override val center: Point
        get() = (begin + end) * p1_2

    override fun compileSvg(): String {
        val (x1, y1) = begin
        val (x2, y2) = end
        return """<line x1="${x1.decimalString()}" y1="${y1.decimalString()}" x2="${x2.decimalString()}" y2="${y2.decimalString()}" ${style.compileAttributes()}/>"""
    }

    override fun scale(factor: NumQ): Line = copy(begin = begin * factor, end = end * factor)

    override fun translate(delta: Point): Line = copy(begin = begin + delta, end = end + delta)

    override fun flipY(height: NumQ): Line = copy(begin = begin.flipY(height), end = end.flipY(height))
}