package biodivine.algebra.svg

import biodivine.algebra.NumQ
import biodivine.algebra.ia.decimalString

data class Text(
    val value: String,
    val anchor: Point,
    override val style: Style = Style.FILL
) : SvgPrimitive<Text> {

    override val center: Point = anchor
    override val anchors: Sequence<Point> = sequenceOf(anchor)

    override fun compileSvg(): String = """<text x="${anchor.x.decimalString()}" y="${anchor.y.decimalString()}" font-size="5" ${style.compileAttributes()}>$value</text>"""

    override fun scale(factor: NumQ): Text = this.copy(anchor = anchor * factor)

    override fun translate(delta: Point): Text = this.copy(anchor = anchor + delta)

    override fun flipY(height: NumQ): Text = this.copy(anchor = anchor.flipY(height))
}