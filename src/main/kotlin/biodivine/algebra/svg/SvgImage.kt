package biodivine.algebra.svg

import biodivine.algebra.NumQ
import biodivine.algebra.ia.decimalString
import biodivine.algebra.ia.div
import biodivine.algebra.ia.times
import java.io.Writer

data class SvgImage(
    val primitives: List<SvgPrimitive<*>>,
    val arrowSize: NumQ
) {

    private val dimensions: Point
    private val down: Point
    private val up: Point

    init {
        var minX = zero ; var minY = zero
        var maxX = zero ; var maxY = zero
        primitives.asSequence().flatMap { it.anchors }.forEach { (x,y) ->
            if (x > maxX) maxX = x
            if (x < minX) minX = x
            if (y > maxY) maxY = y
            if (y < minY) minY = y
        }
        up = xy(maxX, maxY)
        down = xy(minX, minY)
        dimensions = up - down
    }

    /**
     * Scale and offset this SVG image to target width, preserving the aspect ratio and making everything visible.
     */
    fun normalize(targetWidth: NumQ? = null): SvgImage {
        val delta = down * n1
        val factor = if (targetWidth != null) targetWidth / dimensions.y else p1
        // the order of transformations is important!
        return copy(primitives = primitives.map { it.translate(delta).flipY(dimensions.y).scale(factor) }, arrowSize = factor * arrowSize)
    }

    fun writeTo(writer: Writer) = writer.run {
        appendln(svgPrefix())
        primitives.forEach {
            appendln(it.compileSvg())
        }
        appendln(svgSuffix())
    }

    fun compileSvg(): String = buildString {
        appendln(svgPrefix())
        primitives.forEach {
            appendln(it.compileSvg())
        }
        appendln(svgSuffix())
    }

    private fun svgPrefix(): String =
"""<?xml version="1.0" standalone="no"?>
<svg width="${dimensions.x}" height="${dimensions.y}" version="1.1" xmlns="http://www.w3.org/2000/svg">
    <defs>
        <marker id="arrow" markerWidth="${arrowSize.decimalString()}" markerHeight="${arrowSize.decimalString()}" refX="${arrowSize.decimalString()}" refY="${(arrowSize/2).decimalString()}" orient="auto" markerUnits="strokeWidth">
            <path d="M0,0 L0,${arrowSize.decimalString()} L${arrowSize.decimalString()},${(arrowSize/2).decimalString()} z" fill="#000" />
        </marker>
    </defs>
"""

    private fun svgSuffix(): String =
"""
</svg>
"""

}