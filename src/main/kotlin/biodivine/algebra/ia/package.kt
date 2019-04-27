package biodivine.algebra.ia

import biodivine.algebra.MPoly
import biodivine.algebra.NumQ
import biodivine.algebra.UPoly
import biodivine.algebra.params.SemiAlgSet
import biodivine.algebra.project
import cc.redberry.rings.Rational
import cc.redberry.rings.Rings
import java.awt.Color
import java.awt.image.BufferedImage
import java.io.File
import javax.imageio.ImageIO


fun UPoly.evaluate(x: Interval): Interval {
    var result = Interval(cc(),cc())
    var point = x
    for (i in 1..degree()) {
        result += point * this.get(i)
        point *= x
    }
    return result
}

fun MPoly.evaluate(vararg x: Interval): Interval {
    val powCache = HashMap<Pair<Interval, Int>, Interval>()

    fun pow(x: Interval, e: Int): Interval {
        if (e == 0) error("Unsupported")
        return when (e) {
            1 -> x
            2 -> x*x
            else -> {
                val exp = powCache[x to e]
                if (exp == null) {
                    val a = pow(x, e/2)
                    val computed = if (e % 2 == 0) a*a else a*a*x
                    powCache[x to e] = computed
                    computed
                } else {
                    exp
                }
                /*powCache.computeIfAbsent(x to e) { _ ->
                    val a = pow(x, e/2)
                    if (e % 2 == 0) a*a else a*a*x
                }*/
            }
        }
    }

    var result = Interval(Rings.Q.zero, Rings.Q.zero)
    for (monomial in this) {
        // evaluate the variables
        val vars = monomial.exponents.foldIndexed(Interval(Rings.Q.one, Rings.Q.one)) { i, acc, degree ->
            if (degree == 0) acc else {
                acc * pow(x[i], degree)
            }
        }
        result += vars * monomial.coefficient
    }

    return result
}

fun Rational<MPoly>.evaluate(vars: Array<Interval>, params: Array<Interval>): Interval {
    return this.evaluate(*(params + vars))
}

fun Rational<MPoly>.evaluate(vararg x: Interval): Interval {
    val numerator = this.numerator().evaluate(*x)
    val denominator = this.denominator().evaluate(*x)
    return numerator / denominator
}

fun MPoly.draw(x: Interval, y: Interval, pixelsX: Int, pixelsY: Int,
               image: BufferedImage = BufferedImage(pixelsX, pixelsY, BufferedImage.TYPE_INT_RGB)
): BufferedImage {
    val stepX = x.size / pixelsX
    val stepY = y.size / pixelsY
    var min = Rings.Q.zero
    var max = Rings.Q.zero
    val values: Array<Array<Interval>> = Array(pixelsX) { i ->
        Array(pixelsY) { j ->
            val iX = Interval(x.low + stepX * i, x.low + stepX * (i+1))
            val iY = Interval(y.low + stepY * j, y.low + stepY * (j+1))
            val res = evaluate(iX, iY)
            if (res.low < min) min = res.low
            if (res.high > max) max = res.high
            res
        }
    }
    /*for (i in 0 until pixelsX) {
        for (j in 0 until pixelsY) {
            val res = values[i][j]
            when {
                res.high < Rings.Q.zero -> {
                    val fraction = res.low / min    // should be [0,1]
                    val colorFraction = (fraction * 255)
                    val colorResult = colorFraction.numerator().divide(colorFraction.denominator()).toLong().toInt()
                    image.setRGB(i, j, 0xff0000.and(colorResult.shl(16)))
                }
                res.low > Rings.Q.zero -> {
                    val fraction = res.high / max    // should be [0,1]
                    val colorFraction = (fraction * 255)
                    val colorResult = colorFraction.numerator().divide(colorFraction.denominator()).toLong().toInt()
                    image.setRGB(i, j, 0x00ff00.and(colorResult.shl(8)))
                }
                else -> {
                    image.setRGB(i, j, Color.white.rgb)
                }
            }
        }
    }*/
    for (i in 0 until pixelsX) {
        for (j in 0 until pixelsX) {
            val res = values[i][j]
            if (res.low < Rings.Q.zero && res.high > Rings.Q.zero) {
                image.setRGB(i, j, Color.white.rgb)
            }
        }
    }
    return image
}

fun Rational<MPoly>.draw(x: Interval, y: Interval, pixelsX: Int, pixelsY: Int): BufferedImage {
    val stepX = x.size / pixelsX
    val stepY = y.size / pixelsY
    val image = BufferedImage(pixelsX, pixelsY, BufferedImage.TYPE_INT_RGB);
    var min = Rings.Q.zero
    var max = Rings.Q.zero
    val values: Array<Array<Interval>> = Array(pixelsX) { i ->
        Array(pixelsY) { j ->
            val iX = Interval(x.low + stepX * i, x.low + stepX * (i+1))
            val iY = Interval(y.low + stepY * (pixelsY - j - 1), y.low + stepY * (pixelsY - j))
            val res = evaluate(iX, iY)
            if (res.low < min) min = res.low
            if (res.high > max) max = res.high
            res
        }
    }
    for (i in 0 until pixelsX) {
        for (j in 0 until pixelsY) {
            val res = values[i][j]
            when {
                res.high < Rings.Q.zero -> {
                    val fraction = res.low / min    // should be [0,1]
                    val colorFraction = (fraction * 255)
                    val colorResult = 255 - colorFraction.numerator().divide(colorFraction.denominator()).toLong().toInt()
                    val color = 0xff0000.or(colorResult.shl(8)).or(colorResult)
                    image.setRGB(i, j, color)
                }
                res.low > Rings.Q.zero -> {
                    val fraction = res.high / max    // should be [0,1]
                    val colorFraction = (fraction * 255)
                    val colorResult = 255 - colorFraction.numerator().divide(colorFraction.denominator()).toLong().toInt()
                    val color = 0x00ff00.or(colorResult.shl(16)).or(colorResult)
                    image.setRGB(i, j, color)
                }
                else -> {
                    image.setRGB(i, j, Color.black.rgb)
                }
            }
        }
    }
    return image
}

fun NumQ.mapInto(bounds: Interval, max: Int): Int {
    // get a number in [0..max] which corresponds to the position in bounds interval
    val scale = ((this - bounds.low) / (bounds.high - bounds.low)) * max
    return scale.numerator().divide(scale.denominator()).toLong().toInt()
}

fun SemiAlgSet.draw(x: Interval, y: Interval, pixelsX: Int, pixelsY: Int): BufferedImage {
    val stepX = x.size / pixelsX
    val stepY = y.size / pixelsY
    val image = BufferedImage(pixelsX, pixelsY, BufferedImage.TYPE_INT_RGB);
    this.levelGraph.basis.forEach {
        it.draw(x, y, pixelsX, pixelsY, image)
    }
    this.levelGraph.walkCells().forEach { (point, cell) ->
        println("Cell $cell for $point is ${cell in validCells}")
        if (cell in validCells) {
            image.setRGB(point[0].mapInto(x, pixelsX), point[1].mapInto(y, pixelsY), Color.green.rgb)
        }
    }
    /*for (i in 0 until pixelsX) {
        println("Pixel $i/$pixelsX")
        for (j in 0 until pixelsY) {
            val point = listOf(x.low + stepX*i, y.low + stepY*j)
            val cell = this.levelGraph.cellForPoint(point)
            if (cell in validCells) {
                image.setRGB(i, j, Color.BLACK.rgb)
            }
        }
    }*/
    return image
}

fun main() {
    projectPolyY()
}

fun projectPolyY() {
    //p1 = x
    //p2 = y
    // x = z
    val ring = Rings.MultivariateRingQ(3)   // p1, p2 and x (y = 1)
    val p1 = ring.parse("5/100")
    val p2 = ring.parse("625/10000 * 4/100 * 4/100")
    val p3 = ring.parse("x^2")
    val p4 = ring.parse("1 + x^2")
    val p5 = ring.parse("y")
    val p6 = ring.parse("z + y")
    val p7 = ring.parse("1")
    val p8 = ring.parse("1/10")

    val r1 = Rational(ring, p3.copy(), p4.copy())
    val r2 = Rational(ring, p5.copy(), p6.copy())
    val r3 = Rational(ring, p7.copy(), p4.copy())

    val m1 = p1.copy()
    val m2 = r1.multiply(r2).multiply(p2)
    val m3 = r3.multiply(r2)
    val m4 = p8

    println("M: $m1 | $m2 | $m3 | $m4")
    println("${m2.add(m1).add(m3).subtract(m4)}")
    val b1 = ring.parse("z-1")
    val b2 = ring.parse("z-2")
    val poly = m2.add(m1).add(m3).subtract(m4)
    val num = poly.numerator()
    println("Numerator: $num")
    val proj = project(setOf(num,b1,b2).map { it.asUnivariate(2) }.toSet())

    proj.forEach {
        println("Projection : $it")
    }
}

fun showPolyX() {
    val ring = Rings.MultivariateRingQ(2)
    val p1 = ring.parse("1/10 * x")
    val p2 = ring.parse("1/2")
    val p3 = ring.parse("x + 1/2")
    val p4 = ring.parse("y")
    val p5 = ring.parse("1/2 + y")
    val r1 = Rational(ring, p4, p5)
    val r2 = Rational(ring, p2, p3)
    val polyX = r1.multiply(r2).subtract(p1)
    println("Poly: $polyX")
    val image = polyX.draw(Interval(0, 4), Interval(0, 4), 500, 500)
    ImageIO.write(image, "PNG", File("out.png"))
}