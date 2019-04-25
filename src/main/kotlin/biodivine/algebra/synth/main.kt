package biodivine.algebra.synth

import biodivine.algebra.ia.Interval
import cc.redberry.rings.Rational
import cc.redberry.rings.Rings
import java.io.File

fun main() {
    val ring = Rings.MultivariateRingQ(4)
    val parser = ring.mkCoder("x", "y", "p", "q")
    val polyX = run {
        val p1 = parser.parse("1/10 * x")
        val p2 = parser.parse("1/2")
        val p3 = parser.parse("x + 1/2")
        val p4 = parser.parse("y")
        val p5 = parser.parse("1/2 + y")
        val r1 = Rational(ring, p4, p5)
        val r2 = Rational(ring, p2, p3)
        // k_1 * (y)/(K_1 + y) * (K_2)/(K_2 + x) - y_PRB * x
        r1.multiply(r2).subtract(p1)
    }

    val polyY = run {
        val kp = parser.parse("5/100")
        val c1 = parser.parse("1 * 4/100 * 4/100 * 625/10000")
        val num1 = parser.parse("p^2")
        val den1 = parser.parse("y^2 + p^2")
        val num2 = parser.parse("q")
        val den2 = parser.parse("x + q")
        val num3 = parser.parse("y^2")
        val r1 = Rational(ring, num1, den1)
        val r2 = Rational(ring, num2, den2)
        val r3 = Rational(ring, num3, den1)
        val p1 = parser.parse("1/10 * y")
        val m2 = r1.multiply(r2).multiply(c1)
        val m3 = r2.multiply(r3)
        m2.add(m3).add(kp).subtract(p1)
    }
    
    val model = Model(
        ring = ring,
        varNum = 2, varBounds = Box(Interval(Rings.Q.parse("1/10"), Rings.Q.parse("10")), Interval(Rings.Q.parse("1/10"), Rings.Q.parse("10"))),
        paramNum = 2, paramBounds = Box(Interval(3, 5), Interval(4,6)),
        equations = listOf(polyX, polyY)
    )

    val ss = model.computeStateSpace(/*listOf(
        Box(Interval(Rings.Q.parse("1/10"), Rings.Q.parse("10")), Interval(Rings.Q.parse("1/10"), Rings.Q.parse("2"))),
        Box(Interval(Rings.Q.parse("1/10"), Rings.Q.parse("10")), Interval(Rings.Q.parse("2"), Rings.Q.parse("10")))
    )*/listOf(model.varBounds), Rings.Q.parse("1/100000"))
    println("States: ${ss.size}")

    val two = Rings.Q.parse("2")
    val notSmall = ss.indices.filter { i ->
        val b = ss[i]
        b.data[1].high > two
    }.toSet()

    val imageProp = PropertySpace(ss, notSmall).draw().normalize(Rings.Q.parse("1000"))
    File("notSmall.svg").bufferedWriter().use { imageProp.writeTo(it) }

    val image = ss.draw().normalize(Rings.Q.parse("1000"))
    File("ss.svg").bufferedWriter().use { image.writeTo(it) }

    val ts = ss.makeTransitions(model)
    val tsImage = ts.draw().normalize(Rings.Q.parse("1000"))
    File("ts.svg").bufferedWriter().use { tsImage.writeTo(it) }

    val canReachNotSmall = ts.reachBackward(notSmall)
    val imageCanReach = PropertySpace(ss, canReachNotSmall).draw().normalize(Rings.Q.parse("1000"))
    File("canReach.svg").bufferedWriter().use { imageCanReach.writeTo(it) }
}