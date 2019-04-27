package biodivine.algebra.synth

import biodivine.algebra.ia.Interval
import biodivine.algebra.ia.draw
import cc.redberry.rings.Rational
import cc.redberry.rings.Rings
import java.io.File
import javax.imageio.ImageIO
import kotlin.system.measureTimeMillis

val two = Rings.Q.parse("2")
val ten = Rings.Q.parse("10")

fun main() {
    val ring = Rings.MultivariateRingQ(4)
    val parser = ring.mkCoder("p", "q", "x", "y")
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

    println("PolyY: ${polyY}")
    
    val model = Model(
        ring = ring,
        varNum = 2, varBounds = Box(Interval(Rings.Q.parse("1/2"), Rings.Q.parse("10")), Interval(Rings.Q.parse("1/2"), Rings.Q.parse("10"))),
        paramNum = 2, paramBounds = Box(Interval(3, 5), Interval(4,6)),
        equations = listOf(polyX, polyY)
    )

    val ss = model.computeStateSpace(/*listOf(
        Box(Interval(Rings.Q.parse("0"), Rings.Q.parse("10")), Interval(Rings.Q.parse("0"), Rings.Q.parse("2"))),
        Box(Interval(Rings.Q.parse("0"), Rings.Q.parse("10")), Interval(Rings.Q.parse("2"), Rings.Q.parse("10")))
    )*/listOf(model.varBounds), Rings.Q.parse("1/10000"), Rings.Q.parse("1/1000"))
    println("States: ${ss.size}")

    val notSmall = ss.indices.filter { i ->
        val b = ss[i]
        //b.data[1].high > two
        b.data[1].low < two
    }.toSet()

    /*val imageProp = PropertySpace(ss, notSmall).draw().normalize(Rings.Q.parse("1000"))
    File("notSmall.svg").bufferedWriter().use { imageProp.writeTo(it) }

    val image = ss.draw().normalize(Rings.Q.parse("1000"))
    File("ss.svg").bufferedWriter().use { image.writeTo(it) }

    val ts = ss.makeTransitions(model)
    val tsImage = ts.draw().normalize(Rings.Q.parse("1000"))
    File("ts.svg").bufferedWriter().use { tsImage.writeTo(it) }

    val canReachNotSmall = ts.reachBackward(notSmall)
    val imageCanReach = PropertySpace(ss, canReachNotSmall).draw().normalize(Rings.Q.parse("1000"))
    File("canReach.svg").bufferedWriter().use { imageCanReach.writeTo(it) }*/

    val elapsed = measureTimeMillis {
        val ts = model.makeSemiAlgTransitions(ss)
        ts.run {
            val initial = ConcurrentArrayStateMap(ts.states.size, solver)
            for (s in notSmall) { initial.union(s, solver.one) }

            // AG small = ! EF ! small
            val efNotSmall = initial.reachBackward()

            println("Reachability done.")

            solver.run {
                var notAG = zero
                for (s in 0 until ts.states.size) {
                    println("$s / ${ts.states.size}")
                    //println("In state ${ts.states[s]} can reach large ${efNotSmall.get(s)}")
                    notAG = notAG or efNotSmall.get(s).not()
                }
                println("Not AG: $notAG")

                val image = notAG.draw(model.paramBounds.data[0], model.paramBounds.data[1], 500, 500)
                ImageIO.write(image, "PNG", File("out.png"))
            }
        }
    }

    println("Making transitions: $elapsed")
}