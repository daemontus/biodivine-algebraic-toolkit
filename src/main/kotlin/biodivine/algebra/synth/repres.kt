package biodivine.algebra.synth

import biodivine.algebra.NumQ
import cc.redberry.rings.Rational
import cc.redberry.rings.Rings
import biodivine.algebra.ia.Interval
import cc.redberry.rings.Rings.Q
import java.io.File
import javax.imageio.ImageIO
import kotlin.system.measureTimeMillis

fun represilator(paramIntervals: Map<String, Interval>): Model {

    val constants: Map<String, NumQ> = mapOf(
        "k1" to Q.one,
        "k2" to Q.one,
        "K1" to Q.parse("5"),
        "K2" to Q.parse("5"),
        "f1" to Q.parse("1/10"),
        "f2" to Q.parse("1/10")
    ).filterKeys { it !in paramIntervals.keys }

    val vars = listOf("k1", "k2", "K1", "K2", "f1", "f2", "x1", "x2")
    val ring = Rings.MultivariateRingQ(vars.size)
    val variableIndices = vars.mapIndexed { index, v -> v to index }.toMap()
    val parser = ring.mkCoder(*vars.toTypedArray())
    fun p(s: String) = parser.parse(s)
    val r1 = Rational(ring, p("K1^5"), p("K1^5 + x2^5"))
    val r2 = Rational(ring, p("K2^5"), p("K2^5 + x1^5"))

    val eq1 = r1.multiply(p("k1")).subtract(p("f1 * x1"))
    val eq2 = r2.multiply(p("k2")).subtract(p("f2 * x2"))

    val constantValuations = constants.toList()
    val constantIndexArray = constantValuations.map { variableIndices.getValue(it.first) }.toIntArray()
    val constantValues = constantValuations.map { it.second }.toTypedArray()
    val reducedRing = Rings.MultivariateRingQ(vars.size - constants.size)
    return Model(
        ring = reducedRing,
        varNum = 2, varBounds = Box(
            Interval(Q.parse("0"), Q.parse("15")),
            Interval(Q.parse("0"), Q.parse("15"))
        ),
        paramNum = paramIntervals.size, paramBounds = Box(
            *vars.mapNotNull { paramIntervals[it] }.toTypedArray()
        ),
        equations = listOf(
            Rational(reducedRing,
                reducedRing.valueOf(eq1.numerator().eliminate(constantIndexArray, constantValues)),
                reducedRing.valueOf(eq1.denominator().eliminate(constantIndexArray, constantValues))
            ),
            Rational(reducedRing,
                reducedRing.valueOf(eq2.numerator().eliminate(constantIndexArray, constantValues)),
                reducedRing.valueOf(eq2.denominator().eliminate(constantIndexArray, constantValues))
            )
        )
    )
}

fun main() {

    val model = represilator(mapOf(
        "K1" to Interval(Q.parse("1/100"), Q.parse("10")),
        "f1" to Interval(Q.parse("1/100"), Q.parse("1"))
    ))

    println("Model: $model")

    val ss = model.computeStateSpace(listOf(model.varBounds), Rings.Q.parse("1/100000"), Rings.Q.parse("1/10000"))
    println("States: ${ss.size}")

    val high = false
    val four = Q.parse("4")

    val notSmall = ss.indices.filter { i ->
        val b = ss[i]
        if (high) b.data[1].low < four else b.data[1].high > four
    }.toSet()

    println("Not small ${notSmall.size}")

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
            for (s in notSmall) {
                /*successors[s]!!.forEach { (t, p) ->
                    if (t !in notSmall) println("Outgoing params: $p")
                }*/
                initial.union(s, solver.one)
            }

            // AG small = ! EF ! small
            val efNotSmall = initial.reachBackward()

            println("Reachability done.")

            solver.run {
                var notAG = zero
                for (s in 0 until ts.states.size) {
                    if (s%1000 == 0) println("$s / ${ts.states.size}")
                    notAG = notAG or (efNotSmall.get(s).not())
                }
                println("Not AG: $notAG")

                val image = notAG.draw(500, 500)
                ImageIO.write(image, "PNG", File("AG_${if (high) "high" else "low"}.png"))
            }
        }
    }

    println("Synthesis: $elapsed")
}