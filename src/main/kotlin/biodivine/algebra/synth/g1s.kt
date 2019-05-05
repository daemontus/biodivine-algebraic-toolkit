package biodivine.algebra.synth

import biodivine.algebra.NumQ
import cc.redberry.rings.Rational
import cc.redberry.rings.Rings
import biodivine.algebra.ia.Interval
import cc.redberry.rings.Rings.Q
import java.io.File
import javax.imageio.ImageIO
import kotlin.system.measureTimeMillis

fun G1S(paramIntervals: Map<String, Interval>): Model {

    val constants: Map<String, NumQ> = mapOf(
        "a" to Q.parse("4/100"),
        "k1" to Q.parse("1"),
        "k2" to Q.parse("16/10"),
        //"k2" to Q.parse("1"),
        "kp" to Q.parse("5/100"),
        "fpRB" to Q.parse("5/1000"),
        //"fpRB" to Q.parse("1/10"),
        "fE2F1" to Q.parse("1/10"),
        "J11" to Q.parse("1/2"),
        "J12" to Q.parse("5"),
        "Km1" to Q.parse("1/2"),
        "Km2" to Q.parse("4")
    ).filterKeys { it !in paramIntervals.keys }

    val ring = Rings.MultivariateRingQ(12)
    val vars = listOf("k1", "Km1", "J11", "fpRB", "kp", "k2", "a", "Km2", "J12", "fE2F1", "pRB", "E2F1")
    val variableIndices = vars.mapIndexed { index, v -> v to index }.toMap()
    val parser = ring.mkCoder(*vars.toTypedArray())
    fun p(s: String) = parser.parse(s)
    val r1 = Rational(ring, p("E2F1"), p("Km1 + E2F1"))
    val r2 = Rational(ring, p("J11"), p("J11 + pRB"))
    val r3 = Rational(ring, p("a^2 + E2F1^2"), p("Km2^2 + E2F1^2"))
    val r4 = Rational(ring, p("J12"), p("J12 + pRB"))

    val eq1 = r1.multiply(r2).multiply(p("k1")).subtract(p("fpRB * pRB"))
    val eq2 = r3.multiply(r4).multiply(p("k2")).subtract(p("fE2F1 * E2F1")).add(p("kp"))

    println(parser.encode(eq1.numerator()))
    println(parser.encode(eq2.numerator()))
    val constantValuations = constants.toList()
    val constantIndexArray = constantValuations.map { variableIndices.getValue(it.first) }.toIntArray()
    val constantValues = constantValuations.map { it.second }.toTypedArray()
    val reducedRing = Rings.MultivariateRingQ(12 - constants.size)
    return Model(
        ring = reducedRing,
        varNum = 2, varBounds = Box(
            Interval(Q.parse("0"), Q.parse("5")),
            Interval(Q.parse("0"), Q.parse("5"))
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

    val model = G1S(mapOf(
        //"Km1" to Interval(Q.parse("1/100"), Q.one),
        //"J11" to Interval(Q.parse("1/100"), Q.one),
        "Km2" to Interval(Q.parse("1/100"), Q.parse("6")),    // 3;5
        "J12" to Interval(Q.parse("1/100"), Q.parse("6"))     // 4;6
        //"fpRB" to Interval(Q.parse("1/1000"), Q.parse("1"))
        //"k1" to Interval(1, 10)
    ))

    val ss = model.computeStateSpace(listOf(model.varBounds), Rings.Q.parse("1/40000"), Rings.Q.parse("1/8000"))
    println("States: ${ss.size}")

    val high = true

    val notSmall = ss.indices.filter { i ->
        val b = ss[i]
        if (high) b.data[1].low < two else b.data[1].high > two
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
                println("One is: ${solver.one.isEmpty()}")
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