package biodivine.algebra.synth

import biodivine.algebra.ia.Interval
import cc.redberry.rings.Rational
import cc.redberry.rings.Rings
import cc.redberry.rings.Rings.Q
import java.io.File
import javax.imageio.ImageIO
import kotlin.system.measureTimeMillis

fun clark(paramIntervals: Map<String, Interval>): Model {
    val constants = mapOf(
        "fast" to Q.parse("100"),
        "Vcell" to Q.parse("52/100000000000000000"),
        "Vcarb" to Q.parse("17/1000000000000000000")
    ).filterKeys { it !in paramIntervals.keys }

    val vars = listOf("fast", "Vcell", "Vcarb", "CO2_cyt", "HCO3_cyt", "CO2_carb", "HCO3_carb")
    val variableIndices = vars.mapIndexed { index, v -> v to index }.toMap()
    val ring = Rings.MultivariateRingQ(vars.size)
    val parser = ring.mkCoder(*vars.toTypedArray())
    fun p(s: String) = parser.parse(s)
    // parameters fast, Vcell, Vcarb
    val eq1 = p("-1*27128*(1/100000000000000)*100000000*CO2_cyt  +  84/100000*HCO3_cyt  +  fast*(9163470612/1000000000000000)  +  -1*fast*CO2_cyt")
    val eq2 = Rational(ring, p("(1/1000000000000000000)*(6807094932/10000000000)"), p("Vcell")).multiply(p("27128*(1/100000000000000)*100000000*CO2_cyt  +  -1*(84/100000)*HCO3_cyt  +  -1*fast*HCO3_cyt  +  fast*HCO3_carb"))
    val eq3 = Rational(ring, p("-1*(19/1000000000000000000000000)*13400"), p("Vcarb")).multiply(Rational(ring, p("CO2_carb"), p("CO2_carb + (2484615385/10000000000000)"))).add(
        p("-1*(1/1000)*fast*27128*(1/100000000000000)*100000000*CO2_carb  +  (119047619/10000)*fast*(84/100000)*HCO3_carb")
    )
    val eq4 = p("fast*HCO3_cyt  +  -1*fast*HCO3_carb  +   (1/1000)*fast*27128*(1/100000000000000)*100000000*CO2_carb  +  -1*(119047619/10000)*fast*(84/100000)*HCO3_carb")

    val constantValuations = constants.toList()
    val constantIndexArray = constantValuations.map { variableIndices.getValue(it.first) }.toIntArray()
    val constantValues = constantValuations.map { it.second }.toTypedArray()
    val reducedRing = Rings.MultivariateRingQ(vars.size - constants.size)
    return Model(
        ring = reducedRing,
        varNum = 4, varBounds = Box(
            Interval(Q.parse("0"), Q.parse("2/1000000")),
            Interval(Q.parse("0"), Q.parse("2/1000000")),
            Interval(Q.parse("0"), Q.parse("2/1000000")),
            Interval(Q.parse("0"), Q.parse("2/1000000"))
        ),
        paramNum = paramIntervals.size, paramBounds = Box(
            *vars.mapNotNull { paramIntervals[it] }.toTypedArray()
        ),
        equations = listOf(
            Rational(reducedRing,
                reducedRing.valueOf(eq1.eliminate(constantIndexArray, constantValues))
            ),
            Rational(reducedRing,
                reducedRing.valueOf(eq2.numerator().eliminate(constantIndexArray, constantValues)),
                reducedRing.valueOf(eq2.denominator().eliminate(constantIndexArray, constantValues))
            ),
            Rational(reducedRing,
                reducedRing.valueOf(eq3.numerator().eliminate(constantIndexArray, constantValues)),
                reducedRing.valueOf(eq3.denominator().eliminate(constantIndexArray, constantValues))
            ),
            Rational(reducedRing,
                reducedRing.valueOf(eq4.eliminate(constantIndexArray, constantValues))
            )
        )
    )
}

fun main() {

    val model = clark(mapOf(
        "fast" to Interval(Q.parse("1"), Q.parse("200")),
        "Vcell" to Interval(Q.parse("26/100000000000000000"), Q.parse("104/100000000000000000"))
        //"fpRB" to Interval(Q.parse("1/1000"), Q.parse("1"))
        //"k1" to Interval(1, 10)
    ))

    val ss = model.computeStateSpace(listOf(model.varBounds), Rings.Q.parse("1/20000"), Rings.Q.parse("1/2000"))
    println("States: ${ss.size}")

    val high = true

    val notSmall = ss.indices.filter { i ->
        val b = ss[i]
        b.data[3].high > Q.parse("1/1000000")
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
            for (s in notSmall) { initial.union(s, solver.one) }

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