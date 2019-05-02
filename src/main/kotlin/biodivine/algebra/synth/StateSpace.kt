package biodivine.algebra.synth

import biodivine.algebra.MPoly
import biodivine.algebra.NumQ
import biodivine.algebra.ia.div
import biodivine.algebra.ia.evaluate
import biodivine.algebra.ia.plus
import biodivine.algebra.ia.times
import biodivine.algebra.params.LevelGraph
import biodivine.algebra.params.SemiAlgSet
import biodivine.algebra.params.SemiAlgSolver
import biodivine.algebra.params.SemiAlgTreeSolver
import biodivine.algebra.svg.*
import cc.redberry.rings.Rings
import cc.redberry.rings.Rings.Q
import kotlin.system.exitProcess

fun Model.computeStateSpace(initial: List<Box>, maxDivision: NumQ, minDivision: NumQ): List<Box> {
    val result = ArrayList<Box>()
    val bounds = this.varBounds
    val smallest = bounds.volume * maxDivision
    val largest = bounds.volume * minDivision
    var toResolve = initial
    var safeVolume = Rings.Q.zero
    var unsafeVolume = Rings.Q.zero
    while (toResolve.isNotEmpty()) {
        val nextResolve = ArrayList<Box>()
        for (box in toResolve) {
            when {
                box.volume > largest -> {
                    nextResolve.addAll(box.subdivide())
                }
                box.volume < smallest -> {
                    //println("Unsafe: $box")
                    result.add(box)
                    unsafeVolume += box.volume
                } // the box is too small, tough luck
                equations.any { !it.evaluate(box.data, paramBounds.data).hasZero } -> {
                    //println("Is safe $box: ${equations.map { it.evaluate(box.data, paramBounds.data) }}")
                    result.add(box)
                    safeVolume += box.volume
                } // some of the variables is non-zero, no equilibrium in here!
                else -> {
                    nextResolve.addAll(box.subdivide())
                }
            }
        }
        toResolve = nextResolve
        println("To resolve: ${toResolve.size}")
    }

    println("Safe: $safeVolume, unsafe: $unsafeVolume. Safe ratio: ${safeVolume / bounds.volume}")
    return result
}

fun List<Box>.makeTransitions(model: Model): RectTransitionSystem {
    val dimensions = this.first().data.size
    val transitions = ArrayList<Pair<Int, Int>>()
    val boxToIndex = this.mapIndexed { i, box -> box to i }.toMap()
    for (state in this) {
        neighbour@ for (other in this) {
            if (state === other) continue@neighbour
            var transitionDimension = -1
            for (d in 0 until dimensions) {
                if (!state.data[d].intersects(other.data[d])) {
                    // if they do not intersect, they are not related!
                    continue@neighbour
                }
                if (state.data[d].low == other.data[d].high || state.data[d].high == other.data[d].low) {
                    if (transitionDimension > -1) {
                        // these two boxes are already touching in other dimension, they can't touch in multiple
                        continue@neighbour
                    } else {
                        transitionDimension = d
                    }
                }
            }
            // now check if the transition even can be enabled under some circumstances:
            val vars = Array(dimensions) { d ->
                (state.data[d] intersect other.data[d])!!
            }
            if (vars[transitionDimension].low != vars[transitionDimension].high) error("WTF?")
            val possibleDerivatives = model.equations[transitionDimension].evaluate(vars, model.paramBounds.data)
            if (vars[transitionDimension].low == state.data[transitionDimension].high) {
                // we are going up from the state, so the derivative should be positive
                if (possibleDerivatives.high > zero) {
                    transitions.add(boxToIndex.getValue(state) to boxToIndex.getValue(other))
                }
            } else {
                // we are going down from the state, so the derivative should be negative
                if (possibleDerivatives.low < zero) {
                    transitions.add(boxToIndex.getValue(state) to boxToIndex.getValue(other))
                }
            }
        }
    }

    println("Transitions: ${transitions.size}")

    return RectTransitionSystem(this, transitions)
}

fun Model.makeSemiAlgTransitions(states: List<Box>): SemiAlgTransitionSystem {
    // Important: parameters are the first dimensions since we want to project to them later
    val reducedRing = Rings.MultivariateRingQ(paramNum + varNum - 1)
    val parser = reducedRing.mkCoder(*(0 until (paramNum + varNum - 1)).map { "x$it" }.toTypedArray())
    val paramRing = Rings.MultivariateRingQ(paramNum)
    val paramSolver = SemiAlgSolver(paramBounds, paramRing)

    val paramBoundsExtended = (0 until paramNum).map { xi ->
        val low = parser.parse("x$xi - ${paramBounds.data[xi].low}")
        val high = parser.parse("x$xi - ${paramBounds.data[xi].high}")
        low to high
    }

    val transitions = ArrayList<Triple<Int, Int, SemiAlgSet>>()
    val boxToIndex = states.mapIndexed { i, box -> box to i }.toMap()
    var i = 0
    for (state in states) {
        i += 1
        if (i % 100 == 0) println("Resolve transitions $i/${states.size}")
        neighbour@ for (other in states) {
            if (state === other) continue@neighbour
            var print = false//boxToIndex.getValue(state) == 1 && boxToIndex.getValue(other) == 0
            var transitionDimension = -1
            for (d in 0 until varNum) {
                if (!state.data[d].intersects(other.data[d])) {
                    // if they do not intersect, they are not related!
                    continue@neighbour
                }
                if (state.data[d].low == other.data[d].high || state.data[d].high == other.data[d].low) {
                    if (transitionDimension > -1) {
                        // these two boxes are already touching in other dimension, they can't touch in multiple
                        continue@neighbour
                    } else {
                        transitionDimension = d
                    }
                }
            }
            // now check if the transition even can be enabled under some circumstances:
            val facetBounds = Array(varNum) { d ->
                (state.data[d] intersect other.data[d])!!
            }
            if (facetBounds[transitionDimension].low != facetBounds[transitionDimension].high) error("WTF?")
            val possibleDerivatives = equations[transitionDimension].evaluate(facetBounds, paramBounds.data)
            if (facetBounds[transitionDimension].low == state.data[transitionDimension].high) {
                // we are going up from the state, so the derivative should be positive
                if (possibleDerivatives.high > zero) {
                    if (possibleDerivatives.low > zero) {
                        // All parameters lead this way, so no CAD is needed
                        transitions.add(
                            Triple(boxToIndex.getValue(state), boxToIndex.getValue(other), paramSolver.one)
                        )
                    } else {
                        // Now we actually create the set which enables the transition:
                        // at this point, we assume the sign of the denominator is fixed, so we only care about the numerator
                        val numerator = equations[transitionDimension].numerator()
                        // this is the exact value of the transition dimension which we need to substitute
                        val commonValue = facetBounds[transitionDimension].low
                        // this is the equation that needs to be positive for our derivatives
                        val condition = numerator.eliminate(paramNum + transitionDimension, commonValue)
                        if (print) println("Condition is $condition")
                        // make a decomposition of the whole system
                        // TODO: We can make this more efficient by stripping away unused dimensions in sparse equations
                        /*val cad = LevelGraph(
                            listOf(condition),
                            reducedRing,
                            paramBounds.extend(Box(*facetBounds).eliminate(transitionDimension))
                        )
                        val cadProjected = LevelGraph(cad, paramNum, paramRing)
                        val validCells = cad.walkCells().mapNotNullTo(HashSet()) { (point, cell) ->
                            val result = condition.evaluate(*point.toTypedArray())
                            if (result <= Q.zero) null else {
                                cell.project(paramNum)
                            }
                        }
                        print = validCells.isNotEmpty()
                        if (print) println("CAD: $cad")
                        if (print) println("CAD projected: $cadProjected")*/
                        val extendedBounds = Box(*facetBounds).eliminate(transitionDimension)
                        val extendedBoundPolynomials = extendedBounds.data.mapIndexed { xi, interval ->
                            val low = parser.parse("x${xi + paramNum} - ${interval.low}")
                            val high = parser.parse("x${xi + paramNum} - ${interval.high}")
                            low to high
                        }
                        if (print) paramSolver.debug = true
                        val positiveProjected = paramSolver.positiveExtended(condition, extendedBounds, paramBoundsExtended + extendedBoundPolynomials)
                        paramSolver.run {
                            val params = positiveProjected
                            if (print) println("Resulting params: $params")
                            if (print) exitProcess(0)
                            if (params.isNotEmpty()) {
                                transitions.add(
                                    Triple(boxToIndex.getValue(state), boxToIndex.getValue(other),params)
                                )
                                /*if (params.validCells.size > 1) {
                                    println("Create a partial transition, because intervals give us $possibleDerivatives")
                                    println("Common value $commonValue")
                                    println("Transition from $state to $other with $params")
                                    //exitProcess(0)
                                }*/
                            }
                        }
                    }
                }
            } else {
                // we are going down from the state, so the derivative should be negative
                if (possibleDerivatives.low < zero) {
                    if (possibleDerivatives.high < zero) {
                        // All parameters lead this way, so no CAD is needed
                        transitions.add(
                            Triple(boxToIndex.getValue(state), boxToIndex.getValue(other), paramSolver.one)
                        )
                    } else {
                        val numerator = equations[transitionDimension].numerator()
                        val commonValue = facetBounds[transitionDimension].low
                        // this is the equation that needs to be negative for our derivatives
                        val condition = numerator.eliminate(paramNum + transitionDimension, commonValue)
                        // TODO: We can make this more efficient by stripping away unused dimensions in sparse equations
                        if (print) println("Condition: $condition")
                        val extendedBounds = Box(*facetBounds).eliminate(transitionDimension)
                        val extendedBoundPolynomials = extendedBounds.data.mapIndexed { xi, interval ->
                            val low = parser.parse("x${xi + paramNum} - ${interval.low}")
                            val high = parser.parse("x${xi + paramNum} - ${interval.high}")
                            low to high
                        }
                        if (print) paramSolver.debug = true
                        // note that we can't just negate the set, because its a different type of condition (quantifiers)
                        val positiveProjected = paramSolver.negativeExtended(condition, extendedBounds, paramBoundsExtended + extendedBoundPolynomials)
                        paramSolver.run {
                            val params = positiveProjected
                            if (print) println("Resulting params: $params")
                            if (print) exitProcess(0)
                            if (params.isNotEmpty()) {
                                transitions.add(
                                    Triple(boxToIndex.getValue(state), boxToIndex.getValue(other),params)
                                )
                                /*if (params.validCells.size > 1) {
                                    println("1Create a partial transition, because intervals give us $possibleDerivatives")
                                    println("Common value $commonValue")
                                    println("Transition from $state to $other with $params")
                                    //exitProcess(0)
                                }*/
                            }
                        }
                    }
                }
            }
        }
    }

    println("Transitions: ${transitions.size}")

    return SemiAlgTransitionSystem(paramSolver, states, transitions)
}

fun List<Box>.draw(): SvgImage {
    val rectangles = this.map {
        val (ix, iy) = it.data
        Rectangle(xy(ix.center, iy.center), xy(ix.size, iy.size), Style.STROKE)
    }
    return SvgImage(rectangles, p1)
}


fun PropertySpace.draw(): SvgImage {
    val rectangles = this.ss.mapIndexed { i, it ->
        val (ix, iy) = it.data
        Rectangle(xy(ix.center, iy.center), xy(ix.size, iy.size), if (i in valid) Style.FILL else Style.STROKE)
    }
    return SvgImage(rectangles, p1)
}

fun RectTransitionSystem.draw(): SvgImage {
    val rectangles = this.states.map {
        val (ix, iy) = it.data
        Rectangle(xy(ix.center, iy.center), xy(ix.size, iy.size), Style.STROKE)
    }
    val arrows = this.transitions.map { (iFrom, iTo) ->
        val from = rectangles[iFrom]; val to = rectangles[iTo]
        Line(from.center, to.center, Style.ARROW)
    }
    return SvgImage(rectangles + arrows, Rings.Q.parse("1/20"))
}