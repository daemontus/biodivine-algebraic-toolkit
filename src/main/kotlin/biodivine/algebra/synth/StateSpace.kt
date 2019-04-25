package biodivine.algebra.synth

import biodivine.algebra.NumQ
import biodivine.algebra.ia.div
import biodivine.algebra.ia.evaluate
import biodivine.algebra.ia.plus
import biodivine.algebra.ia.times
import biodivine.algebra.svg.*
import cc.redberry.rings.Rings

fun Model.computeStateSpace(initial: List<Box>, maxDivision: NumQ): List<Box> {
    val result = ArrayList<Box>()
    val bounds = this.varBounds
    val smallest = bounds.volume * maxDivision
    var toResolve = initial
    var safeVolume = Rings.Q.zero
    var unsafeVolume = Rings.Q.zero
    while (toResolve.isNotEmpty()) {
        val nextResolve = ArrayList<Box>()
        for (box in toResolve) {
            when {
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