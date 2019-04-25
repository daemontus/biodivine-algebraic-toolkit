package biodivine.algebra.params

import biodivine.algebra.MPoly
import biodivine.algebra.NumQ
import biodivine.algebra.ia.*
import biodivine.algebra.rootisolation.DescartWithSquareFreeFactorisationRootIsolation
import biodivine.algebra.synth.Box
import cc.redberry.rings.Rings
import cc.redberry.rings.poly.IPolynomialRing
import cc.redberry.rings.poly.PolynomialMethods.Factor
import cc.redberry.rings.poly.univar.UnivariateResultants

/**
 * CADPolynomials represents partial information obtained when creating the cells of a CAD. It does
 * not contain all the cells, but we can iterate over them easily.
 *
 * Note that we assume that the bound polynomials are also added to the structure at this point. It is a little
 * wasteful in terms of memory, but it saves a lot of implementation/architecture time for now.
 */
data class CADPolynomials(
    val bounds: Box,
    val levels: List<List<MPoly>>
) {

    companion object {

        /**
         * Create a CAD polynomials structure from a set of basis polynomials and a bounding box.
         *
         * The bound polynomials will be added.
         */
        fun make(ring: IPolynomialRing<MPoly>, basis: List<MPoly>, bounds: Box): CADPolynomials {
            val dimensions = bounds.data.size
            val parser = ring.mkCoder(*(0 until dimensions).map { "x$it" }.toTypedArray())
            val boundPolynomials = (0 until dimensions).flatMap { d ->
                val low = parser.parse("x$d - ${bounds.data[d].low}")
                val high = parser.parse("x$d - ${bounds.data[d].high}")
                listOf(low, high)
            }
            val levels = ArrayList<List<MPoly>>()
            levels.add((basis + boundPolynomials).filter {
                it.evaluate(*bounds.data).hasZero
            })
            for (d in 0 until (dimensions - 1)) {
                val levelPolynomials = levels.last()
                val projection1 = levelPolynomials.flatMap { p ->
                    val pUni = p.asUnivariate(d)
                    UnivariateResultants.Subresultants(pUni, pUni.derivative())
                }
                .filter { it.evaluate(*bounds.data).hasZero }
                .normalize()
                .filter { it.evaluate(*bounds.data).hasZero }
                val projection2 = ArrayList<MPoly>()
                for (i in levelPolynomials.indices) {
                    for (j in (i+1) until levelPolynomials.size) {
                        val res = UnivariateResultants.Subresultants(
                            levelPolynomials[i].asUnivariate(d), levelPolynomials[j].asUnivariate(d)
                        )
                            .filter {it.evaluate(*bounds.data).hasZero}
                            .normalize()
                            .filter {it.evaluate(*bounds.data).hasZero}
                        projection2.addAll(res)
                    }
                }
                levels.add((projection1 + projection2).toSet().toList())
            }
            return CADPolynomials(bounds, levels.reversed())
        }

    }

    /**
     * For a given (rational) point, find the coordinates of the cell where this point belongs,
     * or error if the point is too close to a root.
     */
    fun cellForPoint(point: List<NumQ>): List<Int> {
        val dimensions = bounds.data.size
        val result = ArrayList<Int>()
        val levelsCopy = levels.map { l -> l.map { p -> p.copy() }.toMutableList() }
        for (d in 0 until dimensions) {
            val dimBounds = bounds.data[d]
            val levelPolynomials = levelsCopy[d].map { it.asUnivariate() }  // it this point, level d should be univariate
            val roots = DescartWithSquareFreeFactorisationRootIsolation
                .isolateRoots(levelPolynomials, Rings.Q.parse("1/1000"))
                .filter { root ->
                    if (root.low < dimBounds.low && root.high > dimBounds.low) error("Root crosses boundary")
                    if (root.low < dimBounds.high && root.high > dimBounds.high) error("Root crosses boundary")
                    root.low >= dimBounds.low && root.high <= dimBounds.high
                }
                .sortedBy { it.low }
            val pointValue = point[d]
            // in fact, since bound polynomials should be part of this, we should always get dimBounds as roots too
            if (roots.first().low != roots.first().high && roots.first().low != dimBounds.low) error("WTF: bounds missing")
            if (roots.last().low != roots.last().high && roots.last().low != dimBounds.high) error("WTF: bounds missing")
            if (roots.any { pointValue in it }) error("Point inside root interval :(")
            for (cell in 0 until (roots.size - 1))  {   // for two roots, there is one cell, for three, there are two
                val cellLow = roots[cell]
                val cellHigh = roots[cell + 1]
                if (cellLow.high <= pointValue && pointValue <= cellHigh.low) {
                    // WE FOUND IT! Save the cell id and update all levels polynomials
                    result.add(cell)
                    levelsCopy.forEach { level ->
                        level.replaceAll { p -> p.evaluate(d, pointValue) }
                    }
                    break
                }
            }
        }
        return result
    }

    /**
     * Provides an iterator over the cells of this (bounded) CAD. The cell is represented by its rational sample point
     * and its coordinates.
     */
    fun walkCells(): Iterator<Pair<List<NumQ>, List<Int>>> = iterator {
        val workStack = ArrayList<Triple<List<List<MPoly>>, List<Int>, List<NumQ>>>()
        workStack.add(Triple(levels, emptyList(), emptyList()))
        while (workStack.isNotEmpty()) {
            val (levels, coordinates, samplePoint) = workStack.removeAt(workStack.lastIndex)
            if (levels.isEmpty()) {
                yield(samplePoint to coordinates); continue
            }
            val d = bounds.data.size - coordinates.size - 1
            val dimBounds = bounds.data[d]
            val levelPolynomials = levels.first().map { it.asUnivariate() }
            val roots = DescartWithSquareFreeFactorisationRootIsolation
                .isolateRoots(levelPolynomials, Rings.Q.parse("1/1000"))
                .toSet()
                .filter { root ->
                    if (root.low < dimBounds.low && root.high > dimBounds.low) error("Root crosses boundary")
                    if (root.low < dimBounds.high && root.high > dimBounds.high) error("Root crosses boundary")
                    root.low >= dimBounds.low && root.high <= dimBounds.high
                }
                .sortedBy { it.low }
            println("Roots at level $d: $roots")
            // in fact, since bound polynomials should be part of this, we should always get dimBounds as roots too
            if (roots.first().low != roots.first().high && roots.first().low != dimBounds.low) error("WTF: bounds missing")
            if (roots.last().low != roots.last().high && roots.last().low != dimBounds.high) error("WTF: bounds missing")
            for (cell in 0 until (roots.size - 1)) {
                val cellMiddle = roots[cell].high + ((roots[cell + 1].low - roots[cell].high) / 2)
                workStack.add(
                    Triple(
                        levels.drop(1).map { it.map { it.evaluate(d, cellMiddle) } },
                        coordinates + cell,
                        samplePoint + cellMiddle
                    )
                )
            }
        }
    }

}

fun main() {
    val ring = Rings.MultivariateRingQ(2)
    val p = ring.parse("x^2-y")
    val q = ring.parse("x^2 - 4*x + 4 - y")
    val bounds = Box(Interval(0, 2), Interval(0, 2))
    val cad = CADPolynomials.make(ring, listOf(p, q), bounds)
    println("CAD Polynomials: $cad")
    cad.walkCells().forEach { (point, coordinates) ->
        println("Point $point with id $coordinates")
    }

    println("Cell ${cad.cellForPoint(listOf(Rings.Q.parse("1/2"), Rings.Q.parse("3/2")))}")
}

fun List<MPoly>.normalize(): List<MPoly> = this.flatMap { Factor(it) }.map { it.divideByLC(it) }