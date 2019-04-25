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

val ISOLATE_PRECISION = Rings.Q.parse("1/1000")

/**
 *
 * CADPolynomials represents partial information obtained when creating the cells of a CAD. It does
 * not contain all the cells, but we can iterate over them easily.
 *
 * Note that we assume that the bound polynomials are also added to the structure at this point. It is a little
 * wasteful in terms of memory, but it saves a lot of implementation/architecture time for now.
 */
data class CADPolynomials(
    /*
        The parameter bounds.
     */
    val bounds: Box,
    /*
        The (interesting) polynomials in each projection stage of the CAD. i-th element depends on all variables
        up to i. (So polynomials at position 0 contain only x1, position 3 contains x1,x2,x3,x4, etc.).

        Note that not all variables need to be necessarily used, but certainly no extra can be present (so x3 - 2
        can be on position 4, but can't be on position 1)
     */
    val levels: List<List<MPoly>>   // The polynomials in each dimension of the CAD. First element
) {

    companion object {

        /**
         * Create a CAD polynomials structure from a set of basis polynomials and a bounding box.
         *
         * The bound polynomials will be added.
         */
        fun make(ring: IPolynomialRing<MPoly>, basis: List<MPoly>, bounds: Box): CADPolynomials {
            val dimensions = bounds.data.size
            // create bound polynomials (we need these to compute intersections with bounds of the parameter space)
            val parser = ring.mkCoder(*(0 until dimensions).map { "x$it" }.toTypedArray())
            val boundPolynomials = (0 until dimensions).flatMap { d ->
                val low = parser.parse("x$d - ${bounds.data[d].low}")
                val high = parser.parse("x$d - ${bounds.data[d].high}")
                listOf(low, high)
            }
            val levels = ArrayList<List<MPoly>>()
            // we construct the levels in reversed order:
            // first, we create the top level as basis polynomials plus bound polynomials
            // at each point, we remove all polynomials which are guaranteed to be zero
            levels.add((basis + boundPolynomials).filter { it.canHaveZero(bounds) })
            // now we do projection in each dimension, but we start from the topmost dimension
            // this decision is pretty much arbitrary, but it matches what the users expect
            // (that x1 is the first variable)
            // Note that we don't project in dimension zero as there we just doo root isolation afterwards.
            for (d in (dimensions - 1)..1) {
                val levelPolynomials = levels.last()
                // Compute projection of the polynomials themselves - i.e. their outline defined by their first derivative.
                val projection1 = levelPolynomials.flatMap { p ->
                    val pUni = p.asUnivariate(d)
                    UnivariateResultants.Subresultants(pUni, pUni.derivative())
                }
                .normalize().filter { it.canHaveZero(bounds) }
                // Here we compute projection of pair-wise polynomial intersections.
                val projection2 = ArrayList<MPoly>()
                for (i in levelPolynomials.indices) {
                    for (j in (i+1) until levelPolynomials.size) {
                        val res = UnivariateResultants.Subresultants(
                            levelPolynomials[i].asUnivariate(d), levelPolynomials[j].asUnivariate(d)
                        )
                        .normalize().filter { it.canHaveZero(bounds) }
                        projection2.addAll(res)
                    }
                }
                // the final projection is the sum of the two partial projections - there can be duplicates, so remove those
                levels.add((projection1 + projection2).toSet().toList())
            }
            // finally, we reverse the levels here
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
        // we make a copy of all the polynomials because we are going to partially-evaluate them.
        val levelsCopy = levels.map { l -> l.map { p -> p.copy() }.toMutableList() }
        // here, we really go from dimension 0 up, as level 0 only contains variable x1 and level DIM all variables.
        for (d in 0 until dimensions) {
            val dimBounds = bounds.data[d]
            val levelPolynomials = levelsCopy[d].map { it.asUnivariate() }  // it this point, level d should be univariate
            val roots = DescartWithSquareFreeFactorisationRootIsolation
                .isolateRoots(levelPolynomials, ISOLATE_PRECISION)
                .toSet()    // there might be duplicates - shouldn't due to factoring, but can be
                .filter { root ->
                    // Here, we handle some edge cases - ideally, you would refine these roots further, but at this point
                    // we just hope we don't have to. If you run into this issue, just increase ISOLATE_PRECISION
                    if (root.low < dimBounds.low && root.high > dimBounds.low) error("Root crosses boundary")
                    if (root.low < dimBounds.high && root.high > dimBounds.high) error("Root crosses boundary")
                    root.low >= dimBounds.low && root.high <= dimBounds.high    // exclude roots which fall outside of the interval
                }
                .sortedBy { it.low }
            val pointValue = point[d]
            // in fact, since bound polynomials should be part of this, we should always get dimBounds as roots too
            if (roots.first().low != roots.first().high && roots.first().low != dimBounds.low) error("WTF: bounds missing")
            if (roots.last().low != roots.last().high && roots.last().low != dimBounds.high) error("WTF: bounds missing")
            // If the requested point falls too close to a root, we also fail, but we should refine the interval instead.
            if (roots.any { pointValue in it }) error("Point inside root interval :(")
            // TODO: Here, we should be using binary search because the number of roots can be high
            for (cell in 0 until (roots.size - 1))  {   // n thresholds defines n-1 cells
                val cellLow = roots[cell]
                val cellHigh = roots[cell + 1]
                if (cellLow.high <= pointValue && pointValue <= cellHigh.low) {
                    // WE FOUND IT! Save the cell id and update all levels polynomials
                    result.add(cell)
                    // partially evaluate all level polynomials at this point
                    levelsCopy.forEach { level ->
                        level.replaceAll { p -> p.evaluate(d, pointValue) }
                    }
                    // skip all other roots
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
        val dimensions = bounds.data.size
        // work stack entry contains three items:
        // 1. Partially evaluated list of levels (the ones that still need to be applied)
        // 2. Partial call coordinates
        // 3. Partial sample point
        // Note that each stack entry is actually a cell in some lower dimensional CAD, but we use this to
        // avoid storing all cells at the same time.
        val workStack = ArrayList<Triple<List<List<MPoly>>, List<Int>, List<NumQ>>>()
        workStack.add(Triple(levels, emptyList(), emptyList())) // resolve all dimensions with no initial information
        while (workStack.isNotEmpty()) {
            val (remainingLevels, partialCell, partialPoint) = workStack.removeAt(workStack.lastIndex)
            if (remainingLevels.isEmpty()) { // if there are no dimensions to process, output and continue
                yield(partialPoint to partialCell); continue
            }
            // The dimension we should be working on is the one right after the last processed
            val d = partialCell.size
            val dimBounds = bounds.data[d]
            // The first item of levels should depend only on dimension d
            val levelPolynomials = remainingLevels.first().map { it.asUnivariate() }
            val roots = DescartWithSquareFreeFactorisationRootIsolation
                .isolateRoots(levelPolynomials, ISOLATE_PRECISION)
                .toSet()    // there might be duplicates - shouldn't due to factoring, but can be
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
            for (cell in 0 until (roots.size - 1)) { // n thresholds defines n-1 cells
                val cellMiddle = roots[cell].high + ((roots[cell + 1].low - roots[cell].high) / 2)
                workStack.add(
                    Triple(
                        // partially evaluate all level polynomials at this point
                        remainingLevels.drop(1).map { level -> level.map { it.copy().evaluate(d, cellMiddle) } },
                        // add cell index to partial coordinates
                        partialCell + cell,
                        // add cell middle point to partial sample point
                        partialPoint + cellMiddle
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
fun MPoly.canHaveZero(bounds: Box): Boolean = this.evaluate(*bounds.data).hasZero