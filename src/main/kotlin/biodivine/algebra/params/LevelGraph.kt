package biodivine.algebra.params

import biodivine.algebra.MPoly
import biodivine.algebra.NumQ
import biodivine.algebra.ia.div
import biodivine.algebra.ia.minus
import biodivine.algebra.ia.plus
import biodivine.algebra.rootisolation.DescartWithSquareFreeFactorisationRootIsolation
import biodivine.algebra.svg.zero
import biodivine.algebra.synth.Box
import cc.redberry.rings.poly.IPolynomialRing
import cc.redberry.rings.poly.univar.UnivariateResultants

/**
 * LevelGraph is a dependency graph of polynomials, where each polynomial has a level, i.e.
 * the index of the highest variable which is still present in the polynomial (level of x1-1 is
 * 0, level of x3-x2+1 is 2, etc. and constant polynomial has level -1)
 *
 * Additionally, each polynomial can be added due to either
 * a) direct projection, in which case it has one dependency
 * b) combination projection, in which case it has two dependencies
 * c) combination of a) and b), in which case it has combination of dependencies
 *
 * Addition: add a polynomial in which case the direct projections and combination
 * projections (with the current level) are immediately computed as well.
 * Subtraction: removal of a polynomial, in which case all polynomials which depend on
 * it are removed also (as long as they don't depend on anything else).
 *
 * Note that we also automatically prune all polynomials that we can prove are sign-invariant
 * in the given bounding box.
 */
class LevelGraph(
    private val ring: IPolynomialRing<MPoly>,
    private val bounds: Box
) {

    private val dimensions: Int = bounds.data.size

    private val levels: List<MutableSet<MPoly>> = (1..dimensions).map { HashSet<MPoly>() }
    private val dependencies: MutableList<Dependency> = ArrayList()

    init {
        // Ensure dimension bounds are included
        val parser = ring.mkCoder(*(0 until dimensions).map { "x$it" }.toTypedArray())
        (0 until dimensions).forEach { d ->
            val low = parser.parse("x$d - ${bounds.data[d].low}")
            val high = parser.parse("x$d - ${bounds.data[d].high}")
            insert(low); insert(high)
        }
    }

    constructor(basis: List<MPoly>, ring: IPolynomialRing<MPoly>, bounds: Box) : this(ring, bounds) {
        basis.forEach { insert(it) }
    }

    constructor(graph1: LevelGraph, graph2: LevelGraph): this(graph1.ring, graph1.bounds) {
        if (graph1.ring != graph2.ring || graph1.bounds != graph2.bounds) error("Incompatible graphs: $graph1, $graph2")
        this.dependencies.addAll((graph1.dependencies + graph2.dependencies).toSet())   // copy dependencies
        graph1.levels.forEach { level -> level.forEach { insert(it) } }
        graph2.levels.forEach { level -> level.forEach { insert(it) } }
    }

    constructor(graph: LevelGraph, exclude: MPoly): this(graph.ring, graph.bounds) {
        if (exclude !in graph.levels[exclude.level]) error("Removing non-existent polynomial")
        // make copy of original graph
        graph.levels.forEachIndexed { l, level ->
            this.levels[l].addAll(level)
        }
        dependencies.addAll(graph.dependencies)
        // remove excluded polynomial with its dependencies
        remove(exclude)
    }

    /**
     * Computes the polynomials which have no dependencies in the level graph, i.e. the original
     * graph can be generated using this set.
     */
    val basis: List<MPoly>
        get() = levels.toList().flatten().filter { poly -> dependencies.all { it.target != poly } }

    private fun insert(poly: MPoly) {
        // note that the function is recursive, but the depth is always small (number of dimensions)
        val level = poly.level
        if (levels[level].add(poly) && level > 0) {
            // If new polynomial was added to the level and the level is more than 0, compute also dependent polynomials.
            // (if the poly was not added, dependencies are already computed and if the level is zero,
            // then there is nothing to compute)

            // Also note that some dependencies can be transferred from previous graphs, so it is possible that
            // the results are already present.

            // 1. Discriminant
            if (dependencies.find { it.checkDiscriminant(poly) } == null) {
                val uni = poly.asUnivariate(level)
                UnivariateResultants.Subresultants(uni, uni.derivative())
                    .normalize().filter { it.canHaveZero(bounds) }
                    .forEach { projectionPoly ->
                        insert(projectionPoly)
                        dependencies.add(Dependency.Single(poly, projectionPoly))
                    }
            }

            // 2. Resultants
            for (other in levels[level]) {
                if (other == poly) continue
                if (dependencies.find { it.checkCombination(poly, other) } == null) {
                    UnivariateResultants.Subresultants(poly.asUnivariate(level), other.asUnivariate(level))
                        .normalize().filter { it.canHaveZero(bounds) }
                        .forEach { projectionPoly ->
                            insert(projectionPoly)
                            dependencies.add(Dependency.Combination(poly, other, projectionPoly))
                        }
                }
            }
        }
    }

    private fun remove(poly: MPoly) {
        val level = poly.level
        if (levels[level].remove(poly) && level > 0) {
            // Similar to insert: Only continue if something was really removed. If level is zero, there should be no
            // dependencies to remove.

            val ourDependencies = this.dependencies.filter { it.isDependencyOf(poly) }
            dependencies.removeIf { it.isDependencyOf(poly) }

            ourDependencies.forEach { removed ->
                if (dependencies.all { it.target != removed.target }) {
                    // All dependencies are removed, we should also get rid of this target
                    remove(removed.target)
                }
            }
        }
    }

    fun cellForPoint(point: List<NumQ>): Cell {
        val result = IntArray(dimensions)
        val evalLevels = levels.map { l -> l.map { it.copy() }.toMutableList() }
        for (d in 0 until dimensions) {
            val interval = bounds.data[d]
            val levelPolynomials = evalLevels[d].map { it.asUnivariate() }  // it this point, level d should be univariate
            val roots = DescartWithSquareFreeFactorisationRootIsolation
                .isolateRoots(levelPolynomials, ISOLATE_PRECISION)
                .toSet()    // there might be duplicates - shouldn't due to factoring, but can be
                .filter { root ->
                    // Here, we handle some edge cases - ideally, you would refine these roots further, but at this point
                    // we just hope we don't have to. If you run into this issue, just increase ISOLATE_PRECISION
                    if (root.low < interval.low && root.high > interval.low) error("Root crosses boundary")
                    if (root.low < interval.high && root.high > interval.high) error("Root crosses boundary")
                    root.low >= interval.low && root.high <= interval.high    // exclude roots which fall outside of the interval
                }
                .sortedBy { it.low }
            val pointValue = point[d]
            // in fact, since bound polynomials should be part of this, we should always get dimBounds as roots too
            if (roots.first().low != roots.first().high && roots.first().low != interval.low) error("WTF: bounds missing")
            if (roots.last().low != roots.last().high && roots.last().low != interval.high) error("WTF: bounds missing")
            // If the requested point falls too close to a root, we also fail, but we should refine the interval instead.
            if (roots.any { pointValue in it }) error("Point inside root interval :(")
            // TODO: Here, we should be using binary search because the number of roots can be high
            for (cell in 0 until (roots.size - 1))  {   // n thresholds defines n-1 cells
                val cellLow = roots[cell]
                val cellHigh = roots[cell + 1]
                if (cellLow.high <= pointValue && pointValue <= cellHigh.low) {
                    // WE FOUND IT! Save the cell id and update all levels polynomials
                    result[d] = cell
                    // partially evaluate all level polynomials at this point
                    evalLevels.forEach { level ->
                        level.replaceAll { p -> p.evaluate(d, pointValue) }
                    }
                    // skip all other roots
                    break
                }
            }
        }
        return Cell(result)
    }

    /**
     * Provides an iterator over the cells of this (bounded) CAD. The cell is represented by its rational sample point
     * and its coordinates.
     */
    fun walkCells(): Sequence<Pair<List<NumQ>, Cell>> = sequence {
        // work stack entry contains three items:
        // 1. Partially evaluated list of levels (the ones that still need to be applied)
        // 2. Partial call coordinates
        // 3. Partial sample point
        // Note that each stack entry is actually a cell in some lower dimensional CAD, but we use this to
        // avoid storing all cells at the same time.
        val workStack = ArrayList<Triple<List<List<MPoly>>, List<Int>, List<NumQ>>>()
        workStack.add(Triple(levels.map { it.toList() }, emptyList(), emptyList()))
        while (workStack.isNotEmpty()) {
            val (remainingLevels, partialCell, partialPoint) = workStack.removeAt(workStack.lastIndex)
            if (remainingLevels.isEmpty()) { // if there are no dimensions to process, output and continue
                yield(partialPoint to Cell(partialCell.toIntArray())); continue
            }
            // The dimension we should be working on is the one right after the last processed
            val d = dimensions - remainingLevels.size
            val interval = bounds.data[d]
            // The first item of levels should depend only on dimension d
            val levelPolynomials = remainingLevels.first().map { it.asUnivariate() }
            val roots = DescartWithSquareFreeFactorisationRootIsolation
                .isolateRoots(levelPolynomials, ISOLATE_PRECISION)
                .toSet()    // there might be duplicates - shouldn't due to factoring, but can be
                .filter { root ->
                    if (root.low < interval.low && root.high > interval.low) error("Root crosses boundary")
                    if (root.low < interval.high && root.high > interval.high) error("Root crosses boundary")
                    root.low >= interval.low && root.high <= interval.high
                }
                .sortedBy { it.low }
            // in fact, since bound polynomials should be part of this, we should always get dimBounds as roots too
            if (roots.first().isNumber() && roots.first().low != interval.low) error("WTF: bounds missing")
            if (roots.last().isNumber() && roots.last().low != interval.high) error("WTF: bounds missing")
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


    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as LevelGraph

        if (bounds != other.bounds) return false
        if (levels != other.levels) return false
        if (dependencies != other.dependencies) return false

        return true
    }

    override fun hashCode(): Int {
        var result = bounds.hashCode()
        result = 31 * result + levels.hashCode()
        result = 31 * result + dependencies.hashCode()
        return result
    }

    override fun toString(): String {
        return "LevelGraph(levels=$levels)"
    }

    /**
     * The highest variable which has a non-zero degree in the polynomial. Note that this is also the index
     * of a variable where we want to do projections!
     */
    private val MPoly.level: Int
        get() {
            for (d in (nVariables - 1) downTo 0) {
                if (degree(d) > 0) {
                    // this is the highest d such that there is a non-zero degree there, so this is the level
                    return d
                }
            }
            return -1
        }



    sealed class Dependency {
        data class Single(private val source: MPoly, override val target: MPoly) : Dependency() {
            override fun checkDiscriminant(source: MPoly): Boolean = this.source == source
            override fun checkCombination(source1: MPoly, source2: MPoly): Boolean = false
            override fun isDependencyOf(poly: MPoly): Boolean = source == poly
        }

        data class Combination(private val source1: MPoly, private val source2: MPoly, override val target: MPoly) : Dependency() {
            override fun checkDiscriminant(source: MPoly): Boolean = false
            override fun checkCombination(source1: MPoly, source2: MPoly): Boolean {
                return (this.source1 == source1 && this.source2 == source2) || (this.source1 == source2 && this.source2 == source1)
            }

            override fun isDependencyOf(poly: MPoly): Boolean = source1 == poly || source2 == poly
        }

        abstract fun checkDiscriminant(source: MPoly): Boolean
        abstract fun checkCombination(source1: MPoly, source2: MPoly): Boolean
        abstract fun isDependencyOf(poly: MPoly): Boolean
        abstract val target: MPoly
    }

}

operator fun LevelGraph.plus(other: LevelGraph) = LevelGraph(this, other)
operator fun LevelGraph.minus(exclude: MPoly) = LevelGraph(this, exclude)