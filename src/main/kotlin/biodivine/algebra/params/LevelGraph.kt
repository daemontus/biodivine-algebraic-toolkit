package biodivine.algebra.params

import biodivine.algebra.MPoly
import biodivine.algebra.NumQ
import biodivine.algebra.rootisolation.AdaptiveRootIsolation
import biodivine.algebra.synth.Box
import cc.redberry.rings.Rings
import cc.redberry.rings.poly.IPolynomialRing
import cc.redberry.rings.poly.univar.UnivariateResultants

var maxLevels = 0

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
class LevelGraph private constructor(
    private val ring: IPolynomialRing<MPoly>,
    private val bounds: Box
) {

    private val dimensions: Int = bounds.data.size

    private val levels: List<MutableSet<MPoly>> = (1..dimensions).map { HashSet<MPoly>() }
    private val dependencies: MutableList<Dependency> = ArrayList()

    private val boundPolynomials: MutableSet<MPoly> = HashSet()

    init {
        // Ensure dimension bounds are included
        val parser = ring.mkCoder(*(0 until dimensions).map { "x$it" }.toTypedArray())
        (0 until dimensions).forEach { d ->
            val low = parser.parse("x$d - ${bounds.data[d].low}")
            val high = parser.parse("x$d - ${bounds.data[d].high}")
            insert(low); insert(high)
            boundPolynomials.add(low)
            boundPolynomials.add(high)
        }
    }

    constructor(basis: List<MPoly>, ring: IPolynomialRing<MPoly>, bounds: Box) : this(ring, bounds) {
        basis.forEach { insert(it) }

        if (levels[0].size > maxLevels) {
            maxLevels = levels[0].size
            println("New Max: $maxLevels $this")
        }
    }

    constructor(graph1: LevelGraph, graph2: LevelGraph): this(graph1.ring, graph1.bounds) {
        if (graph1.ring != graph2.ring || graph1.bounds != graph2.bounds) error("Incompatible graphs: $graph1, $graph2; ${graph1.ring} ${graph2.ring}; ${graph1.bounds} ${graph2.bounds}")
        this.dependencies.addAll((graph1.dependencies + graph2.dependencies).toSet())   // copy dependencies
        graph1.levels.forEach { level -> level.forEach { insert(it) } }
        graph2.levels.forEach { level -> level.forEach { insert(it) } }

        if (levels[0].size > maxLevels) {
            maxLevels = levels[0].size
            println("New Max: $maxLevels $this")
        }
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

        if (levels[0].size > maxLevels) {
            maxLevels = levels[0].size
            println("New Max: $maxLevels $this")
        }
    }

    /**
     * Produce a level graph which is a projection of this graph to the first [retain] variables.
     */
    constructor(graph: LevelGraph, retain: Int, reducedRing: IPolynomialRing<MPoly>? = null): this(
        reducedRing ?: Rings.MultivariateRingQ(retain), graph.bounds.project(retain)
    ) {
        // copy levels, but only the ones we want to retain and set their number of variables accordingly
        this.levels.forEachIndexed { l, level ->
            graph.levels[l].forEach { poly ->
                level.add(poly.setNVariables(retain))
            }
        }

        // Now this is a little nasty, because for each dependency, we also have to drop the number of variables
        // in each polynome
        // Only keep dependencies where the source level is still applicable to projected data set.
        val newDependencies = graph.dependencies.mapNotNull {
            when (it) {
                is Dependency.Single -> {
                    if (it.source.level >= retain) null else {
                        Dependency.Single(it.source.setNVariables(retain), it.target.setNVariables(retain))
                    }
                }
                is Dependency.Combination -> {
                    if (it.source1.level >= retain || it.source2.level >= retain) null else {
                        Dependency.Combination(
                            it.source1.setNVariables(retain),
                            it.source2.setNVariables(retain),
                            it.target.setNVariables(retain))
                    }
                }
            }
        }
        dependencies.addAll(newDependencies)

        if (levels[0].size > maxLevels) {
            maxLevels = levels[0].size
            println("New Max: $maxLevels $this")
        }
    }

    /**
     * Computes the polynomials which have no dependencies in the level graph, i.e. the original
     * graph can be generated using this set.
     */
    val basis: List<MPoly>
        get() = levels.toList().flatten().filter { poly ->
            poly !in boundPolynomials && dependencies.all { it.target != poly }
        }

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
            val levelPolynomials = evalLevels[d].map { it.asUnivariate() }  // d-level has always one variable remaining
            val pointValue = point[d]
            val roots = AdaptiveRootIsolation.isolateRootsInBounds(levelPolynomials, interval).toList() // values are already sorted
            val searchResult = roots.binarySearch { root -> root.compareTo(pointValue) }
            if (searchResult >= 0) {
                error("The point $pointValue of $point is also a root in $roots. This should not happen!")
            } else {
                // (-result - 1) is index of first element which is larger than pointValue, hence the index
                // of our cell is the index of one root below (we should never get that insertion index
                // is 0 or roots.size, because that means the point is outside the variable bounds
                result[d] = (-searchResult - 1) - 1
                for (e in (d+1) until dimensions) {
                    evalLevels[e].replaceAll { p -> p.evaluate(d, pointValue) }
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
            val roots = AdaptiveRootIsolation.isolateRootsInBounds(levelPolynomials, interval).toList()
            for (cell in 0 until (roots.size - 1)) { // n thresholds defines n-1 cells
                val cellMiddle = roots[cell].middleValue(roots[cell + 1])
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
        //return "LevelGraph(levels=$levels, dependencies=$dependencies)"
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
        data class Single(val source: MPoly, override val target: MPoly) : Dependency() {
            override fun checkDiscriminant(source: MPoly): Boolean = this.source == source
            override fun checkCombination(source1: MPoly, source2: MPoly): Boolean = false
            override fun isDependencyOf(poly: MPoly): Boolean = source == poly
        }

        data class Combination(val source1: MPoly, val source2: MPoly, override val target: MPoly) : Dependency() {
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