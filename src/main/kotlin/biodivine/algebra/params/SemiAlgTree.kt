package biodivine.algebra.params

import biodivine.algebra.MPoly
import biodivine.algebra.NumQ
import biodivine.algebra.ia.Interval
import biodivine.algebra.ia.div
import biodivine.algebra.ia.minus
import biodivine.algebra.ia.plus
import biodivine.algebra.rootisolation.AdaptiveRootIsolation
import biodivine.algebra.rootisolation.Root
import biodivine.algebra.synth.Box
import cc.redberry.rings.Rings
import cc.redberry.rings.Rings.Q
import java.util.concurrent.atomic.AtomicReference
import kotlin.math.max

typealias LevelList = List<Set<MPoly>>
/*
/**
 * Merge levels takes a set of bound polynomials, two level lists and a set of new polynomials
 * and merges them at dimension d.
 *
 * This consists of:
 * 1) polynomials with smaller level are transferred from add to the appropriate level
 * 2) All projections and intersections with bounds in f and g are already computed
 * (if they were removed, they were unnecessary), so we don't do any projections there.
 * 3) However, by merging f and g, new intersections can be created - we compute these.
 * 4) Finally, the newly added polynomials need to be projected and intersected with the
 * rest of f and g.
 * We recursively continue, adding into the lower level new polynomials which were created in
 * these projections and returning as the current level the merger of f,g, and new polynomials.
 */
fun mergeLevels(
    boundPolynomials: List<Pair<MPoly, MPoly>>,
    f: LevelList, g: LevelList, add: Set<MPoly>, d: Int
): LevelList {
    if (d == 0) return listOf(f[0] + g[0] + add)
    // set of polynomials which should be transferred to a lower level
    val transfer = add.filter { it.level < d }
    // set of polynomials which should be added to this level
    val new = add - transfer
    // differences between the two levels
    val extraF = f[d] - g[d]
    val extraG = g[d] - f[d]
    val projection = HashSet<MPoly>(transfer)
    // combinations of polynomials from F and G
    for (p in extraF) {
        for (q in g[d]) {
            projection.addAll(Projection.resultant(p, q, d))
        }
    }
    for (p in extraG) {
        for (q in f[d]) {
            projection.addAll(Projection.resultant(p, q, d))
        }
    }
    val old = f[d] + g[d]
    // projection of new polynomials
    new.forEach { projection.addAll(Projection.discriminant(it, d)) }
    // combination of new polynomials with old and bounds
    val (low, high) = boundPolynomials[d]
    for (p in new) {
        for (q in old) {
            projection.addAll(Projection.resultant(p, q, d))
        }
        projection.addAll(Projection.resultant(p, low, d))
        projection.addAll(Projection.resultant(p, high, d))
    }
    // remove elements which are already there
    projection.removeAll(f[d-1])
    projection.removeAll(g[d-1])
    return mergeLevels(boundPolynomials, f, g, projection, d-1) + listOf(f[d] + g[d] + new)
}

/**
 * Compute root pairs for this map of partially evaluated polynomials
 */
fun Map<MPoly, MPoly>.roots(bounds: Interval): List<Pair<Root, MPoly>> {
    return this.flatMap { (eval, original) ->
        AdaptiveRootIsolation.isolateRootsInBounds(listOf(eval.asUnivariate()), bounds)
            .map { it to original }
    }.sortedBy { it.first }
}

/**
 * Given a list of root pairs, compute the sample points defined by the system bounds and the given roots.
 */
fun List<Pair<Root, MPoly>>.samplePoints(bounds: Interval): List<NumQ> {
    val result = ArrayList<NumQ>(this.size + 1)
    val (l, h) = bounds
    if (this.isEmpty()) {
        result.add(l + (h - l)/2)
    } else {
        result.add(this.first().first.middleValue(l))
        for (i in 0 until (size - 1)) {
            result.add(this[i].first.middleValue(this[i+1].first))
        }
        result.add(this.last().first.middleValue(h))
    }
    return result
}

sealed class SemiAlgTree {

    companion object {
        fun positiveSet(poly: MPoly, bounds: Box, boundPoly: List<Pair<MPoly, MPoly>>): SemiAlgTree {
            val dimensions = bounds.data.size
            val levels = mergeLevels(
                boundPoly,
                (0..dimensions).map { emptySet<MPoly>() },
                (0..dimensions).map { emptySet<MPoly>() },
                setOf(poly),
                dimensions - 1
            )

            return makeTree(levels.map { level ->
                level.associateBy({it}, {it})
            }, bounds, emptyList(), 0, poly)
        }

        private fun makeTree(levels: List<Map<MPoly, MPoly>>, bounds: Box, point: List<NumQ>, d: Int, poly: MPoly): SemiAlgTree {
            if (levels.isEmpty()) return Leaf(d, poly.evaluate(*point.toTypedArray()) > Q.zero)
            val roots = levels.first().roots(bounds.data[d])
            val cells = roots.samplePoints(bounds.data[d]).map { sample ->
                makeTree(
                    levels.drop(1).map { level -> level.mapKeys { it.key.evaluate(d, sample) } },
                    bounds, point + sample, d + 1, poly
                )
            }

            return Cylinder(d, roots, cells)
        }
    }

    data class Leaf(override val level: Int, val member: Boolean) : SemiAlgTree() {
        override val projection: Set<MPoly>
            get() = emptySet()
        override val roots: List<Pair<Root, MPoly>>
            get() = emptyList()
        override fun similar(other: SemiAlgTree): Boolean = this == other
        override fun prune(): SemiAlgTree = this

        override fun levelPolynomials(remaining: Int): LevelList {
            return emptyLevelList (remaining)
        }

        override fun lookup(point: NumQ): SemiAlgTree = this

    }

    data class Cylinder(
        override val level: Int,
        override val roots: List<Pair<Root, MPoly>>,
        val cells: List<SemiAlgTree>
    ) : SemiAlgTree() {

        private var projectionCache = AtomicReference<Set<MPoly>?>()

        override val projection: Set<MPoly>
            get() {
                return projectionCache.get() ?: run {
                    val transfer = cells.flatMapTo(HashSet()) { cell ->
                        cell.projection.filter { it.level < this.level }
                    }
                    val myPolynomials = roots.mapTo(HashSet()) { it.second }.toList()
                    val discriminants = myPolynomials.flatMapTo(HashSet()) { poly ->
                        Projection.discriminant(poly, level)
                    }
                    val resultants = HashSet<MPoly>()
                    for (i in myPolynomials.indices) {
                        for (j in (i+1) until myPolynomials.size) {
                            Projection.resultant(myPolynomials[i], myPolynomials[j], level).forEach {
                                resultants.add(it)
                            }
                        }
                    }
                    val result = transfer + discriminants + resultants
                    projectionCache.set(result)
                    result
                }
            }

        override fun similar(other: SemiAlgTree): Boolean {
            if (other !is Cylinder) return false
            if (this.cells.size != other.cells.size) return false
            return this.cells.asSequence().zip(other.cells.asSequence()).all { (a, b) -> a == b }
        }

        override fun prune(): SemiAlgTree {
            val pruned = cells.map { it.prune() }
            if (pruned.all { it is Leaf && it.member }) return Leaf(level, true)
            if (pruned.all { it is Leaf && !it.member }) return Leaf(level, false)
            val prunedRoots = ArrayList<Pair<Root, MPoly>>(roots.size)
            val prunedCells = ArrayList<SemiAlgTree>(cells.size)
            for (i in roots.indices) {
                val r = roots[i]
                val a = pruned[i]; val b = pruned[i+1]
                if (r.second in a.projection || r.second in b.projection || !a.similar(b)) {
                    prunedRoots.add(r); prunedCells.add(a)
                }
            }
            val last = pruned.last()
            if (prunedCells.isEmpty() && last is Leaf) return Leaf(level, last.member)
            return Cylinder(level, prunedRoots, prunedCells + last)
        }

        override fun levelPolynomials(remaining: Int): LevelList {
            val parentLevels: LevelList = cells.map { it.levelPolynomials(remaining - 1) }
                .fold(emptyLevelList(remaining - 1)) { a, b -> a.zip(b) }
            return listOf<Set<MPoly>>(this.roots.mapTo(HashSet()) { it.second }) + parentLevels
        }

        override fun lookup(point: NumQ): SemiAlgTree {
            val searchResult = roots.binarySearch { (root, _) -> root.compareTo(point) }
            if (searchResult >= 0) {
                error("The point $point is also a root in $roots. This should not happen!")
            } else {
                return cells[(-searchResult - 1)]
            }
        }

    }

    abstract val projection: Set<MPoly>
    abstract val level: Int
    abstract val roots: List<Pair<Root, MPoly>>

    abstract fun levelPolynomials(remaining: Int): LevelList

    abstract fun prune(): SemiAlgTree

    abstract fun similar(other: SemiAlgTree): Boolean

    abstract fun lookup(point: NumQ): SemiAlgTree

}

fun emptyLevelList(size: Int): LevelList = (1..size).map { emptySet<MPoly>() }
*/
fun LevelList.zip(other: LevelList): LevelList {
    if (this.size < other.size) return other.zip(this)
    // this is larger
    val result = ArrayList(this)
    for (i in other.indices) {
        result[i] = result[i] + other[i]
    }
    return result
}