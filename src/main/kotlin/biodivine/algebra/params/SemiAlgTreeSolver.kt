package biodivine.algebra.params

import biodivine.algebra.MPoly
import biodivine.algebra.NumQ
import biodivine.algebra.ia.Interval
import biodivine.algebra.rootisolation.AdaptiveRootIsolation
import biodivine.algebra.rootisolation.Root
import biodivine.algebra.synth.Box
import cc.redberry.rings.Rings
import cc.redberry.rings.Rings.Q
import cc.redberry.rings.poly.IPolynomialRing
import cc.redberry.rings.poly.multivar.MultivariatePolynomial

class SemiAlgTreeSolver(
    private val boundBox: Box,
    private val ring: IPolynomialRing<MPoly>
) {

    private val dimensions = boundBox.data.size
    private val parser = ring.mkCoder(*(0 until dimensions).map { "x$it" }.toTypedArray())

    private val boundPolynomials: List<Pair<MPoly, MPoly>> = (0 until dimensions).map { d ->
        val low = parser.parse("x$d - ${boundBox.data[d].low}")
        val high = parser.parse("x$d - ${boundBox.data[d].high}")
        low to high
    }

    private val boundRoots: List<Pair<Root, Root>> = boundPolynomials.map { (l, h) ->
        Root.linear(l.asUnivariate()) to Root.linear(h.asUnivariate())
    }

    infix fun Tree.and(other: Tree): Tree {
        val merged = apply(emptyList(), this, other) { a, b -> a && b }
        return merged
    }

    private fun makePositive(point: List<NumQ>, poly: MPoly, levels: LevelList): Tree {
        val level = point.size
        if (level == dimensions) return Tree.Leaf(level, poly.evaluate(*point.toTypedArray()) > Q.zero)
        val roots = levels[level].flatMapTo(HashSet()) { poly ->
                AdaptiveRootIsolation.isolateRootsInBounds(
                    listOf(poly.evaluate(point).asUnivariate()), boundBox.data[level]
                ).map { it to poly }
        }.sortedBy { it.first }
        val cells = ArrayList<Tree>(roots.size + 1)
        roots.rootPairs(level) { (l, h) ->
            val sample = samplePoint(l.first, h.first)
            cells.add(makePositive(point + sample, poly, levels))
        }
        return Tree.Cylinder(level, roots, cells)
    }

    private fun apply(point: List<NumQ>, a: Tree, b: Tree, op: (Boolean, Boolean) -> Boolean): Tree {
        val level = point.size
        if (a is Tree.Leaf && b is Tree.Leaf) return Tree.Leaf(level, op(a.member, b.member))
        val polyA = a.roots.map { it.second }.toSet()
        val polyB = b.roots.map { it.second }.toSet()
        val mergedPolynomials = polyA + polyB
        // due to possible split in lower dimensions, we have to re-evaluate these roots
        val mergedRoots = mergedPolynomials.flatMapTo(HashSet()) { poly ->
            AdaptiveRootIsolation.isolateRootsInBounds(
                listOf(poly.evaluate(point).asUnivariate()), boundBox.data[level]
            ).map { it to poly }
        }.sortedBy { it.first }
        val (boundLow, _) = boundRoots[level]
        val roots = ArrayList<Pair<Root, MPoly>>()
        val cells = ArrayList<Tree>()
        var iA = 0
        var iB = 0
        mergedRoots.rootPairs(level) { (l, h) ->
            val childA = a.lookup(iA)
            val childB = b.lookup(iB)
            if (h.second in polyA) iA += 1
            if (h.second in polyB) iB += 1
            if (childA is Tree.Leaf && childB is Tree.Leaf) {
                // we know there are no level polynomials above us - can make the op right now
                if (l.first != boundLow) roots.add(l)
                cells.add(Tree.Leaf(level + 1, op(childA.member, childB.member)))
            } else {
                val intersectionPolynomials = intersectCylinder(0, childA.levelList, childB.levelList)
                val newRoots = intersectionPolynomials.flatMapTo(HashSet()) { poly ->
                    AdaptiveRootIsolation.isolateRootsInBounds(
                        listOf(poly.evaluate(point).asUnivariate()), l.first, h.first
                    ).map { it to poly }
                }.sortedBy { it.first }
                newRoots.rootPairs(l, h) { (ll, hh) ->
                    val s = samplePoint(ll.first, hh.first)
                    if (ll.first != boundLow) roots.add(ll)
                    cells.add(apply(point + s, childA, childB, op))
                }
            }
        }
        return Tree.Cylinder(level, roots, cells)
    }

    private fun intersectCylinder(rootLevel: Int, a: LevelList, b: LevelList, i: Int = 0): Set<MPoly> {
        return when {
            i >= a.size && i >= b.size -> emptySet()
            i >= a.size -> b[i] // from this level on, there is only b, so b is correct at this point
            i >= b.size -> a[i]
            else -> {
                val level = rootLevel + i
                val (low, high) = boundPolynomials[level]
                val left = a[i]
                val right = b[i]
                val current = left + right
                val intersections = intersectCylinder( rootLevel, a, b, i+1) - current
                val carry = intersections.filter { it.level < level }
                val insert = intersections - carry
                val result = HashSet<MPoly>(carry)  // carry goes automatically down
                // project new polynomials into lower dimension (including bounds)
                for (p in insert) {
                    result.addAll(Projection.discriminant(p, level))
                    result.addAll(Projection.resultant(p, low, level))
                    result.addAll(Projection.resultant(p, high, level))
                    // find intersections of insert and current polynomials
                    for (c in current) {
                        result.addAll(Projection.resultant(p, c, level))
                    }
                }
                // find intersections of polynomials which are extra on the left with everything on the right
                for (l in left) {
                    if (l in right) continue
                    for (r in right) {
                        result.addAll(Projection.resultant(l, r, level))
                    }
                }
                for (r in right) {
                    if (r in left) continue
                    for (l in left) {
                        result.addAll(Projection.resultant(l, r, level))
                    }
                }
                result
            }
        }
    }

    private inline fun List<Pair<Root, MPoly>>.rootPairs(
        low: Pair<Root, MPoly>, high: Pair<Root, MPoly>,
        action: (Pair<Pair<Root, MPoly>, Pair<Root, MPoly>>) -> Unit
    ) {
        val (rootL, polyL) = low
        val (rootH, polyH) = high
        if (this.isEmpty()) {
            action((rootL to polyL) to (rootH to polyH))
        } else {
            action((rootL to polyL) to first())
            for (i in 0 until (size - 1)) {
                action(this[i] to this[i+1])
            }
            action(last() to (rootH to polyH))
        }
    }

    private inline fun List<Pair<Root, MPoly>>.rootPairs(
        level: Int, action: (Pair<Pair<Root, MPoly>, Pair<Root, MPoly>>) -> Unit
    ) {
        val (rootL, rootH) = boundRoots[level]
        val (polyL, polyH) = boundPolynomials[level]
        if (this.isEmpty()) {
            action((rootL to polyL) to (rootH to polyH))
        } else {
            action((rootL to polyL) to first())
            for (i in 0 until (size - 1)) {
                action(this[i] to this[i+1])
            }
            action(last() to (rootH to polyH))
        }
    }

    private fun samplePoint(a: Root, b: Root): NumQ {
        return a.middleValue(b)
    }

    private fun makeLevelList(set: Set<MPoly>, level: Int = dimensions - 1): LevelList {
        if (level == 0) return listOf(set)
        val carry = set.filter { it.level < level }
        val thisLevel = set - carry
        val list = thisLevel.toList()
        val next = HashSet(carry)
        val (low, high) = boundPolynomials[level]
        for (p in thisLevel) {
            next.addAll(Projection.discriminant(p, level))
            next.addAll(Projection.resultant(p, low, level))
            next.addAll(Projection.resultant(p, high, level))
        }
        for (i in list.indices) {
            for (j in (i+1)..list.lastIndex) {
                next.addAll(Projection.resultant(list[i], list[j], level))
            }
        }
        return makeLevelList(next, level - 1) + listOf(thisLevel.toSet())
    }

    fun positive(poly: MPoly): Tree {
        val levels = makeLevelList(setOf(poly))
        return makePositive(emptyList(), poly, levels)
    }

    sealed class Tree {

        abstract val level: Int
        abstract val roots: List<Pair<Root, MPoly>>
        abstract val levelList: List<Set<MPoly>>

        abstract fun lookup(point: NumQ): Tree
        abstract fun lookup(cell: Int): Tree

        data class Leaf(override val level: Int, val member: Boolean) : Tree() {

            override val levelList: List<Set<MPoly>>
                get() = emptyList()

            override val roots: List<Pair<Root, MPoly>>
                get() = emptyList()

            override fun lookup(point: NumQ): Tree = this
            override fun lookup(cell: Int): Tree = this
        }

        data class Cylinder(
            override val level: Int,
            override val roots: List<Pair<Root, MPoly>>,
            val cells: List<Tree>
        ) : Tree() {

            override val levelList: List<Set<MPoly>> =
                listOf<Set<MPoly>>(roots.mapTo(HashSet()) { it.second }) + cells
                    .map { it.levelList }
                    .fold1 { a, b -> a.zip(b) }

            override fun lookup(point: NumQ): Tree {
                val searchResult = roots.binarySearch { (root, _) -> root.compareTo(point) }
                if (searchResult >= 0) {
                    error("The point $point is also a root in $roots. This should not happen!")
                } else {
                    return cells[(-searchResult - 1)]
                }
            }

            override fun lookup(cell: Int): Tree = cells[cell]

        }

    }

}

fun <T> List<T>.fold1(action: (T, T) -> T): T {
    var result = first()
    for (i in 1 until size) {
        result = action(result, this[i])
    }
    return result
}

private val evalArrays: List<IntArray> = (0..20).map { vars -> IntArray(vars) { it } }

fun MPoly.evaluate(point: List<NumQ>): MultivariatePolynomial<NumQ> {
    if (point.isEmpty()) return this
    return this.evaluate(evalArrays[point.size], point.toTypedArray())
}

fun main() {
    val dimensions = 2
    val ring = Rings.MultivariateRingQ(dimensions)
    val bounds = Box(
        Interval(Q.parse("0"), Q.parse("3")),
        Interval(Q.parse("0"), Q.parse("2"))
    )
    SemiAlgTreeSolver(bounds, ring).run {
        val poly = ring.parse("(y-1/2)^2 - x + 1/2")
        val poly2 = ring.parse("(x-1/2)^2 - y - 1/2")
        val set1 = positive(poly)
        val set2 = positive(poly2)
        println("Set: $set1")
        //println("Simplified: ${set1.prune()}")
        println("Set2: $set2")
        //println("Simplified: ${set2.prune()}")
        println("Intersection: ${set1 and set2}")
    }
}