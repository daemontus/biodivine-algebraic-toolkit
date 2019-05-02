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
import cc.redberry.rings.poly.univar.UnivariateResultants

typealias SemiAlgSet = SemiAlgTreeSolver.Tree
typealias SemiAlgSolver = SemiAlgTreeSolver

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

    val one: Tree = Tree.Leaf(0, true)
    val zero: Tree = Tree.Leaf(0, false)

    infix fun Tree.subset(other: Tree): Boolean {
        return apply(emptyList(), this, other) { a, b -> !a || b }.all()
    }

    var debug = false

    infix fun Tree.and(other: Tree): Tree {
        return apply(emptyList(), this, other) { a, b -> a && b }.prune()/*.also {
            if (!(it subset this) || !(it subset other)) {
                println("First: $this")
                println("Second: $other")
                println("Result: $it")
                val full = apply(emptyList(), this, other) { a, b -> a && b }
                println("Unpruned: ${full}")
                debug = true
                val pruned = full.prune()
                println("Pruned: $pruned")
                error("")
            }
        }*/
    }

    infix fun Tree.or(other: Tree): Tree {
        return apply(emptyList(), this, other) { a, b -> a || b }.prune()/*.also {
            if (!(this subset it) || !(other subset it)) error("")
        }*/
    }

    fun Tree.not(): Tree = when (this) {
        is Tree.Leaf -> Tree.Leaf(level, member.not())
        is Tree.Cylinder -> Tree.Cylinder(level, roots, cells.map { it.not() })
    }

    fun Tree.isEmpty(): Boolean = !this.any()
    fun Tree.isNotEmpty(): Boolean = this.any()

    private fun makePositive(point: List<NumQ>, poly: MPoly, levels: LevelList): Tree {
        val level = point.size
        if (level == dimensions) return Tree.Leaf(level, poly.evaluate(*point.toTypedArray()) > Q.zero)
        val roots = levels[level].findMRoots(point, boundBox.data[level])
        val cells = ArrayList<Tree>(roots.size + 1)
        roots.rootPairs(level) { (l, h) ->
            val sample = samplePoint(l.root, h.root)
            cells.add(makePositive(point + sample, poly, levels))
        }
        return Tree.Cylinder(level, roots, cells)
    }

    private fun apply(point: List<NumQ>, a: Tree, b: Tree, op: (Boolean, Boolean) -> Boolean): Tree {
        val level = point.size
        if (a is Tree.Leaf && b is Tree.Leaf) return Tree.Leaf(level, op(a.member, b.member))
        val polyA = a.roots.map { it.poly }.toSet()
        val polyB = b.roots.map { it.poly }.toSet()
        // due to possible split in lower dimensions, we have to re-evaluate these roots
        val mergedRoots = (a.roots.moveTo(point, boundBox.data[level]) + b.roots.moveTo(point, boundBox.data[level]))
            .toSet().sortedBy { it.root }
        if (debug) println("Merged roots: ${mergedRoots}")
        val (boundLow, _) = boundRoots[level]
        val roots = ArrayList<MRoot>()
        val cells = ArrayList<Tree>()
        var iA = 0
        var iB = 0
        mergedRoots.rootPairs(level) { (l, h) ->
            if (debug) println("Merged cell $l $h")
            val childA = a.lookup(iA)
            val childB = b.lookup(iB)
            if (h.poly in polyA) iA += 1
            if (h.poly in polyB) iB += 1
            if (childA is Tree.Leaf && childB is Tree.Leaf) {
                // we know there are no level polynomials above us - can make the op right now
                if (l.root != boundLow) roots.add(l)
                cells.add(Tree.Leaf(level + 1, op(childA.member, childB.member)))
            } else {
                val intersectionPolynomials = intersectCylinder(level + 1, childA.levelList, childB.levelList)
                if (debug) println("Intersection: $intersectionPolynomials")
                val newRoots = intersectionPolynomials.findMRoots(point, boundBox.data[level], l.root, h.root)
                newRoots.rootPairs(l, h) { (ll, hh) ->
                    if (debug) println("New cell $ll $hh")
                    val s = samplePoint(ll.root, hh.root)
                    if (ll.root != boundLow) roots.add(ll)
                    cells.add(apply(point + s, childA, childB, op))
                }
            }
        }
        return Tree.Cylinder(level, roots, cells)
    }

    private fun Tree.prune(): Tree = when (this) {
        is Tree.Leaf -> this
        is Tree.Cylinder -> {
            val pruned = this.cells.map { it.prune() }
            when {
                pruned.all { it is Tree.Leaf && it.member } -> Tree.Leaf(level, true)
                pruned.all { it is Tree.Leaf && !it.member } -> Tree.Leaf(level, false)
                else -> {
                    val prunedRoots = ArrayList<MRoot>(roots.size)
                    val prunedCells = ArrayList<Tree>(cells.size)
                    if (debug) println("Pruned: $pruned")
                    for (i in roots.indices) {
                        val r = roots[i]
                        if (debug) println("($level) after root $r")
                        val a = pruned[i]; val b = pruned[i+1]
                        val necessary = makeLevelList(a.polynomials + b.polynomials, baseLevel = level)[0]
                        if (debug) println("($level) Necessary: $necessary")
                        if (r.poly in necessary || !a.similar(b)) {
                            prunedRoots.add(r); prunedCells.add(a)
                        }
                    }
                    val last = pruned.last()
                    if (prunedCells.isEmpty() && last is Tree.Leaf) {
                        Tree.Leaf(level, last.member)
                    } else {
                        Tree.Cylinder(level, prunedRoots, prunedCells + last)
                    }
                }
            }
        }
    }

    private fun intersectCylinder(rootLevel: Int, a: LevelList, b: LevelList, i: Int = 0): Set<MPoly> {
        return when {
            i >= a.size && i >= b.size -> emptySet()
            i >= a.size -> emptySet()//b[i] // from this level on, there is only b, so b is correct at this point
            i >= b.size -> emptySet()//a[i]
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

    private inline fun List<MRoot>.rootPairs(
        low: MRoot, high: MRoot,
        action: (Pair<MRoot, MRoot>) -> Unit
    ) {
        if (this.isEmpty()) {
            action(low to high)
        } else {
            action(low to first())
            for (i in 0 until (size - 1)) {
                action(this[i] to this[i+1])
            }
            action(last() to high)
        }
    }

    private inline fun List<MRoot>.rootPairs(
        level: Int, action: (Pair<MRoot, MRoot>) -> Unit
    ) {
        val (rootL, rootH) = boundRoots[level]
        val (polyL, polyH) = boundPolynomials[level]
        if (this.isEmpty()) {
            action(MRoot(rootL, 0, polyL) to MRoot (rootH, 0, polyH))
        } else {
            action(MRoot(rootL, 0, polyL) to first())
            for (i in 0 until (size - 1)) {
                action(this[i] to this[i+1])
            }
            action(last() to MRoot (rootH, 0, polyH))
        }
    }

    private fun samplePoint(a: Root, b: Root): NumQ {
        return a.middleValue(b)
    }

    private fun makeLevelList(set: Set<MPoly>, baseLevel: Int = 0, level: Int = dimensions - 1): LevelList {
        if (level == baseLevel) return listOf(set)
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
        return makeLevelList(next, baseLevel,level - 1) + listOf(thisLevel.toSet())
    }

    fun positive(poly: MPoly): Tree {
        val levels = makeLevelList(setOf(poly))
        return makePositive(emptyList(), poly, levels).prune()
    }

    private fun makeLevelListExtended(
        set: Set<MPoly>, level: Int,
        boundPolynomials: List<Pair<MPoly, MPoly>>
    ): LevelList {
        //println("Make levels, starting with $set and doing $level")
        if (level == 0) return listOf(set)
        val carry = set.filter { it.level < level }
        val thisLevel = set - carry
        val list = thisLevel.toList()
        val next = HashSet(carry)
        val (low, high) = boundPolynomials[level]
        for (p in thisLevel) {
            val pu = p.asUnivariate(level)
            next.addAll(UnivariateResultants.Subresultants(pu, pu.derivative()).normalize())
            next.addAll(UnivariateResultants.Subresultants(pu, low.asUnivariate(level)).normalize())
            //println("Projection low $p with $low: ${UnivariateResultants.Subresultants(pu, low.asUnivariate(level)).normalize()}")
            next.addAll(UnivariateResultants.Subresultants(pu, high.asUnivariate(level)).normalize())
            //println("Projection low $p with $high: ${UnivariateResultants.Subresultants(pu, high.asUnivariate(level)).normalize()}")
        }
        for (i in list.indices) {
            for (j in (i+1)..list.lastIndex) {
                next.addAll(UnivariateResultants.Subresultants(list[i].asUnivariate(level), list[j].asUnivariate(level)).normalize())
            }
        }
        return makeLevelListExtended(next, level - 1, boundPolynomials) + listOf(thisLevel.toSet())
    }

    fun positiveExtended(poly: MPoly, extraBounds: Box, extraBoundPoly: List<Pair<MPoly, MPoly>>): Tree {
        val extraDimensions = extraBounds.data.size
        val levels = makeLevelListExtended(setOf(poly),
            level = dimensions + extraDimensions - 1,
            boundPolynomials = extraBoundPoly
        )
        if (debug) println("Extended levels: $levels")
        return makePositiveExtended(
            positive = true,
            point = emptyList(),
            poly = poly,
            levels = levels,
            dimensions = dimensions + extraDimensions,
            boundBox = boundBox.extend(extraBounds),
            boundRoots = extraBoundPoly.map { (a, b) ->
                Root.linear(a.asUnivariate()) to Root.linear(b.asUnivariate())
            }
        ).inRing(ring).also {
            if (debug) println("Unpruned: $it")
        }.prune()
    }

    fun negativeExtended(poly: MPoly, extraBounds: Box, extraBoundPoly: List<Pair<MPoly, MPoly>>): Tree {
        val extraDimensions = extraBounds.data.size
        val levels = makeLevelListExtended(setOf(poly),
            level = dimensions + extraDimensions - 1,
            boundPolynomials = extraBoundPoly
        )
        if (debug) println("Extended levels: $levels")
        return makePositiveExtended(
            positive = false,
            point = emptyList(),
            poly = poly,
            levels = levels,
            dimensions = dimensions + extraDimensions,
            boundBox = boundBox.extend(extraBounds),
            boundRoots = extraBoundPoly.map { (a, b) ->
                Root.linear(a.asUnivariate()) to Root.linear(b.asUnivariate())
            }
        ).inRing(ring).also {
            if (debug) println("Unpruned: $it")
        }.prune()
    }

    private fun makePositiveExtended(
        positive: Boolean,
        point: List<NumQ>,
        poly: MPoly,
        levels: LevelList,
        dimensions: Int,
        boundBox: Box,
        boundRoots: List<Pair<Root, Root>>
    ): Tree {
        val level = point.size
        if (level == dimensions) {
            if (debug) println("Sample point $point and value ${poly.evaluate(*point.toTypedArray())}")
            return Tree.Leaf(level, poly.evaluate(*point.toTypedArray()) > Q.zero == positive)
        }
        val roots = levels[level].findMRoots(point, boundBox.data[level])
        if (debug) println("Roots: $roots in ${boundBox.data[level]}")
        val cells = ArrayList<Tree>(roots.size + 1)
        roots.run {
            val (rootL, rootH) = boundRoots[level]
            if (this.isEmpty()) {
                val sample = samplePoint(rootL, rootH)
                cells.add(makePositiveExtended(positive, point + sample, poly, levels, dimensions, boundBox, boundRoots))
            } else {
                run {
                    val sample = samplePoint(rootL, first().root)
                    cells.add(makePositiveExtended(positive, point + sample, poly, levels, dimensions, boundBox, boundRoots))
                }
                for (i in 0 until (size - 1)) {
                    val sample = samplePoint(this[i].root, this[i+1].root)
                    cells.add(makePositiveExtended(positive, point + sample, poly, levels, dimensions, boundBox, boundRoots))
                }
                run {
                    val sample = samplePoint(last().root, rootH)
                    cells.add(makePositiveExtended(positive, point + sample, poly, levels, dimensions, boundBox, boundRoots))
                }
            }
        }
        return if (level >= this.dimensions) {
            Tree.Leaf(level, cells.any { it.any() })
        } else {
            Tree.Cylinder(level, roots, cells)
        }
    }

    sealed class Tree {

        abstract val level: Int
        abstract val roots: List<MRoot>
        abstract val levelList: List<Set<MPoly>>
        abstract val polynomials: Set<MPoly>

        abstract fun any(): Boolean
        abstract fun all(): Boolean

        //abstract fun lookup(point: NumQ): Tree
        abstract fun lookup(cell: Int): Tree
        abstract fun similar(other: Tree): Boolean
        abstract fun inRing(ring: IPolynomialRing<MPoly>): Tree

        data class Leaf(override val level: Int, val member: Boolean) : Tree() {

            override val levelList: List<Set<MPoly>>
                get() = emptyList()

            override val roots: List<MRoot>
                get() = emptyList()

            override val polynomials: Set<MPoly>
                get() = emptySet()

            //override fun lookup(point: NumQ): Tree = this
            override fun lookup(cell: Int): Tree = this
            override fun similar(other: Tree): Boolean = this == other
            override fun any(): Boolean = member
            override fun all(): Boolean = member
            override fun inRing(ring: IPolynomialRing<MPoly>): Tree = this

        }

        data class Cylinder(
            override val level: Int,
            override val roots: List<MRoot>,
            val cells: List<Tree>
        ) : Tree() {

            override val levelList: List<Set<MPoly>> =
                listOf<Set<MPoly>>(roots.mapTo(HashSet()) { it.poly }) + cells
                    .map { it.levelList }
                    .fold1 { a, b -> a.zip(b) }

            override val polynomials: Set<MPoly> = roots.mapTo(HashSet()) { it.poly } + cells
                .map { it.polynomials }.fold1 { a, b -> a + b }

            override fun inRing(ring: IPolynomialRing<MPoly>): Tree {
                if (roots.isEmpty()) return Cylinder(level, roots, cells.map { it.inRing(ring) })
                val eliminate = roots.first().poly.nVariables - ring.nVariables()
                val array = IntArray(eliminate) { ring.nVariables() + it }
                return Cylinder(level, roots.map { (root, n, poly) ->
                    MRoot(root, n, poly.eliminate(array, Array(eliminate) { Q.zero }))
                }, cells.map { it.inRing(ring) })
            }

            /*override fun lookup(point: NumQ): Tree {
                val searchResult = roots.binarySearch { (root, _) -> root.compareTo(point) }
                if (searchResult >= 0) {
                    error("The point $point is also a root in $roots. This should not happen!")
                } else {
                    return cells[(-searchResult - 1)]
                }
            }*/

            override fun lookup(cell: Int): Tree = cells[cell]

            override fun similar(other: Tree): Boolean {
                if (other !is Cylinder) return false
                if (this.cells.size != other.cells.size) return false
                return this.cells.asSequence().zip(other.cells.asSequence()).all { (a, b) -> a == b }
            }

            override fun any(): Boolean = cells.any { it.any() }
            override fun all(): Boolean = cells.all { it.all() }

        }

    }

}

fun Iterable<MRoot>.moveTo(point: List<NumQ>, bounds: Interval): List<MRoot> {
    val set = this.toSet()
    return this.map { it.poly }.findMRoots(point, bounds).filter { it in set }.sortedBy { it.root }
}

fun Iterable<MPoly>.findMRoots(point: List<NumQ>, bounds: Interval): List<MRoot> = this.flatMap { it.findMRoots(point, bounds) }.sortedBy { it.root }
fun Iterable<MPoly>.findMRoots(point: List<NumQ>, bounds: Interval, l: Root, h: Root): List<MRoot> = this
    .flatMap { it.findMRoots(point, bounds) }
    .filter { l < it.root && it.root < h }
    .sortedBy { it.root }

fun MPoly.findMRoots(point: List<NumQ>, bounds: Interval): List<MRoot> = AdaptiveRootIsolation.isolateRootsInBounds(listOf(this.evaluate(point).asUnivariate()), bounds)
    .sorted().mapIndexed { index, root -> MRoot(root, index, this) }

data class MRoot(
    val root: Root, val n: Int, val poly: MPoly
) {

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as MRoot

        if (n != other.n) return false
        if (poly != other.poly) return false

        return true
    }

    override fun hashCode(): Int {
        var result = n
        result = 31 * result + poly.hashCode()
        return result
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
        val circPoly = ring.parse("(x-8/10)^2 + (y-1/2)^2 - 1/10")
        val set1 = positive(poly).not()
        val set2 = positive(poly2).not()
        val circ = positive(circPoly).not()
        println("Set: $set1")
        //println("Simplified: ${set1.prune()}")
        println("Set2: $set2")
        //println("Simplified: ${set2.prune()}")
        val or = set1 and set2
        println("Intersection: $or")
        println("Circle: $circ")
        println("Is subset: ${circ subset or}")
    }
}