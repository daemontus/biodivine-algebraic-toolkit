package biodivine.algebra.params

import biodivine.algebra.MPoly
import biodivine.algebra.ia.Interval
import biodivine.algebra.synth.Box
import cc.redberry.rings.Rings
import cc.redberry.rings.poly.IPolynomialRing
import java.lang.IllegalStateException
import java.util.*
import kotlin.collections.HashSet
import kotlin.system.measureTimeMillis

data class Cell(
    private val coordinates: IntArray
) {

    fun project(retain: Int) = Cell(coordinates.take(retain).toIntArray())

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as Cell

        if (!coordinates.contentEquals(other.coordinates)) return false

        return true
    }

    override fun hashCode(): Int {
        return coordinates.contentHashCode()
    }

    override fun toString(): String {
        return "Cell(coordinates=${Arrays.toString(coordinates)})"
    }


}

/**
 * Represents a semi-algebraic set defined by the given set of algebraic varieties (0 = polynomial).
 * The valid cells of the decomposition are then given in the bit set.
 *
 */
data class SemiAlgSet(
    val levelGraph: LevelGraph,
    val validCells: Set<Cell>
)

class SemiAlgSolver(
    private val bounds: Box,
    private val ring: IPolynomialRing<MPoly>
) {

    private val dimensions = bounds.data.size

    val zero: SemiAlgSet = SemiAlgSet(LevelGraph(emptyList(), ring, bounds), emptySet())
    // if there are no bounds, there is only one cell: the whole box - all coordinates zero
    val one: SemiAlgSet = SemiAlgSet(LevelGraph(emptyList(), ring, bounds), setOf(Cell(IntArray(dimensions)))).also { println("One: $it") }

    /*
        General overview of the operations:

        1.  Create a "common" CADPolynomials - a CADPolynomials of the union of the two bound polynomial sets.
        2.  Translate each cell set into this common CADPolynomials domain - i.e. for each cell/sample point in the new cad,
            check if it belongs to a cell in the old cad.
        3.  Perform the operation cell-wise, since now we have a common representation for both sets.
        4.  Remove redundant polynomials. How do you do this? Well, there might be other, more efficient options,
        but current idea is like this: We try to remove a polynomial from the CADPolynomials, we then compute the new CADPolynomials and
        for each cell in the original CADPolynomials, we compute the ID in the new CADPolynomials - if two incompatible cells are merged,
        we know we removed a wrong polynomial. All polynomials are independent in this sense - but remember that
        there might be 3 polynomials which are redundant, but by removing first two, the third is no longer redundant
        (because they have common roots!). So always check redundancy of other polynomials on the already reduced CADPolynomials.
     */

    infix fun SemiAlgSet.subset(that: SemiAlgSet): Boolean {
        val levelUnion = LevelGraph(this.levelGraph, that.levelGraph)
        return levelUnion.walkCells().all { (point, _) ->
            val cellInThis = this.levelGraph.cellForPoint(point)
            if (cellInThis !in this.validCells) true else {
                val cellInThat = that.levelGraph.cellForPoint(point)
                cellInThat in that.validCells
            }
        }
    }

    infix fun SemiAlgSet.or(that: SemiAlgSet): SemiAlgSet {
        val levelUnion = LevelGraph(this.levelGraph, that.levelGraph)
        val validDisjunction = levelUnion.walkCells().mapNotNull { (point, cell) ->
            val cellInThis = this.levelGraph.cellForPoint(point)
            if (cellInThis in this.validCells) cell else {
                val cellInThat = that.levelGraph.cellForPoint(point)
                if (cellInThat in that.validCells) cell else null
            }
        }.toSet()
        return SemiAlgSet(levelUnion, validDisjunction).simplify()
    }

    infix fun SemiAlgSet.and(that: SemiAlgSet): SemiAlgSet {
        if (this.validCells.isEmpty() || that.validCells.isEmpty()) return zero
        val levelUnion = LevelGraph(this.levelGraph, that.levelGraph)
        val validIntersection = levelUnion.walkCells().mapNotNull { (point, cell) ->
            val cellInThis = this.levelGraph.cellForPoint(point)
            if (cellInThis !in this.validCells) null else {
                val cellInThat = that.levelGraph.cellForPoint(point)
                if (cellInThat !in that.validCells) null else cell
            }
        }.toSet()
        return SemiAlgSet(levelUnion, validIntersection).simplify()
    }

    fun SemiAlgSet.not(): SemiAlgSet {
        val negated = levelGraph.walkCells().mapNotNull { (_, cell) -> cell.takeIf { it !in validCells } }.toSet()
        return SemiAlgSet(levelGraph, negated)
    }

    fun SemiAlgSet.simplify(): SemiAlgSet {
        var simplified = this
        val checkCells = levelGraph.walkCells()
        poly@ for (poly in this.levelGraph.basis) {
            val removed = simplified.levelGraph - poly
            val valid = HashSet<Cell>()
            val invalid = HashSet<Cell>()
            for ((point, cell) in checkCells) {
                val cellInRemoved = removed.cellForPoint(point)
                if (cell in validCells) {
                    // cell should be valid - if we found it in invalid before, this simplification is not safe
                    if (cellInRemoved in invalid) continue@poly
                    valid.add(cellInRemoved)
                } else {
                    // cell should be invalid - if we found it in valid before, this simplification is not safe
                    if (cellInRemoved in valid) continue@poly
                    invalid.add(cellInRemoved)
                }
            }
            simplified = SemiAlgSet(removed, valid)
        }
        return simplified
    }

    /**
     * Create a projection of the semi-algebraic set, retaining only the first [retain] variables.
     */
    fun SemiAlgSet.project(retain: Int, reducedRing: IPolynomialRing<MPoly>? = null): SemiAlgSet {
        val levelProjection = LevelGraph(this.levelGraph, retain, reducedRing)
        val validityProjection = this.validCells.mapTo(HashSet()) { cell ->
            cell.project(retain)
        }
        return SemiAlgSet(levelProjection, validityProjection)
    }

    fun SemiAlgSet.isEmpty(): Boolean = this.validCells.isEmpty()
    fun SemiAlgSet.isNotEmpty(): Boolean = this.validCells.isNotEmpty()

    fun positive(poly: MPoly): SemiAlgSet {
        val levels = LevelGraph(listOf(poly), ring, bounds)
        val valid = levels.walkCells().mapNotNull { (point, cell) ->
            cell.takeIf { poly.evaluate(*point.toTypedArray()) > Rings.Q.zero }
        }.toSet()
        return SemiAlgSet(levels, valid)
    }

    fun negative(poly: MPoly): SemiAlgSet {
        val levels = LevelGraph(listOf(poly), ring, bounds)
        val valid = levels.walkCells().mapNotNull { (point, cell) ->
            cell.takeIf { poly.evaluate(*point.toTypedArray()) < Rings.Q.zero }
        }.toSet()
        return SemiAlgSet(levels, valid)
    }

}

fun main() {
    val ring = Rings.MultivariateRingQ(2)
    val bounds = Box(Interval(0, 2), Interval(0, 2))
    val solver = SemiAlgSolver(bounds, ring)
    val p = ring.parse("x^2-y")
    val q = ring.parse("x^2 - 4*x + 4 - y")
    solver.run {
        val pNeg = solver.negative(p)
        println("P: $pNeg")
        val qNeg = solver.negative(q)
        println("Q: $qNeg")
        val a = pNeg and qNeg
        val b = pNeg or qNeg
        val c = pNeg subset qNeg
        val d = pNeg and qNeg.not()
        val e = (pNeg and qNeg.not()) subset pNeg
        println("P and Q: $a")
        println("P or Q: $b")
        println("P subset Q: $c")
        println("P and !Q: $d")
        println("P and !Q subset P: $e")
    }
    val elapsed = measureTimeMillis {
        var v = true
        repeat(2500) { i ->
            val pNeg = solver.negative(p)
            val qNeg = solver.negative(q)
            solver.run {
                val a = pNeg and qNeg
                val b = pNeg or qNeg
                val c = pNeg subset qNeg
                val d = pNeg and qNeg.not()
                val e = (pNeg and qNeg.not()) subset pNeg
                v = v && a.isNotEmpty() && b.isNotEmpty() || c && d.isEmpty() || e
                if (i % 100 == 0) {
                    println("Iter: $i (${(a and b and d).isEmpty()} $c $e")
                }
            }
        }
    }
    println("Elapsed: $elapsed")
}