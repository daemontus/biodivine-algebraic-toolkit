package biodivine.algebra.rootisolation

import biodivine.algebra.LRUCache
import biodivine.algebra.UPoly
import biodivine.algebra.getNumberOfSignChanges
import biodivine.algebra.ia.Interval
import biodivine.algebra.ia.div
import biodivine.algebra.ia.minus
import biodivine.algebra.ia.plus
import biodivine.algebra.transformPolyToInterval
import cc.redberry.rings.Rings
import cc.redberry.rings.poly.univar.UnivariateSquareFreeFactorization
import java.util.*
import kotlin.collections.ArrayList
import kotlin.collections.HashSet

object AdaptiveRootIsolation {

    private val rootIsolationCache = ThreadLocal.withInitial<LRUCache<UPoly, List<Root>>> {
        LRUCache(
            10_000
        )
    }

    fun isolateRootsInBounds(polynomials: Collection<UPoly>, bounds: Interval): NavigableSet<Root> {
        val result = TreeSet<Root>()    // results are saved
        polynomials
            .flatMapTo(HashSet()) { listOf(it)/*UnivariateSquareFreeFactorization.SquareFreeFactorization(it)*/ }    // first, find all unique factors
            .forEach { factor ->                                            // then, find roots for every factor
                if (factor.isLinearExactly) {
                    val root = Root.rational(factor)
                    if (root in bounds) result.add(root)
                } else {
                    // due to initial interval, we know all roots are in bounds
                    val cache = rootIsolationCache.get()
                    val cached = cache.get(factor)
                    if (cached != null) result.addAll(cached) else {
                        val roots = isolateIrrationalRoots(factor, bounds)
                        cache.set(factor, roots)
                        result.addAll(roots)
                    }
                }
            }
        return result
    }

    private fun isolateIrrationalRoots(polynomial: UPoly, initialInterval: Interval): List<Root> {
        // At this point, we know there are no rational roots and that all roots we create are actually
        // unique, but only for this polynomial!
        val roots = ArrayList<Root>()
        val queue = ArrayList<Interval>()
        queue.add(initialInterval)
        while (queue.isNotEmpty()) {
            val (low, high) = queue.removeAt(queue.lastIndex)
            val signChanges = polynomial.transformPolyToInterval(low, high).getNumberOfSignChanges()

            when {
                signChanges == 1 -> roots.add(Root.irrational(low, high, polynomial))
                signChanges > 1 -> {
                    val middle = low + (high - low) / 2

                    if (middle.isRationalRoot(polynomial)) {
                        roots.add(Root.irrational(middle, middle, polynomial))
                    }

                    queue.add(Interval(low, middle))
                    queue.add(Interval(middle, high))
                }
            }
        }
        val (low, high) = initialInterval
        if (low.isRationalRoot(polynomial)) roots.add(Root.irrational(low, low, polynomial))
        if (high.isRationalRoot(polynomial)) roots.add(Root.irrational(high, high, polynomial))
        return roots
    }

    private fun NumQ.isRationalRoot(poly: UPoly): Boolean = poly.evaluate(this) == Rings.Q.zero

}