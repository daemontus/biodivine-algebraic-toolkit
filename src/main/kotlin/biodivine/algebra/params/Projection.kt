package biodivine.algebra.params

import biodivine.algebra.LRUCache
import biodivine.algebra.MPoly
import biodivine.algebra.synth.cache_size
import cc.redberry.rings.poly.univar.UnivariateResultants

object Projection {

    private val discriminantCache = ThreadLocal.withInitial<LRUCache<Pair<MPoly, Int>, List<MPoly>>> {
        LRUCache(
            cache_size
        )
    }

    private val resultantCache = ThreadLocal.withInitial<LRUCache<CacheKey, List<MPoly>>> {
        LRUCache(
            cache_size
        )
    }

    fun discriminant(poly: MPoly, variable: Int): List<MPoly> {
        val cache = discriminantCache.get()
        val key = poly to variable
        val cached = cache.get(key)
        return if (cached != null) cached else {
            val u = poly.asUnivariate(variable)
            val projection = UnivariateResultants.Subresultants(u, u.derivative()).normalize()
            cache.set(key, projection)
            projection
        }
    }

    fun resultant(a: MPoly, b: MPoly, variable: Int): List<MPoly> {
        val cache = resultantCache.get()
        val key = CacheKey(a, b, variable)
        val cached = cache.get(key)
        return if (cached != null) cached else {
            val projection = UnivariateResultants.Subresultants(a.asUnivariate(variable), b.asUnivariate(variable)).normalize()
            cache.set(key, projection)
            projection
        }
    }

    data class CacheKey(
        val a: MPoly, val b: MPoly, val variable: Int
    ) {

        /*TODO
        override fun equals(other: Any?): Boolean {
            if (this === other) return true
            if (javaClass != other?.javaClass) return false

            other as CacheKey

            // Importnat: keys are equivalent if they are symmetric
            if (a != other.a && a != other.b) return false
            if (b != other.b && b != other.a) return false
            if (a == b && other.a != other.b) return false
            if (variable != other.variable) return false

            return true
        }

        override fun hashCode(): Int {
            var result = a.hashCode() * b.hashCode()    // symmetric operation
            result = 31 * result + variable
            return result
        }*/
    }

}