package biodivine.algebra.rootisolation

import biodivine.algebra.ia.Interval
import cc.redberry.rings.Rational
import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.univar.UnivariatePolynomial

typealias NumQ = Rational<BigInteger>

interface RootIsolation {
    fun isolateInBounds(polynomial: UnivariatePolynomial<NumQ>, bounds: Interval, precision: NumQ): List<Interval>
    fun isolateRoots(polynomial: UnivariatePolynomial<NumQ>, precision: NumQ): List<Interval>

    fun isolateRootsInBounds(polynomials: Collection<UnivariatePolynomial<NumQ>>, bounds: Interval, precision: NumQ): List<Interval> {
        return polynomials.flatMap { isolateInBounds(it, bounds, precision) }.toSet().toList()
    }

    fun isolateRoots(polynomials: Collection<UnivariatePolynomial<NumQ>>, precision: NumQ): List<Interval> {
        val result = mutableListOf<Interval>()
        for (poly in polynomials) {
            result.addAll(isolateRoots(poly, precision))
        }
        return result
    }
}

