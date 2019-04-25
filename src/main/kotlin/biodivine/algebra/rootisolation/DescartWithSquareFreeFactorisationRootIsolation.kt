package biodivine.algebra.rootisolation

import biodivine.algebra.ia.Interval
import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.univar.UnivariatePolynomial
import cc.redberry.rings.poly.univar.UnivariateSquareFreeFactorization
import biodivine.algebra.getDefaultBoundForDescartMethod
import biodivine.algebra.getNumberOfSignChanges
import biodivine.algebra.transformPolyToInterval

object DescartWithSquareFreeFactorisationRootIsolation : RootIsolation {

    override fun isolateRoots(polynomial: UnivariatePolynomial<NumQ>, precision: NumQ): List<Interval> {
        val computedRoots = mutableListOf<Interval>()
        for (poly in UnivariateSquareFreeFactorization.SquareFreeFactorization(polynomial)) {
            computedRoots += isolateAllRoots(poly, precision)
        }
        return computedRoots
    }

    private fun isolateAllRoots(polynomial: UnivariatePolynomial<NumQ>, precision: NumQ): List<Interval> {
        val defaultBound = polynomial.getDefaultBoundForDescartMethod()
        val boundRoots = listOf(defaultBound, defaultBound.negate()).filter { polynomial.evaluate(it).isZero }.map { q -> Interval(q,q) }
        return boundRoots + isolateAllRootsRecursively(polynomial, defaultBound.negate(), defaultBound, precision)
    }


    private fun isolateAllRootsRecursively(polynomial: UnivariatePolynomial<NumQ>, lowerBound: NumQ, upperBound: NumQ, precision: NumQ): List<Interval> {
        val numberOfSignChangesInInterval = polynomial.transformPolyToInterval(lowerBound, upperBound).getNumberOfSignChanges()
        val middleValue = (upperBound.add(lowerBound)).divide(BigInteger.TWO)
        val actualPrecision = (upperBound.subtract(lowerBound)).divide(BigInteger.TWO)

        if (numberOfSignChangesInInterval == 0) {
            return listOf()
        }

        if (numberOfSignChangesInInterval == 1 && precision > actualPrecision) {
            return listOf(Interval(lowerBound, upperBound))
        }

        return isolateAllRootsRecursively(polynomial, lowerBound, middleValue, precision) +
                isolateAllRootsRecursively(polynomial, middleValue, upperBound, precision) +
                listOf(middleValue).filter { polynomial.evaluate(it).isZero }.map { q ->
                    Interval(q, q)
                }
    }
}