package biodivine.algebra.rootisolation

import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.univar.UnivariateFactorization
import cc.redberry.rings.poly.univar.UnivariatePolynomial
import biodivine.algebra.getDefaultBoundForDescartMethod
import biodivine.algebra.getNumberOfSignChanges
import biodivine.algebra.ia.Interval
import biodivine.algebra.transformPolyToInterval

object DescartWithPolyFactorisationRootIsolation : RootIsolation {

    override fun isolateRoots(polynomial: UnivariatePolynomial<NumQ>, precision: NumQ): List<Interval> {
        val computedRoots = mutableListOf<Interval>()
        for (poly in UnivariateFactorization.Factor(polynomial)) {
            if (poly.isLinearExactly) {
                computedRoots += Interval(poly[0].negate(), poly[0].negate())
            } else {
                computedRoots += isolateIrationalRoots(poly, precision)
            }
        }
        return computedRoots
    }


    private fun isolateIrationalRoots(polynomial: UnivariatePolynomial<NumQ>, precision: NumQ): List<Interval> {
        val defaultBound = polynomial.getDefaultBoundForDescartMethod()
        return isolateIrracionalRootsRecursively(polynomial, defaultBound.negate(), defaultBound, precision)
    }

    private fun isolateIrracionalRootsRecursively(polynomial: UnivariatePolynomial<NumQ>, lowerBound: NumQ, upperBound: NumQ, precision: NumQ): List<Interval> {
        val numberOfSignChangesInInterval = polynomial.transformPolyToInterval(lowerBound, upperBound).getNumberOfSignChanges()
        val middleValue = (upperBound.add(lowerBound)).divide(BigInteger.TWO)
        val actualPrecision = (upperBound.subtract(lowerBound)).divide(BigInteger.TWO)

        if (numberOfSignChangesInInterval == 0) {
            return listOf()
        }

        if (numberOfSignChangesInInterval == 1 && precision > actualPrecision) {
            return listOf(Interval(lowerBound, upperBound))
        }

        return isolateIrracionalRootsRecursively(polynomial, lowerBound, middleValue, precision) +
                isolateIrracionalRootsRecursively(polynomial, middleValue, upperBound, precision)
    }

}