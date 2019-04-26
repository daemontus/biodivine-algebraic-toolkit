package biodivine.algebra.rootisolation

import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.univar.UnivariateFactorization
import cc.redberry.rings.poly.univar.UnivariatePolynomial
import biodivine.algebra.getDefaultBoundForDescartMethod
import biodivine.algebra.getNumberOfSignChanges
import biodivine.algebra.ia.Interval
import biodivine.algebra.transformPolyToInterval

object DescartWithPolyFactorisationRootIsolation : RootIsolation {

    override fun isolateInBounds(
        polynomial: UnivariatePolynomial<NumQ>,
        bounds: Interval,
        precision: NumQ
    ): List<Interval> {
        val computedRoots = mutableListOf<Interval>()
        for (poly in UnivariateFactorization.Factor(polynomial)) {
            if (poly.isLinearExactly) {
                if (poly[0].negate() in bounds) {
                    computedRoots += Interval(poly[0].negate(), poly[0].negate())
                }
            } else {
                computedRoots += isolateIrracionalRootsRecursively(poly, bounds, precision)
            }
        }
        return computedRoots
    }

    override fun isolateRoots(polynomial: UnivariatePolynomial<NumQ>, precision: NumQ): List<Interval> {
        val defaultBound = Interval(polynomial.getDefaultBoundForDescartMethod().negate(), polynomial.getDefaultBoundForDescartMethod())
        return isolateInBounds(polynomial, defaultBound, precision)
    }

    private fun isolateIrracionalRootsRecursively(polynomial: UnivariatePolynomial<NumQ>, bounds: Interval, precision: NumQ): List<Interval> {
        val result = ArrayList<Interval>()
        val workQueue = ArrayList<Interval>()
        workQueue.add(bounds)
        while (workQueue.isNotEmpty()) {
            val (lowerBound, upperBound) = workQueue.removeAt(workQueue.lastIndex)
            val numberOfSignChangesInInterval = polynomial.transformPolyToInterval(lowerBound, upperBound).getNumberOfSignChanges()
            val middleValue = (upperBound.add(lowerBound)).divide(BigInteger.TWO)
            val actualPrecision = (upperBound.subtract(lowerBound)).divide(BigInteger.TWO)

            when {
                numberOfSignChangesInInterval == 0 -> Unit
                numberOfSignChangesInInterval == 1 && precision > actualPrecision -> {
                    result.add(Interval(lowerBound, upperBound))
                }
                else -> {
                    workQueue.add(Interval(lowerBound, middleValue))
                    workQueue.add(Interval(middleValue, upperBound))
                }
            }
        }
        return result
    }

}