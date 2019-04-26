package biodivine.algebra.rootisolation

import biodivine.algebra.ia.Interval
import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.univar.UnivariatePolynomial
import cc.redberry.rings.poly.univar.UnivariateSquareFreeFactorization
import biodivine.algebra.getDefaultBoundForDescartMethod
import biodivine.algebra.getNumberOfSignChanges
import biodivine.algebra.ia.evaluate
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


    private fun isolateAllRootsRecursively(polynomial: UnivariatePolynomial<NumQ>, mainLowerBound: NumQ, mainUpperBound: NumQ, precision: NumQ): List<Interval> {
        val result = ArrayList<Interval>()
        val workQueue = ArrayList<Interval>()
        workQueue.add(Interval(mainLowerBound, mainUpperBound))
        while (workQueue.isNotEmpty()) {
            val interval = workQueue.removeAt(workQueue.lastIndex)
            if (!polynomial.evaluate(interval).hasZero) continue    // fast-skip non-zero intervals
            val (lowerBound, upperBound) = interval
            val numberOfSignChangesInInterval = polynomial.transformPolyToInterval(lowerBound, upperBound).getNumberOfSignChanges()
            val middleValue = (upperBound.add(lowerBound)).divide(BigInteger.TWO)
            val actualPrecision = (upperBound.subtract(lowerBound)).divide(BigInteger.TWO)

            if (polynomial.evaluate(middleValue).isZero) {
                result.add(Interval(middleValue, middleValue))
            }

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