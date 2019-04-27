package biodivine.algebra.rootisolation

import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.univar.UnivariateFactorization
import cc.redberry.rings.poly.univar.UnivariatePolynomial
import biodivine.algebra.getDefaultBoundForDescartMethod
import biodivine.algebra.getNumberOfSignChanges
import biodivine.algebra.ia.Interval
import biodivine.algebra.transformPolyToInterval
import cc.redberry.rings.poly.PolynomialMethods.Factor

object DescartWithPolyFactorisationRootIsolation : RootIsolation {

    override fun isolateInBounds(
        polynomial: UnivariatePolynomial<NumQ>,
        bounds: Interval,
        precision: NumQ
    ): List<Interval> {
        val computedRoots = mutableListOf<Interval>()
        for (poly in UnivariateFactorization.Factor(polynomial)) {
            if (poly.isLinearExactly) {
                val root = poly[0].negate().divide(poly[1])
                //println("Poly $poly is linear and root is $root")
                if (root in bounds) {
                    computedRoots += Interval(root, root)
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

            if (polynomial.evaluate(middleValue).isZero) {
                result.add(Interval(middleValue, middleValue))
            }

            when {
                numberOfSignChangesInInterval == 0 -> Unit
                numberOfSignChangesInInterval == 1 && precision > actualPrecision && lowerBound != bounds.low && upperBound != bounds.high -> {
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

    override fun isolateRootsInBounds(
        polynomials: Collection<UnivariatePolynomial<NumQ>>,
        bounds: Interval,
        precision: NumQ
    ): List<Interval> {
        val factors = polynomials.flatMapTo(HashSet()) { Factor(it) }
        val linearRoots = factors.filter { it.isLinearExactly }.mapNotNull { poly ->
            val root = poly[0].negate().divide(poly[1])
            if (root !in bounds) null else Interval(root, root)
        }
        //val nonLinearCombination = factors.filter { !it.isLinearExactly }
            //.fold(UnivariatePolynomial.one(polynomials.first().ring)) { a, b -> a.multiply(b) }
        val nonLinearRoots = factors.filter { !it.isLinearExactly }.flatMap {
            isolateIrracionalRootsRecursively(it, bounds, precision)
        }.toHashSet()
        return linearRoots + nonLinearRoots
    }
}