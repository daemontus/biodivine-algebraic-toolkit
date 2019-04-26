package biodivine.algebra

import cc.redberry.rings.Rings
import cc.redberry.rings.Rings.Q
import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.univar.UnivariatePolynomial

val coder = Rings.UnivariateRingQ.mkCoder("x")

private val onePlusX = (0..10).map { coder.parse("(1 + x)^$it") }

/***
 * transform polynomial p by formula : (1 + x)^n * p( (ax + b) / (1 + x))
 */

fun UnivariatePolynomial<NumQ>.transformPolyToInterval(lowerBound: NumQ, upperBound: NumQ): UnivariatePolynomial<NumQ> {
    val result = UnivariatePolynomial.zero(Rings.Q)
    var exponent = degree()

    for (coef in this) {
        val multiplyPoly = onePlusX[exponent].copy()
        val subst = coder.parse("($lowerBound * x + $upperBound)^${degree() - exponent}")
        result.add(multiplyPoly.multiply(subst.multiply(coef)))
        exponent -= 1
    }
    return result
}

/**
 *  get number of sign changes in coefficients of input polynomial
 */

fun UnivariatePolynomial<NumQ>.getNumberOfSignChanges(): Int {
    var numberOfChanges = 0
    var previousCoef = this.getFirstNonZeroCoefficient()
    for (coef in this) {
        if (!coef.isZero) {
            if (isChangeInSign(previousCoef, coef)) {
                numberOfChanges++
            }
            previousCoef = coef
        }

    }
    return numberOfChanges
}

private fun isChangeInSign(first: NumQ, second: NumQ): Boolean {
    return first.signum() != second.signum()
}

fun UnivariatePolynomial<NumQ>.getFirstNonZeroCoefficient(): NumQ {
    for (coef in this) {
        if (!coef.isZero) {
            return coef
        }
    }
    return Rings.Q.zero
}

fun UnivariatePolynomial<NumQ>.getDefaultBoundForDescartMethod(): NumQ {
    if (this.isZero)
        return Q.zero

    var bound = Q.zero
    for (coefficient in this) {
        val fraction = (coefficient.divide(lc())).abs()
        if (fraction > bound) {
            bound = fraction
        }
    }
    return bound.multiply(BigInteger.TWO)
}

fun <T> UnivariatePolynomial<T>.reductum(): UnivariatePolynomial<T> {
    return copy().shiftLeft(1)
}

fun <T> UnivariatePolynomial<T>.reductaSet(): Set<UnivariatePolynomial<T>> {
    if (isZero)
        return emptySet()

    val result = HashSet<UnivariatePolynomial<T>>(this.degree())
    var poly = this
    while (!poly.isZero) {
        result.add(poly)
        poly = poly.reductum()
    }
    return result
}




