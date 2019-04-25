package biodivine.algebra

import cc.redberry.rings.Rational
import cc.redberry.rings.Rings.*
import cc.redberry.rings.bigint.BigInteger
import cc.redberry.rings.poly.PolynomialMethods.Factor
import cc.redberry.rings.poly.multivar.GroebnerBases
import cc.redberry.rings.poly.multivar.MonomialOrder
import cc.redberry.rings.poly.multivar.MultivariatePolynomial
import cc.redberry.rings.poly.multivar.MultivariateResultants
import cc.redberry.rings.poly.univar.UnivariatePolynomial
import cc.redberry.rings.poly.univar.UnivariateResultants

typealias NumQ = Rational<BigInteger>
typealias MPoly = MultivariatePolynomial<NumQ>
typealias UPoly = UnivariatePolynomial<NumQ>
typealias UMPoly = UnivariatePolynomial<MPoly>

fun main(args: Array<String>) {
    val vars = arrayOf("x", "y")    // variable names are inferred automatically, but lets make sure
    val ring = MultivariateRingQ(2)
    val coder = ring.mkCoder(*vars)

    val p1 = coder.parse("-x + y")
    val p2 = coder.parse("-x - y")
    //val pInit = coder.parse("(x+3)^2 + (y-3)^2 - 1/4")
    val pFinal = coder.parse("x^2 + y^2 - 1")

    val system = arrayOf(p1, p2)
    //val basePolynomials = setOf(p1, p2, pInit, pFinal)
    //val boundPolynomials = basePolynomials.map { it.boundary(system) }.toSet()
    //val allPolynomials = basePolynomials + boundPolynomials

    /*println("All polynomials: ${allPolynomials.size}")
    println("Base polynomials: ${basePolynomials.size}")
    basePolynomials.forEach { println(it) }
    println("Bound polynomials: ${boundPolynomials.size}")
    boundPolynomials.forEach { println(it) }

    //val equilibriumPoly = (p1.copy().multiply(p1).add(p2.copy().multiply(p2)).subtract(ring.parse("1/100")))
    //println("EQ poly: $equilibriumPoly")

    val b1 = coder.parse("12*y - 2*y^2 - 2*x^2")
    val b2 = coder.parse("x^2 + y^2 - 1")
    val q = coder.parse("5041-3456*y+576*y^2")
    println("Q bound: ${Factor(q.boundary(system)).toList()}")
    val order = MonomialOrder.EliminationOrder(MonomialOrder.GREVLEX, 1)
    val bases = GroebnerBases.GroebnerBasis(listOf(pInit), order)
    bases.forEach { println("Bases: $it") }
    */

    //project(allPolynomials.map { it.asUnivariate(1) }.toSet())

    /*val eSystem = system.map { it.insertVariable(2) }.toTypedArray()
    val f1 = pFinal.insertVariable(2)
    val f2 = p1.insertVariable(2)

    val eRing = MultivariateRingQ(3)
    val bound = f1.copy().add(f2.copy().multiply(eRing.parse("z")))
    val border = bound.boundary(eSystem)

    println("F1: $f1")
    println("F2: $f2")
    println("Bound: $bound")
    println("Bound border: $border")

    println("REs: ${UnivariateResultants.Subresultants(bound.asUnivariate(0), border.asUnivariate(0))}")

    val z = eRing.parse("1+1*y*z^3+1*y^2*z^2")
    println("Res: ${UnivariateResultants.Subresultants(z.asUnivariate(1), z.asUnivariate(1).derivative()).map {
        Factor(it)
    }}")
    val proj = project(listOf(f1, f2, bound, border).map { it.asUnivariate(0) }.toSet())
    val proj2 = project(proj.map { it.asUnivariate(1) }.toSet())
    */

    val eSystem = system.map { it.insertVariable(2).insertVariable(3) }.toTypedArray()
    val f1 = pFinal.insertVariable(2).insertVariable(3)
    val f2 = p1.insertVariable(2).insertVariable(3)

    val eRing = MultivariateRingQ(4)
    // f1 + a(f2 + b*f1*f2))
    val bound = f1.copy().add(
        eRing.parse("x3").multiply(
            f2.copy().add(
                eRing.parse("x4").multiply(f1).multiply(f2)
            )
        )
    )
    val border = bound.boundary(eSystem)

    println("Bound: $bound")
    println("Bound border: $border")

    val proj1 = project(listOf(f1, f2, bound, border).map { it.asUnivariate(0) }.toSet())
    println("Projection 1: ${proj1.size}")
    proj1.forEach { println(it) }
    val proj2 = project(proj1.map { it.asUnivariate(1) }.toSet())

    println("Projection 2: ${proj2.size}")
    proj2.forEach { println(it) }

    val proje3 = project(proj2.map { it.asUnivariate(2) }.toSet())
    println("Projection 3: ${proje3.size}")
    proje3.forEach { println(it) }
}

/*
fun <T> UnivariatePolynomial<T>.reductum(): UnivariatePolynomial<T> =
    if (isZero) this else copy().shiftLeft(1)
 */

fun <T> UnivariatePolynomial<T>.reductumSet(): Set<UnivariatePolynomial<T>> {
    val result = HashSet<UnivariatePolynomial<T>>(this.degree())
    var poly = this
    while (!poly.isZero) {
        result.add(poly)
        poly = poly.reductum()
    }
    return result
}

fun project(polynomials: Set<UMPoly>): Set<MPoly> {
    // compute first projection set using partial derivatives
    val projection1: Set<MPoly> = polynomials
        .mapIndexed { i, polynomial -> i to polynomial }
        .flatMap { (i, p) ->
            println("P1: $i/${polynomials.size}")
            p.reductumSet().flatMap { r ->
                val resultants: List<MPoly> = UnivariateResultants.Subresultants(r, r.derivative())
                val leadingCoefficient: MPoly = r.lc()
                resultants + listOf(leadingCoefficient)
            }
        }
        .filter { !it.isConstant }
        .toSet()
        .flatMap { Factor(it).toList() }
        .map { p -> p.divideByLC(p) }
        .toSet()

    //println("Projection1: ${projection1.size}")
    //projection1.forEach { println(it) }

    // compute second projection set as combination of polynomials
    val projection2: MutableSet<MPoly> = HashSet(polynomials.size * polynomials.size)

    val polynomialList = polynomials.toList()
    for (i in 0 until polynomialList.size) {
        println("P2: $i/${polynomials.size}")
        for (j in (i+1) until polynomialList.size) {
            projection2.addAll(UnivariateResultants.Subresultants(polynomialList[i], polynomialList[j]))
        }
    }

    val normalized = projection2
        .filter { !it.isConstant }
        .flatMap { Factor(it).toList() }
        .map { p -> p.divideByLC(p) }

    projection2.clear()
    projection2.addAll(normalized)

    //println("Projection2: ${projection2.size}")
    //projection2.forEach { println(it) }

    return projection1 + projection2
}

fun MPoly.boundary(system: Array<MPoly>): MPoly {
    var result = MultivariatePolynomial.zero(nVariables, ring, ordering)
    for (variable in system.indices) {
        result = result.add(this.derivative(variable).multiply(system[variable]))
    }
    return result
}

/*
/**
 * Extract leading term with respect to given [variable].
 * Unchanged dimensionality.
 */
fun MPoly.leadingTerm(variable: Int): MPoly {
    val degree = this.degree(variable)
    return MultivariatePolynomial.create(nVariables, ring, ordering, this.filter {
        it.exponents[variable] == degree
    })
}

/**
 * Reductum - a polynomial without its leading term (With respect to [variable]).
 * Unchanged dimensionality.
 */
fun MPoly.reductum(variable: Int): MPoly {
    val degree = this.degree(variable)
    return MultivariatePolynomial.create(nVariables, ring, ordering, this.filter {
        it.exponents[variable] != degree
    })
}

fun RED(polynomial: MPoly, variable: Int): List<MPoly> {
    val result = ArrayList<MPoly>()
    var poly = polynomial
    while (!poly.isZero) {
        result.add(poly)
        poly = poly.reductum(variable)
    }
    return result
}

/**
 * Extract leading coefficient with respect to given [variable].
 * The result is an n-1 dimensional polynomial.
 */
fun MPoly.leadingCoefficient(variable: Int): MPoly {
    val degree = this.degree(variable)
    val newTerms = ArrayList<Term>()
    for (monomial in this) {
        if (monomial.exponents[variable] == degree) {
            newTerms.add(monomial.without(variable))
        }
    }
    return MultivariatePolynomial.create(nVariables - 1, ring, ordering, newTerms)
}

/**
 * Compute a partial derivative of a given multivariate polynomial with respect to [variable].
 * Unchanged dimensionality.
 */
fun MPoly.Derivative(variable: Int): MPoly {
    val newTerms = ArrayList<Term>()
    for (monomial in this) {
        val exponent = monomial.exponents[variable]
        if (exponent > 0) { // constant terms are stripped away
            val newDegreeVector = monomial.dvSet(variable, exponent - 1)
            val newCoefficient = monomial.coefficient.multiply(BigInteger.valueOf(exponent))
            newTerms.add(Monomial(newDegreeVector, newCoefficient))
        }
    }
    return MultivariatePolynomial.create(nVariables, ring, ordering, newTerms)
}
*/