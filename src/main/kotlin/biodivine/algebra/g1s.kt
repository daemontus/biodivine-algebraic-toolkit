package biodivine.algebra

import cc.redberry.rings.Rational
import cc.redberry.rings.Ring
import cc.redberry.rings.Rings
import cc.redberry.rings.poly.PolynomialMethods.Factor
import cc.redberry.rings.poly.multivar.GroebnerBases
import cc.redberry.rings.poly.multivar.MonomialOrder
import cc.redberry.rings.poly.multivar.MultivariatePolynomial

fun main(args: Array<String>) {
    val poly = Rings.MultivariateRingQ(2)

    val one = poly.one
    val k1 =  Rational(poly, poly.one, poly.one)
    val k2 = Rational(poly, poly.parse("16"), poly.parse("10"))
    val yE2F1 = Rational(poly, poly.one, poly.parse("10"))
    val a = Rational(poly, poly.one, poly.parse("25"))
    val kp = Rational(poly, poly.parse("5"), poly.parse("100"))
    val yPRB = Rational(poly, poly.parse("11"), poly.parse("1000"))
    val q = Rational(poly, poly.parse("625"), poly.parse("10000"))

    val h0 = Rational(poly, poly.parse("y"), poly.parse("1/2 + y"))
    val h1 = Rational(poly, poly.parse("1/2"), poly.parse("1/2 + x"))
    val f1 = k1.multiply(h0).multiply(h1).subtract(yPRB.multiply(poly.parse("x")))

    val h2 = Rational(poly, poly.parse("16"), poly.parse("16 + y^2"))
    val h3 = Rational(poly, poly.parse("5"), poly.parse("5 + x"))
    val h4 = Rational(poly, poly.parse("y^2"), poly.parse("16 + y^2"))
    val f2 = kp
                .add(k2.multiply(a.multiply(a)).multiply(q).multiply(h2).multiply(h3))
                .add(k2.multiply(h4).multiply(h3))
                .subtract(yE2F1.multiply(poly.parse("y")))

    println("f1: $f1")
    println("f2: $f2")

    val g1 = Factor(f1.numerator()).first()
    val g2 = Factor(f2.numerator()).first()
    println("g1: $g1")
    println("g2: $g2")

    val eq = f1.multiply(f1).add(f2.multiply(f2)).subtract(Rational(poly, poly.one, poly.parse("100000")))

    println("EQ: $eq")
    val ge = Factor(eq.numerator()).first()
    println("gE: ${ge}")

    println("G1 boundary: ${Factor(g1.boundary(poly, arrayOf(f1, f2)))}")
    println("G2 boundary: ${Factor(g2.boundary(poly, arrayOf(f1, f2)))}")
    println("gE boundary: ${Factor(ge.boundary(poly, arrayOf(f1, f2)))}")

    println("G1 boundary groebner: ${GroebnerBases.GroebnerBasis(
        listOf(g1, g1.boundary(poly, arrayOf(f1, f2))), MonomialOrder.EliminationOrder(MonomialOrder.GREVLEX, 1)
    )}")
    println("G2 boundary groebner: ${GroebnerBases.GroebnerBasis(
        listOf(g2, g2.boundary(poly, arrayOf(f1, f2))), MonomialOrder.EliminationOrder(MonomialOrder.GREVLEX, 1)
    )}")
    println("GE boundary groebner: ${GroebnerBases.GroebnerBasis(
        listOf(ge, ge.boundary(poly, arrayOf(f1, f2))), MonomialOrder.DEFAULT
    )}")

    val eRing = Rings.MultivariateRingQ(3)
    val eG1 = g1.insertVariable(2)
    val eG2 = g2.insertVariable(2)
    val bound = eG1.add(eRing.parse("z").multiply(eG2))
    val eSystem = arrayOf(
        Rational(eRing, f1.numerator().insertVariable(2), f1.denominator().insertVariable(2)),
        Rational(eRing, f2.numerator().insertVariable(2), f2.denominator().insertVariable(2))
    )
    val boundBoundary = bound.boundary(eRing, eSystem)
    println("Bound: ${Factor(bound)}")
    println("Bound boundary: ${Factor(boundBoundary)}")
}

fun MPoly.boundary(ring: Ring<MPoly>, system: Array<Rational<MPoly>>): MPoly {
    var result = Rational(ring, ring.one)
    for (variable in system.indices) {
        result = result.add(system[variable].multiply(this.derivative(variable)))
    }
    return result.numerator()
}

fun Rational<MPoly>.derivative(variable: Int): Rational<MPoly> {
    return Rational(ring,
        denominator().copy().multiply(numerator().derivative(variable)).subtract(numerator().copy().multiply(denominator().derivative(variable))),
        denominator().copy().multiply(denominator())
    )
}

