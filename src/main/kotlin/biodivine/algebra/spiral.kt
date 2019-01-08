package biodivine.algebra

import cc.redberry.rings.Rings
import cc.redberry.rings.poly.PolynomialMethods.Factor
import cc.redberry.rings.poly.multivar.GroebnerBases
import cc.redberry.rings.poly.multivar.MonomialOrder

fun main(args: Array<String>) {
    val vars = arrayOf("x", "y")    // variable names are inferred automatically, but lets make sure
    val ring = Rings.MultivariateRingQ(2)
    val coder = ring.mkCoder(*vars)

    val p1 = coder.parse("-x + 4*y")
    val p2 = coder.parse("-x - y")
    val pFinal = coder.parse("x^2 + y^2 - 1")
    val pEq = p1.copy().multiply(p1).add(p2.copy().multiply(p2)).subtract(coder.parse("1/2"))
    val system = arrayOf(p1, p2)

    val eRing = Rings.MultivariateRingQ(3)
    val eSystem = system.map { it.insertVariable(2) }.toTypedArray()
    val eBound = pFinal.insertVariable(2).add(eRing.parse("z").multiply(p2.insertVariable(2)))
    val eBoundBoundary = eBound.boundary(eSystem)

    println("Final boundary: ${pFinal.boundary(system)} boundary boundary: ${pFinal.boundary(system).boundary(system)}")

    println("Bound: $eBound and boundary: $eBoundBoundary")

    val eBoundFinal = eRing.parse("(x^2 + y^2 - 1) + 12/100*(-x + 4*y)")
    val eBound2 = eBoundFinal.add(eRing.parse("z").multiply(p2.insertVariable(2)))
    val eBound2Boundary = eBound2.boundary(eSystem)

    println("Bound2: $eBound2 and boundary: $eBound2Boundary")

    val p = eRing.parse("(x^2 + y^2 - 1) - 7/10*(-x-y)-1/10*(-x+4*y)")

/*
    val eqBorder = Factor(pEq.boundary(system))[0]
    val pFinalBorder = pFinal.boundary(system)
    println("pEQ = $pEq")
    println("EQ border: $eqBorder")
    println("EQ border border: ${eqBorder.boundary(system)}")
    println("pFinal = $pFinal")
    println("pFinal border: ${pFinal.boundary(system)}")
    println("pFinal border: ${pFinal.boundary(system).boundary(system)}")
    val polynomials = listOf(p1, p2, pFinal, pEq)

    println("Groebner: ${GroebnerBases.GroebnerBasis(listOf(pFinal, p1), MonomialOrder.DEFAULT)}")

    val pBound = pFinalBorder.copy()
    println(pBound)

    val eRing = Rings.MultivariateRingQ(3)
    val eBound = pFinal.insertVariable(2).add(eRing.parse("z").multiply(p1.insertVariable(2)))
    val eBoundBoundary = eBound.boundary(system.map { it.insertVariable(2) }.toTypedArray())

    println("Ebound: $eBound")
    println("Ebound border: $eBoundBoundary")



    val proj = project(listOf(p1.insertVariable(2), p2.insertVariable(2), eBound, eBoundBoundary).map {
        it.asUnivariate(0)
    }.toSet())
    val proj2 = project(proj.map { it.asUnivariate(1) }.toSet())

    println("Proj: ${proj2.size}")
    proj2.forEach { println(it) }
*/


}