package biodivine.algebra.params

import biodivine.algebra.MPoly
import biodivine.algebra.NumQ
import biodivine.algebra.ia.*
import biodivine.algebra.rootisolation.DescartWithSquareFreeFactorisationRootIsolation
import biodivine.algebra.synth.Box
import cc.redberry.rings.Rings
import cc.redberry.rings.Rings.Q
import cc.redberry.rings.poly.IPolynomialRing
import cc.redberry.rings.poly.PolynomialMethods.Factor
import cc.redberry.rings.poly.multivar.MonomialOrder
import cc.redberry.rings.poly.multivar.MultivariatePolynomial
import cc.redberry.rings.poly.univar.UnivariateResultants

fun List<MPoly>.normalize(): List<MPoly> = this.filter { !it.isConstant }.flatMap { Factor(it) }.map { it.divideByLC(it) }
fun MPoly.canHaveZero(bounds: Box): Boolean = this.evaluate(*bounds.data).hasZero