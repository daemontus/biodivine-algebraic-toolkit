package biodivine.algebra.params

import biodivine.algebra.MPoly
import cc.redberry.rings.poly.univar.UnivariateResultants

object Projection {

    fun discriminant(poly: MPoly, variable: Int): List<MPoly> {
        val u = poly.asUnivariate(variable)
        return UnivariateResultants.Subresultants(u, u.derivative()).normalize()
    }

    fun resultant(a: MPoly, b: MPoly, variable: Int): List<MPoly> {
        return UnivariateResultants.Subresultants(
            a.asUnivariate(variable), b.asUnivariate(variable)
        ).normalize()
    }

}