package biodivine.algebra.synth

import biodivine.algebra.MPoly
import cc.redberry.rings.Rational
import cc.redberry.rings.poly.MultivariateRing

data class Model(
    val ring: MultivariateRing<MPoly>,
    val varNum: Int, val varBounds: Box,
    val paramNum: Int, val paramBounds: Box,
    val equations: List<Rational<MPoly>>
)