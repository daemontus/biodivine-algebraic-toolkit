package biodivine.algebra.svg

import biodivine.algebra.NumQ
import biodivine.algebra.svg.Point
import cc.redberry.rings.Rings

enum class Dimension { X, Y }

fun xy(x: NumQ, y: NumQ) = Point(x, y)

val zero = Rings.Q.zero
val p3_10 = Rings.Q.parse("3/10")
val p7_10 = Rings.Q.parse("7/10")
val p1_2 = Rings.Q.parse("1/2")
val n1_2 = Rings.Q.parse("-1/2")
val n1 = Rings.Q.negativeOne
val p1 = Rings.Q.one