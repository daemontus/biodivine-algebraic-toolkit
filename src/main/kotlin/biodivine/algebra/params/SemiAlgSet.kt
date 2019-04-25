package biodivine.algebra.params

import biodivine.algebra.MPoly
import cc.redberry.rings.Ring

/**
 * Represents a semi-algebraic set defined by the given set of algebraic varieties (0 = polynomial).
 * The valid cells of the decomposition are then given in the cells array.
 *
 * Note that the global parameter bounds are not stored in the list of polynomials and the solver
 * should apply them automatically.
 */
data class SemiAlgSet(
    val ring: Ring<MPoly>,
    val bounds: Set<MPoly>,
    val cells: List<IntArray>
)

class SemiAlgSolver(
    private val params: Int,
    private val ring: Ring<MPoly>
) {

    val zero: SemiAlgSet = SemiAlgSet(ring, emptySet(), emptyList())
    // if there are no bounds, there is only one cell: the whole box - all coordinates zero
    val one: SemiAlgSet = SemiAlgSet(ring, emptySet(), listOf(IntArray(params)))

    /*
        General overview of the operations:

        1.  Create a "common" CADPolynomials - a CADPolynomials of the union of the two bound polynomial sets.
        2.  Translate each cell set into this common CADPolynomials domain - i.e. for each cell/sample point in the new cad,
            check if it belongs to a cell in the old cad.
        3.  Perform the operation cell-wise, since now we have a common representation for both sets.
        4.  Remove redundant polynomials. How do you do this? Well, there might be other, more efficient options,
        but current idea is like this: We try to remove a polynomial from the CADPolynomials, we then compute the new CADPolynomials and
        for each cell in the original CADPolynomials, we compute the ID in the new CADPolynomials - if two incompatible cells are merged,
        we know we removed a wrong polynomial. All polynomials are independent in this sense - but remember that
        there might be 3 polynomials which are redundant, but by removing first two, the third is no longer redundant
        (because they have common roots!). So always check redundancy of other polynomials on the already reduced CADPolynomials.
     */



}