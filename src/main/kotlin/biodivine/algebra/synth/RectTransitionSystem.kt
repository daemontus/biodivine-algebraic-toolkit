package biodivine.algebra.synth

data class RectTransitionSystem(
    val states: List<Box>, val transitions: List<Pair<Int, Int>>
) {

    private val successors: Map<Int, List<Int>> = transitions.groupBy({ it.first }, { it.second }).withDefault { emptyList() }
    private val predecessors: Map<Int, List<Int>> = transitions.groupBy({ it.second }, { it.first }).withDefault { emptyList() }

    fun nextForward(set: Set<Int>): Set<Int> = set.flatMapTo(HashSet()) { successors.getValue(it) }
    fun nextBackward(set: Set<Int>): Set<Int> = set.flatMapTo(HashSet()) { predecessors.getValue(it) }

    fun reachForward(set: Set<Int>): Set<Int> {
        val result = HashSet(set)
        var recompute = set
        while (recompute.isNotEmpty()) {
            val next = nextForward(recompute)
            recompute = next - result
            result += next
        }
        return result
    }

    fun reachBackward(set: Set<Int>): Set<Int> {
        val result = HashSet(set)
        var recompute = set
        while (recompute.isNotEmpty()) {
            val next = nextBackward(recompute)
            recompute = next - result
            result += next
        }
        return result
    }

}