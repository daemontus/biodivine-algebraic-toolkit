package biodivine.algebra.synth

import biodivine.algebra.params.SemiAlgSet
import biodivine.algebra.params.SemiAlgSolver

typealias StateMap = ConcurrentArrayStateMap

class SemiAlgTransitionSystem(
    val solver: SemiAlgSolver,
    val states: List<Box>, transitions: List<Triple<Int, Int, SemiAlgSet>>
) {

    private val stateCount = states.size

    private fun newMap() = ConcurrentArrayStateMap(stateCount, solver)

    private val successors: Map<Int, List<Pair<Int, SemiAlgSet>>> = transitions
        .groupBy({ it.first }, { it.second to it.third }).withDefault { emptyList() }

    private val predecessor: Map<Int, List<Pair<Int, SemiAlgSet>>> = transitions
        .groupBy({ it.second }, { it.first to it.third }).withDefault { emptyList() }

    fun StateMap.reachBackward(guard: StateMap? = null): StateMap {
        val shouldUpdate = RepeatingConcurrentStateQueue(stateCount)
        val result = newMap()
        // init reach
        for (s in 0 until stateCount) {
            val c = this.getOrNull(s)
            if (c != null) {
                result.union(s, this.get(s))
                shouldUpdate.set(s)
            }
        }
        println("Start reach backward.")
        // repeat
        var i = 0
        fork {
            var state = shouldUpdate.next(0)
            while (state > -1) {
                while (state > -1) {
                    solver.run {
                        i += 1
                        if (i % 50 == 0) println("Processed $i states")
                        for ((source, edgeParams) in predecessor.getValue(state)) {
                            val bound = if (guard == null) result.get(state) else {
                                result.get(state) and guard.get(source)
                            }
                            val original = result.get(source)
                            // update target -> if changed, mark it as working
                            val changed = result.union(source, edgeParams and bound)
                            if (changed) {
                                //println("Update state $source${states[source]} from $state${states[state]} with value ${edgeParams and bound} where previously the value was: ${original}")
                                //println("After update, the value is ${result.get(source)}")
                                shouldUpdate.set(source)
                            }
                        }
                    }
                    state = shouldUpdate.next(state + 1)
                }
                // double check - maybe someone added another thing
                state = shouldUpdate.next(0)
            }
        }

        return result
    }

}