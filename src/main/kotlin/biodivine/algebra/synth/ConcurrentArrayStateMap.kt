package biodivine.algebra.synth

import biodivine.algebra.params.SemiAlgSet
import biodivine.algebra.params.SemiAlgSolver
import java.util.*
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicReferenceArray

typealias State = Int

class ConcurrentArrayStateMap(
    val capacity: Int, private val solver: SemiAlgSolver
) {

    private val sizeAtomic = AtomicInteger(0)

    val size: Int
        get() = sizeAtomic.get()

    private val data = AtomicReferenceArray<SemiAlgSet?>(capacity)

    fun getOrNull(state: State): SemiAlgSet? = data[state]

    fun get(state: State): SemiAlgSet = data[state] ?: solver.zero

    fun union(state: State, value: SemiAlgSet): Boolean {
        solver.run {
            if (value.isEmpty()) return false
            var current: SemiAlgSet?
            do {
                current = data[state]
                val c = current ?: zero
                if (current != null && value subset current) return false
                val union = c or value
            } while (!data.compareAndSet(state, current, union))
            if (current == null) sizeAtomic.incrementAndGet()
            return true
        }
    }

}