package biodivine.algebra.synth

import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors

val parallelism = Runtime.getRuntime().availableProcessors()
val pool: ExecutorService = Executors.newWorkStealingPool(parallelism)//.newFixedThreadPool(parallelism)

fun fork(action: () -> Unit) {
    (1..parallelism).map {
        pool.submit(action)
    }.map { it.get() }
}

fun <T, R> List<T>.mapParallel(action: (T) -> R): List<R> {
    return this.map {
        pool.submit<R> { action(it) }
    }.map { it.get() }
}