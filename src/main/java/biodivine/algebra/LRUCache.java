package biodivine.algebra;

import java.util.LinkedHashMap;
import java.util.Iterator;

public class LRUCache<K, V> {

    private int capacity;
    private LinkedHashMap<K,V> map;

    public LRUCache(int capacity) {
        this.capacity = capacity;
        this.map = new LinkedHashMap<>();
    }

    private long requests = 0;
    private long hits = 0;

    public V get(K key) {
        requests += 1;
        V value = this.map.get(key);
        if (value != null) {
            hits += 1;
            this.set(key, value);
        }
        if (requests % 1000000 == 0) System.out.println("Hit rate "+hits+" / "+requests);
        return value;
    }

    public void set(K key, V value) {
        if (this.map.containsKey(key)) {
            this.map.remove(key);
        } else if (this.map.size() == this.capacity) {
            Iterator<K> it = this.map.keySet().iterator();
            it.next();
            it.remove();
        }
        map.put(key, value);
    }
}
