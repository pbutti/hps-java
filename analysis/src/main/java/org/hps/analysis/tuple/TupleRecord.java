package org.hps.analysis.tuple;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.lcsim.event.EventHeader;

/**
 * Simple class for representing a single N-tuple record as a map.
 */
public class TupleRecord {
    
    private Map<String, Double> dataMap = new HashMap<String, Double>();
    private String name = "TupleData";
    
    public TupleRecord() {        
    }
    
    public TupleRecord(String name) {
        this.name = name;
    }
    
    public String getName() {
        return this.name;
    }
    
    public void setName(String name) {
        this.name = name;
    }
    
    public final void fill(String field, double value) {
        dataMap.put(field, value);
    }
    
    public double get(String field) {
        return dataMap.get(field);
    }
    
    public void add(String field) {
        dataMap.put(field, null);
    }
    
    public void addAll(String[] fields) {
        for (String field : fields) {                   
            dataMap.put(field, null);
        }
    }
    
    public boolean hasVariable(String field) {
        return dataMap.containsKey(field);
    }
    
    public boolean filled(String field) {
        return dataMap.containsKey(field) && dataMap.get(field) != null;
    }
    
    public Set<String> getVariables() {
        return dataMap.keySet();
    }
    
    public void clear() {
        for (String key : this.dataMap.keySet()) {
            this.dataMap.put(key, null);
        }
    }
    
    public static void write(EventHeader event, TupleRecord rec) {
        event.put(rec.getName(), rec);
    }
    
    public static TupleRecord get(EventHeader event, String name) {
        return (TupleRecord) event.get(name);
    }       
        
}
