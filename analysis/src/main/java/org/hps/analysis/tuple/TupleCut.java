package org.hps.analysis.tuple;

import java.util.HashMap;
import java.util.Map;

import org.lcsim.event.EventHeader;

public abstract class TupleCut {
    
    private Map<String, Object> parameters = new HashMap<String, Object>();
    
    public void setParameter(String name, Object value) {
        parameters.put(name, value);
    }
    
    public <T> T getParameter(String name) {
        return (T) parameters.get(name);
    }
    
    public abstract boolean passes(EventHeader event);
}
