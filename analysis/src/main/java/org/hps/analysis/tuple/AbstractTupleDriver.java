package org.hps.analysis.tuple;

import java.util.ArrayList;
import java.util.List;

import org.lcsim.event.EventHeader;
import org.lcsim.util.Driver;

public abstract class AbstractTupleDriver extends Driver {
    
    private List<TupleCut> cuts = new ArrayList<TupleCut>();
    private String tupleFile = null;
    private TupleWriter writer = null;
    private TupleRecord tuple;

    public TupleRecord getTuple() {
        return tuple;
    }
    
    public void startOfData() {
        tuple = new TupleRecord();
        writer = new TupleWriter(tupleFile);
        addVariables(tuple);
        writer.writeHeader(tuple);
    }
    
    public void setTupleName(String name) {
        this.tuple.setName(name);
    }
    
    public void process(EventHeader event) {
        event.put(tuple.getName(), tuple);
        if (passes(event)) {
            fill(tuple);
        }
    }
    
    public void endOfData() {
        writer.close();
    }
        
    public final void addCuts(TupleCut cut) {
        cuts.add(cut);
    }
    
    public final List<TupleCut> getCuts() {
        return cuts;
    }
    
    public boolean passes(EventHeader event) {
        for (TupleCut cut : cuts) {
            if (!cut.passes(event)) {
                return false;
            }
        }
        return true;
    }
    
    public void setTupleFile(String tupleFile) {
        this.tupleFile = tupleFile;
    }
    
    public abstract void addVariables(TupleRecord record);
    
    public abstract void fill(TupleRecord record);
}
