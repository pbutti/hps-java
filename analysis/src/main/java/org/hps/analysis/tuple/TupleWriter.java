package org.hps.analysis.tuple;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

import org.apache.commons.lang3.StringUtils;

public class TupleWriter {
        
    private PrintWriter printWriter;
    
    TupleWriter(String filename) {        
        try {
            this.printWriter = new PrintWriter(filename);
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
    }
    
    void writeHeader(TupleRecord record) {        
        printWriter.println(StringUtils.join(record.getVariables(), ":"));
    }
    
    void writeRecord(TupleRecord record) {
        for (String variable : record.getVariables()) {
            Double value = record.get(variable);
            if (value == null || Double.isNaN(value)) {
                value = -9999.0;
            }
            if (variable.endsWith("/I") || variable.endsWith("/B")) {
                printWriter.format("%d\t", Math.round(value));
            } else {
                printWriter.format("%g\t", value);
            }
        }
        printWriter.println();
    }
    
    void close() {
        printWriter.close();
    }
}
