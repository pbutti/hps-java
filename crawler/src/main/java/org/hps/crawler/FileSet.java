package org.hps.crawler;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.hps.datacat.client.DatasetFileFormat;

/**
 * Map of file format to a list of files.
 * 
 * @author Jeremy McCormick, SLAC
 */
final class FileSet extends HashMap<DatasetFileFormat, List<File>> {
    
    List<File> get(DatasetFileFormat format) {
        if (super.get(format) == null) {
            this.put(format, new ArrayList<File>());
        }
        return super.get(format);
    }
    
    void addFile(DatasetFileFormat format, File file) {
        this.get(format).add(file);
    }
}
