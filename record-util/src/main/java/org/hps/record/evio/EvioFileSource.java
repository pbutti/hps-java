package org.hps.record.evio;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.freehep.record.source.AbstractRecordSource;
import org.freehep.record.source.NoSuchRecordException;
import org.jlab.coda.jevio.EvioEvent;
import org.jlab.coda.jevio.EvioException;
import org.jlab.coda.jevio.EvioReader;

/**
 * A basic implementation of an <tt>AbstractRecordSource</tt> for supplying <tt>EvioEvent</tt> objects to a loop from
 * EVIO files.
 * <p>
 * Unlike the LCIO record source, it has no rewind or indexing capabilities.
 *
 * @author <a href="mailto:jeremym@slac.stanford.edu">Jeremy McCormick</a>
 */
public final class EvioFileSource extends AbstractRecordSource {

    /**
     * The current event.
     */
    private EvioEvent currentEvent;

    /**
     * The current file index.
     */
    private int fileIndex = 0;

    /**
     * The list of input data files.
     */
    private final List<File> files = new ArrayList<File>();

    /**
     * The reader to use for reading and parsing the EVIO data.
     */
    private EvioReader reader;

    /**
     * Constructor taking a single EVIO file.
     *
     * @param file the EVIO file
     */
    public EvioFileSource(final File file) {
        this.files.add(file);
        openReader();
    }

    /**
     * Constructor taking a list of EVIO files.
     *
     * @param files the list of EVIO files
     */
    public EvioFileSource(final List<File> files) {
        this.files.addAll(files);
        openReader();
    }

    /**
     * Close the current reader.
     */
    private void closeReader() {
        try {
            this.reader.close();
        } catch (final IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Return <code>true</code> if there are no more files to open from the list.
     *
     * @return <code>true</code> if there are no more files in the list
     */
    boolean endOfFiles() {
        return this.fileIndex > this.files.size() - 1;
    }

    /**
     * Get the current record which is an <code>EvioEvent</code>.
     *
     * @return the current record
     */
    @Override
    public Object getCurrentRecord() throws IOException {
        return this.currentEvent;
    }

    /**
     * Return <code>true</code> if there is a current record loaded.
     *
     * @return <code>true</code> if there is a current record loaded.
     */
    @Override
    public boolean hasCurrent() {
        return this.currentEvent != null;
    }

    /**
     * Return <code>true</code> if there is a next record.
     *
     * @return <code>true</code> if there are more records to load
     */
    @Override
    public boolean hasNext() {
        try {
            return this.reader.getNumEventsRemaining() != 0;
        } catch (IOException | EvioException e) {
            throw new RuntimeException("Error getting num remaining events.");
        }
    }

    /**
     * Load the next record.
     *
     * @throws NoSuchRecordException if source is exhausted
     * @throws IOException if there is an error creating the next <code>EvioEvent</code>
     */
    @Override
    public void next() throws IOException, NoSuchRecordException {
        for (;;) {
            try {
                this.currentEvent = this.reader.parseNextEvent();
            } catch (final EvioException e) {
                throw new IOException(e);
            }
            if (this.currentEvent == null) {
                closeReader();
                this.fileIndex++;
                if (!endOfFiles()) {
                    openReader();
                    continue;
                } else {
                    throw new NoSuchRecordException();
                }
            }
            return;
        }
    }

    /**
     * Open the next file in the list with the reader.
     *
     * @throws RuntimeException if an <code>EvioException</code> or <code>IOException</code> occurs while opening file
     */
    private void openReader() {
        try {
            System.out.println("Opening reader for file " + this.files.get(this.fileIndex) + " ...");
            this.reader = new EvioReader(this.files.get(this.fileIndex), false,true);
            System.out.println("Done opening file.");
        } catch (EvioException | IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Returns <code>true</code> to indicate next record capability is supported.
     *
     * @return <code>true</code> to indicate next record capability is supported
     */
    @Override
    public boolean supportsNext() {
        return true;
    }
}