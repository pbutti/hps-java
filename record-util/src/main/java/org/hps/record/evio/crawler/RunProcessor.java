package org.hps.record.evio.crawler;

import java.io.File;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.hps.record.epics.EpicsRunProcessor;
import org.hps.record.evio.EvioFileMetaData;
import org.hps.record.evio.EvioFileMetaDataReader;
import org.hps.record.evio.EvioFileSource;
import org.hps.record.evio.EvioLoop;
import org.hps.record.run.RunSummary;
import org.hps.record.scalers.ScalersEvioProcessor;
import org.lcsim.util.log.DefaultLogFormatter;
import org.lcsim.util.log.LogUtil;

/**
 * Processes EVIO files from a run in order to extract various meta data information including start and end dates.
 * <p>
 * This class is a wrapper for activating different sub-tasks, including optionally caching all files from the JLAB MSS
 * to the cache disk.
 * <p>
 * There is also a list of processors which is run on all events from the run.
 *
 * @author Jeremy McCormick, SLAC
 */
final class RunProcessor {

    /**
     * Setup logger.
     */
    private static final Logger LOGGER = LogUtil.create(RunProcessor.class, new DefaultLogFormatter(), Level.FINE);

    /**
     * Process all the runs that were found.
     *
     * @param runs the run log containing the list of run summaries
     * @throws Exception if there is an error processing one of the runs
     */
    static void processAllRuns(final JCacheManager cacheManager, final RunSummaryMap runs, final CrawlerConfig config)
            throws Exception {

        // Process all of the runs that were found.
        for (final RunSummary runSummary : runs.getRunSummaries()) {

            // Clear the cache manager.
            if (config.useFileCache()) {
                LOGGER.info("clearing file cache");
                cacheManager.clear();
            }

            // Create a processor to process all the EVIO events in the run.
            LOGGER.info("creating run processor for " + runSummary.getRun());
            final RunProcessor runProcessor = new RunProcessor(cacheManager, runSummary, config);

            // Process all of the run's files.
            LOGGER.info("processing run " + runSummary.getRun());
            runProcessor.processRun();
        }
    }

    /**
     * The cache manager.
     */
    private final JCacheManager cacheManager;

    /**
     * Processor for extracting EPICS information.
     */
    private final EpicsRunProcessor epicsLog;

    /**
     * The data source with the list of EVIO files to process.
     */
    private final EvioFileSource evioFileSource;

    /**
     * The EVIO event processing loop.
     */
    private final EvioLoop evioLoop = new EvioLoop();

    /**
     * The run summary information updated by running this processor.
     */
    private final RunSummary runSummary;

    /**
     * Processor for extracting scaler data.
     */
    private final ScalersEvioProcessor scalersProcessor;

    /**
     * Set to <code>true</code> to use file caching.
     */
    private boolean useFileCache;

    /**
     * Create a run processor.
     *
     * @param runSummary the run summary object for the run
     * @return the run processor
     */
    RunProcessor(final JCacheManager cacheManager, final RunSummary runSummary, final CrawlerConfig config) {

        this.runSummary = runSummary;
        this.cacheManager = cacheManager;

        // Setup record loop.
        runSummary.sortFiles();
        evioFileSource = new EvioFileSource(runSummary.getEvioFileList());
        evioLoop.setEvioFileSource(evioFileSource);

        // Add EPICS processor.
        epicsLog = new EpicsRunProcessor();
        evioLoop.addEvioEventProcessor(epicsLog);

        // Add Scaler data processor.
        scalersProcessor = new ScalersEvioProcessor();
        scalersProcessor.setResetEveryEvent(false);
        evioLoop.addEvioEventProcessor(scalersProcessor);

        // Set whether file caching from MSS is enabled.
        this.useFileCache(config.useFileCache());
    }

    /**
     * Cache all files and wait for the operation to complete.
     * <p>
     * Potentially, this operation can take a very long time. This can be managed using the
     * {@link JCacheManager#setWaitTime(long)} method to set a timeout.
     */
    private void cacheFiles() {

        LOGGER.info("caching files from run " + this.runSummary.getRun());

        // Cache all the files and wait for the operation to complete (it will take awhile!).
        this.cacheManager.cache(this.runSummary.getEvioFileList());
        final boolean cached = this.cacheManager.waitForCache();

        // If the files weren't cached then die.
        if (!cached) {
            throw new RuntimeException("The cache process did not complete in time.");
        }

        LOGGER.info("done caching files from run " + this.runSummary.getRun());
    }

    /**
     * Process the run by executing the registered {@link org.hps.record.evio.EvioEventProcessor}s and extracting the
     * start and end dates.
     * <p>
     * This method will also execute file caching from MSS, if enabled by the {@link #useFileCache} option.
     *
     * @throws Exception if there is an error processing a file
     */
    void processRun() throws Exception {

        LOGGER.info("processing " + this.runSummary.getEvioFileList().size() + " files from run "
                + this.runSummary.getRun());

        // Cache files from MSS if this is enabled.
        if (this.useFileCache) {
            LOGGER.info("caching files from MSS");
            this.cacheFiles();
        }

        // Run processors over all files.
        LOGGER.info("looping over all events");
        evioLoop.loop(-1);

        // Get run start date.
        LOGGER.info("setting run start date");
        this.setRunStartDate();

        // Get run end date.
        LOGGER.info("setting run end date");
        this.setRunEndDate();

        // Update run summary from processors.
        LOGGER.info("updating run summary");
        this.updateRunSummary();

        LOGGER.info("done processing run " + this.runSummary.getRun());
    }

    /**
     * Set the run end date by getting meta data from the last file and copying it to the run summary.
     */
    private void setRunEndDate() {
        final File lastEvioFile = runSummary.getEvioFileList().get(runSummary.getEvioFileList().size() - 1);
        LOGGER.info("getting meta data for " + lastEvioFile.getPath());
        final EvioFileMetaDataReader metaDataReader = new EvioFileMetaDataReader();
        final EvioFileMetaData metaData = metaDataReader.getMetaData(lastEvioFile);
        LOGGER.info(metaData.toString());
        LOGGER.info("setting unix end time to " + metaData.getEndDate().getTime() + " from meta data");
        runSummary.setEndTimeUtc(metaData.getEndDate().getTime());
        runSummary.setEndOkay(metaData.hasEnd());
    }

    /**
     * Set the run start date by getting meta data from the first file and copying it to the run summary.
     */
    private void setRunStartDate() {
        final File firstEvioFile = runSummary.getEvioFileList().get(0);
        LOGGER.info("getting meta data for " + firstEvioFile.getPath());
        final EvioFileMetaDataReader metaDataReader = new EvioFileMetaDataReader();
        final EvioFileMetaData metaData = metaDataReader.getMetaData(firstEvioFile);
        LOGGER.info(metaData.toString());
        LOGGER.info("setting unix start time to " + metaData.getStartDate().getTime() + " from meta data");
        runSummary.setStartTimeUtc(metaData.getStartDate().getTime());
    }

    /**
     * Update the current run summary by copying data to it from the EVIO processors.
     */
    private void updateRunSummary() {
        // Put scaler data from EVIO processor into run summary.
        runSummary.setScalerData(this.scalersProcessor.getScalerData());

        // Set total number of events on the run summary from the event counter.
        runSummary.setTotalEvents((int) evioLoop.getTotalCountableConsumed());

        // Set EpicsData for the run.
        runSummary.setEpicsData(this.epicsLog.getEpicsData());
    }

    /**
     * Set whether or not to use the file caching, which copies files from the JLAB MSS to the cache disk.
     * <p>
     * Since EVIO data files at JLAB are primarily kept on the MSS, running without this option enabled there will
     * likely cause the job to fail.
     *
     * @param cacheFiles <code>true</code> to enabled file caching
     */
    void useFileCache(final boolean cacheFiles) {
        this.useFileCache = cacheFiles;
        LOGGER.config("file caching enabled");
    }
}
