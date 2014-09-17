package org.hps;

import hep.aida.IHistogram1D;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import junit.framework.TestCase;

import org.hps.users.jeremym.MockDataChallengeDiagnosticDriver;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Track;
import org.lcsim.job.AidaSaveDriver;
import org.lcsim.job.JobControlManager;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.lcsim.util.cache.FileCache;
import org.lcsim.util.loop.LCSimLoop;

/**
 * This test runs the standard reconstruction on a small set of Mock Data events generated by the
 * trigger simulation. It then generates a set of diagnostic plots which are checked against an
 * answer key.
 * 
 * @author Jeremy McCormick <jeremym@slac.stanford.edu>
 */
public class MockDataReconTest extends TestCase {

    static final String className = MockDataReconTest.class.getSimpleName();
    static final File outputDir = new File("./target/test-output/" + className);
    static final File outputFile = new File(outputDir.getAbsolutePath() + File.separator + className);
    static final File reconFile = new File(outputFile.getAbsolutePath() + ".slcio");
    static final File aidaFile = new File(outputFile.getAbsolutePath() + ".aida");    

    static final String steeringResource = "/org/hps/steering/test/MockDataReconTest.lcsim";
    //static final String steeringResource = "/org/hps/steering/recon/HPS2014OfflineTruthRecon.lcsim";
    
    // TODO: Get some values for these and add test assertions!
    /*
    static final int expectedReconstructedParticles = 0;
    static final int expectedClusters = 0;
    static final int expectedTracks = 0;
    static final int expectedCalorimeterHits = 0;
    static final int expectedMCParticles = 0;
    */

    static final String clusterCollectionName = "EcalClusters";
    static final String reconstructedParticleCollectionName = "FinalStateParticles";
    static final String trackCollectionName = "MatchedTracks";

    AIDA aida = AIDA.defaultInstance();

    static final String mockDataUrl = "http://www.slac.stanford.edu/~meeg/hps2/meeg/mock_data/tritrig-beam-tri_1-10_readout.slcio";

    public void setUp() {
        // Delete files if they already exist.
        if (reconFile.exists())
            reconFile.delete();
        if (aidaFile.exists())
            aidaFile.delete();
        
        // Create output dir.
        outputDir.mkdirs();
        if (!outputDir.exists()) {
            throw new RuntimeException("Failed to create test output dir.");
        }
    }
    
    public void testReconMockData() {

        // Run the reconstruction over input events.
        runRecon();

        // Create the plots.
        createPlots();

        // Check the statistics of the plots.
        //checkPlots();
    }

    private void runRecon() {
        
        System.out.println("caching file ...");
        System.out.println(mockDataUrl);
        
        File mockDataFile = null;
        try {
            FileCache cache = new FileCache();
            mockDataFile = cache.getCachedFile(new URL(mockDataUrl));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        System.out.println("running recon using steering resource " + steeringResource);
        JobControlManager jobManager = new JobControlManager();
        jobManager.addVariableDefinition("outputFile", outputFile.getPath());
        jobManager.addInputFile(mockDataFile);
        jobManager.setup(steeringResource);
        jobManager.run();
    }

    private void createPlots() {
        LCSimLoop loop = new LCSimLoop();
        loop.add(new MockDataChallengeDiagnosticDriver());
        loop.add(new CheckDriver());
        AidaSaveDriver aidaSaveDriver = new AidaSaveDriver();
        aidaSaveDriver.setOutputFileName(aidaFile.getAbsolutePath());
        loop.add(aidaSaveDriver);        
        try {
            loop.setLCIORecordSource(reconFile);
            loop.loop(-1);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /*
    private void checkPlots() {

        IHistogram1D reconParticleCountPlot = aida.histogram1D("/" + reconstructedParticleCollectionName + "/Item Count");
        System.out.println("reconParticleCount sumBinHeights = " + reconParticleCountPlot.sumBinHeights());
        System.out.println("reconParticleCount mean = " + reconParticleCountPlot.mean());
        System.out.println("reconParticleCount itemCount = " + computeItemCount(reconParticleCountPlot));

        IHistogram1D trackCountPlot = aida.histogram1D("/" + trackCollectionName + "/Item Count");
        System.out.println("trackCount Plot sumBinHeights = " + trackCountPlot.sumBinHeights());
        System.out.println("trackCount itemCount = " + computeItemCount(trackCountPlot));
        //assertEquals("Wrong number of tracks counted.", expectedTracks, computeItemCount(trackCountPlot));         
    }
     */
    
    /**
     * Compute an item count based on a histogram with bins of size 1.
     * @param histogram
     * @return
     */
    /*
    private int computeItemCount(IHistogram1D histogram) {
        if (histogram.axis().binUpperEdge(0) - histogram.axis().binLowerEdge(0) != 1.0)
            throw new IllegalArgumentException("The bins are the wrong size for this method.");
        if (!histogram.title().equals("Item Count"))
            throw new IllegalArgumentException("The histogram doesn't appear to contain item counts.");
        int nbins = histogram.axis().bins();
        double total = 0;
        for (int i = 0; i < nbins; i++) {
            total += histogram.binHeight(i) * histogram.axis().binLowerEdge(i);
        }
        if (total % 1 != 0)
            throw new RuntimeException("Nonsensical number of items computed: " + total);
        return (int)total;
    }
    */
    
    static class CheckDriver extends Driver {
        
        int ntracks;
        int nparticles;
        int nclusters;
        int nevents;
        
        public void process(EventHeader event) {
            ++nevents;
            ntracks += event.get(Track.class, trackCollectionName).size(); 
            nparticles += event.get(ReconstructedParticle.class, reconstructedParticleCollectionName).size();
            nclusters += event.get(Cluster.class, clusterCollectionName).size();
        }
        
        public void endOfData() {
            System.out.println("CheckDriver got the following ...");
            System.out.println("  nevents = " + nevents);
            System.out.println("  ntracks = " + ntracks);
            System.out.println("  nparticles = " + nparticles);
            System.out.println("  nclusters = " + nclusters);
            System.out.println("  <ntracks / nevents> = " + ((double)ntracks / (double)nevents));
            System.out.println("  <nparticles / nevents> = " + ((double)nparticles / (double)nevents));
            System.out.println("  <nclusters / nevents> = " + ((double)nclusters / (double)nevents));
        }
        
    }
}
