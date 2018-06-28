package org.hps.recon.tracking;

import java.io.File;
//import java.net.URL;

import junit.framework.TestCase;

import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.detector.svt.SvtDetectorSetup;
//import org.lcsim.util.cache.FileCache;
//import org.hps.job.DatabaseConditionsManagerSetup;
//import org.lcsim.util.cache.FileCache;
//import org.lcsim.util.loop.LCIODriver;
import org.lcsim.util.loop.LCSimLoop;
import org.lcsim.recon.tracking.digitization.sisim.config.RawTrackerHitSensorSetup;
//import org.lcsim.recon.tracking.digitization.sisim.config.RawTrackerHitSensorSetup;
//import org.lcsim.job.ConditionsSetup;
import org.lcsim.recon.tracking.digitization.sisim.config.ReadoutCleanupDriver;
import java.util.ArrayList;
import java.util.List;
//import org.lcsim.util.test.TestUtil.TestOutputFile;

/**
 * Test class to create set of histograms (aida/root) from reco LCIO.
 * 
 * @author Miriam Diamond <mdiamond@slac.stanford.edu> 
 */
public class TrackingReconstructionPlotsTest extends TestCase {

    static final String testInput1 = "/nfs/slac/g/hps3/data/mc_production/beamRotationStudy/physrun2016/x31_yneg0pt5_withCorrectedZTilt/recon/WBT/2pt3/v5-3_globalAlign/4.1_2017Dec28/pairs1/wabv3AF-egsv6-triv2AF_10to1_HPS-PhysicsRun2016-v5-3-fieldmap_globalAlign_4.1_2017Dec28_pairs1_";
    static final String testInput2 = ".slcio";
    static final String testURLBase = "http://www.lcsim.org/test/hps-java";

    private final int nEvents = -1;

    public void testTrackRecoPlots() throws Exception {
        //URL testURL = new URL(testURLBase + "/" + testInput);
        //FileCache cache = new FileCache();
        //File lcioInputFile = cache.getCachedFile(testURL);
        int fileNum=0;
        String aidaOutput = "/nfs/slac/work/mdiamond/WABtilt.root";

        final DatabaseConditionsManager manager = new DatabaseConditionsManager();
        manager.addConditionsListener(new SvtDetectorSetup());

        LCSimLoop loop2 = new LCSimLoop();

        RawTrackerHitSensorSetup rthss = new RawTrackerHitSensorSetup();
        String[] readoutColl = { "SVTRawTrackerHits" };
        rthss.setReadoutCollections(readoutColl);
        loop2.add(rthss);

        TrackingReconstructionPlots trp = new TrackingReconstructionPlots();
        trp.setOutputPlots(aidaOutput);
        loop2.add(trp);

        ReadoutCleanupDriver rcd = new ReadoutCleanupDriver();
        loop2.add(rcd);

        while (fileNum < 1001) {
            String testInput = String.format("%s%d%s", testInput1, fileNum, testInput2);
            //String aidaOutput = String.format("/nfs/slac/work/mdiamond/WABtilt%d.root", fileNum);
            File lcioInputFile = new File(testInput);
            fileNum++;
            if (!lcioInputFile.exists()) {
                continue;
            }

            loop2.setLCIORecordSource(lcioInputFile);
            loop2.loop(nEvents, null);
        }
        loop2.dispose();
    }

}
