package org.hps.analysis.geometry;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Set;
import junit.framework.TestCase;
import org.hps.data.detectors.DetectorDataResources;
import org.lcsim.util.cache.FileCache;
import org.lcsim.util.loop.LCSimLoop;
import org.lcsim.util.test.TestUtil;

/**
 *
 * @author ngraf
 */
public class DrawDetectorTest extends TestCase {

    int nEvents = 1;

    public void testIt() throws Exception {
        Set<String> availableDetectors = DetectorDataResources.getDetectorNames();
        for (String s : availableDetectors) {
            System.out.println(s);
        }
        new TestUtil.TestOutputFile(this.getClass().getSimpleName()).mkdir();

        FileCache cache = new FileCache();
        String fileName = "ap2.2gev150mev_slic-3.1.5_geant4-v9r6p1_QGSP_BERT_HPS-Proposal2014-v8-2pt2-42-1.slcio";
//        fileName = "ap2.2gev150mev_slic-3.1.5_geant4-v9r6p1_QGSP_BERT_HPS-Proposal2014-v8-2pt2.slcio";
//        File inputFile = cache.getCachedFile(new URL("http://www.lcsim.org/test/hps-java/APrimeReconTest/" + fileName));
//        File inputFile = cache.getCachedFile(new URL("http://www.lcsim.org/test/hps-java/HpsGblRefitterTest.slcio"));
        File inputFile = cache.getCachedFile(new URL("http://www.lcsim.org/test/hps-java/ap100mev_HPS-ECalCommissioning-v2_twotrackEvents-0-100.slcio")); //HpsGblRefitterTest.slcio"));
        // Read in the LCIO event file and print out summary information.
        System.out.println("Running DrawTracker ...");
        LCSimLoop loop = new LCSimLoop();
        loop.add(new DrawDetector());
        try {
            loop.setLCIORecordSource(inputFile);
            loop.loop(nEvents);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        System.out.println("Loop processed " + loop.getTotalSupplied() + " events.");
        System.out.println("Done!");
    }
}
