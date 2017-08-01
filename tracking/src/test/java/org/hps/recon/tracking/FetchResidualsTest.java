package org.hps.recon.tracking;

/**
 * Test class for raw->reco LCIO + fetching residuals from reco.
 * @author mdiamond <mdiamond@slac.stanford.edu>
 */
public class FetchResidualsTest extends ReconTestSkeleton {

    static final String inputFileName = "ap_prompt_raw.slcio";

    @Override
    public void testRecon() throws Exception {

        testInputFileName = inputFileName;
        nEvents = 10;
        testTrackingDriver = new FetchResiduals();
        super.testRecon();

    }

}
