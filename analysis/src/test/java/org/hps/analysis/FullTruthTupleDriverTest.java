package org.hps.analysis;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.hps.analysis.MC.MCFullDetectorTruth;
import org.hps.analysis.MC.TrackTruthMatching;
import org.hps.conditions.beam.BeamEnergy;
import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.detector.svt.SvtDetectorSetup;
import org.hps.recon.tracking.TrackUtils;
import org.lcsim.recon.tracking.digitization.sisim.config.RawTrackerHitSensorSetup;
import org.lcsim.recon.tracking.digitization.sisim.config.ReadoutCleanupDriver;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.lcsim.util.loop.LCSimLoop;
import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.event.MCParticle;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.geometry.Detector;

import junit.framework.TestCase;

public class FullTruthTupleDriverTest extends TestCase {
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
        
        loop2.add(new org.hps.recon.filtering.EventFlagFilter());
        
        FullTruthTupleTestDriver ftttd = new FullTruthTupleTestDriver();
        ftttd.setOutputPlots(aidaOutput);
        loop2.add(ftttd);

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
    
    protected class FullTruthTupleTestDriver extends Driver {
        public AIDA aida = null;
        private String outputPlots = null;
        private double ebeam = 0;
        private double bfield = 0;
        double tupleTrkPCut = 0.9;
        double tupleMaxSumCut = 1.3;
        private String V0CollectionName = "UnconstrainedV0Candidates";
        private String trackerHitsCollectionName = "";
        
        @Override
        public void endOfData() {
            if (outputPlots != null) {
                try {
                    aida.saveAs(outputPlots);
                } catch (IOException ex) {
                }
            }
        }
        
        public void setOutputPlots(String input) {
            outputPlots = input;
        }
        
        @Override
        protected void detectorChanged(Detector detector) {            
            if (aida == null)
                aida = AIDA.defaultInstance();
            aida.tree().cd("/");
            
            if (ebeam == 0) {
                try {
                    BeamEnergy.BeamEnergyCollection beamEnergyCollection = this.getConditionsManager().getCachedConditions(BeamEnergy.BeamEnergyCollection.class, "beam_energies").getCachedData();
                    ebeam = beamEnergyCollection.get(0).getBeamEnergy();
                } catch (Exception e) {
                }
            }
            
            this.bfield = Math.abs(TrackUtils.getBField(detector).y());
            
            setupPlots();
        }
        
        private void setupPlots() {
            aida.histogram1D("Purity", 20, 0, 1.0);
            aida.histogram1D("Tracking Strategy", 128, 0, 128);
            aida.histogram1D("Good TruthHit in layer", 13, 0, 13);
            for(int layer = 1; layer < 13; layer++)
                aida.histogram1D(String.format("NTruthParticles in layer %d", layer), 10, 0, 10);
        }
        
        @Override
        public void process(EventHeader event) {
            List<ReconstructedParticle> uncV0List = null;
            if (event.hasCollection(ReconstructedParticle.class, V0CollectionName)) {
                uncV0List = event.get(ReconstructedParticle.class, V0CollectionName);
            } else
                return;
            
            List<SimTrackerHit> trackerHits = null;
            if (event.hasCollection(SimTrackerHit.class, trackerHitsCollectionName))
                event.get(SimTrackerHit.class, trackerHitsCollectionName);
            else
                return;

            
            RelationalTable rawtomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
            if (event.hasCollection(LCRelation.class, "SVTTrueHitRelations")) {
                List<LCRelation> trueHitRelations = event.get(LCRelation.class, "SVTTrueHitRelations");
                for (LCRelation relation : trueHitRelations)
                    if (relation != null && relation.getFrom() != null && relation.getTo() != null)
                        rawtomc.add(relation.getFrom(), relation.getTo());
            }
            
            RelationalTable seedtogbl = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
            if (event.hasCollection(LCRelation.class, "MatchedToGBLTrackRelations")) {
                List<LCRelation> trueHitRelations = event.get(LCRelation.class, "MatchedToGBLTrackRelations");
                for (LCRelation relation : trueHitRelations)
                    if (relation != null && relation.getFrom() != null && relation.getTo() != null)
                        seedtogbl.add(relation.getFrom(), relation.getTo());
            }
            
            for (ReconstructedParticle uncV0 : uncV0List) {
                if (doRecoParticle(uncV0)) {
                    getPurity(uncV0.getParticles().get(0).getTracks().get(0), trackerHits, rawtomc);
                    getStrategy(uncV0.getParticles().get(0).getTracks().get(0), seedtogbl);
                }
            }
        }
        
        private void getStrategy(Track gblTrk, RelationalTable seedtogbl) {
            Set<Track> seedTracks = seedtogbl.allTo(gblTrk);
            for (Track seedTrk : seedTracks) {
                aida.histogram1D("Tracking Strategy").fill(seedTrk.getType());
            }
        }
        
        private void getPurity(Track trk, List<SimTrackerHit> trackerHits, RelationalTable rawtomc) {
            Map<MCParticle, List<SimTrackerHit>> trackerHitMap = MCFullDetectorTruth.BuildTrackerHitMap(trackerHits);
            
            TrackTruthMatching pTruth = new TrackTruthMatching(trk, rawtomc, trackerHits);

            aida.histogram1D("Purity").fill(pTruth.getPurity());
            for(int layer = 1; layer < 13; layer++){
                aida.histogram1D(String.format("NTruthParticles in layer %d", layer)).fill(pTruth.getNumberOfMCParticles(layer));
                if(pTruth.getHitList(layer) != null && pTruth.getHitList(layer))
                    aida.histogram1D("Good TruthHit in layer").fill(layer);
            }
        }
        
        
        private boolean doRecoParticle(ReconstructedParticle uncV0) {
            if (uncV0 == null)
                return false;
            ReconstructedParticle electron = uncV0.getParticles().get(0);
            ReconstructedParticle positron = uncV0.getParticles().get(1);
            if (electron == null || positron == null)
                return false;

            Track eleTrk = electron.getTracks().get(0);
            Track posTrk = positron.getTracks().get(0);
            if (eleTrk == null || posTrk == null)
                return false;

            double ele_pt = Math.abs((1.0 / eleTrk.getTrackStates().get(0).getOmega()) * bfield * 2.99792458e-04);
            double ele_px = ele_pt * Math.sin(eleTrk.getTrackStates().get(0).getPhi());
            double ele_py = ele_pt * eleTrk.getTrackStates().get(0).getTanLambda();
            double ele_pz = ele_pt * Math.cos(eleTrk.getTrackStates().get(0).getPhi());
            double ele_p = Math.sqrt(ele_px * ele_px + ele_py * ele_py + ele_pz * ele_pz);

            double pos_pt = Math.abs((1.0 / posTrk.getTrackStates().get(0).getOmega()) * bfield * 2.99792458e-04);
            double pos_px = pos_pt * Math.sin(posTrk.getTrackStates().get(0).getPhi());
            double pos_py = pos_pt * posTrk.getTrackStates().get(0).getTanLambda();
            double pos_pz = pos_pt * Math.cos(posTrk.getTrackStates().get(0).getPhi());
            double pos_p = Math.sqrt(pos_px * pos_px + pos_py * pos_py + pos_pz * pos_pz);

            if (ele_p > tupleTrkPCut * ebeam)
                return false;
            if (pos_p > tupleTrkPCut * ebeam)
                return false;

            if (ele_p + pos_p > tupleMaxSumCut * ebeam)
                return false;

            if (eleTrk.getTrackStates().get(0).getTanLambda() * posTrk.getTrackStates().get(0).getTanLambda() > 0)
                return false;

            if (ele_p > 0.25)
                return false;
            if (eleTrk.getTrackerHits().get(0).getPosition()[2] > 0)
                return false;

            return true;
            
        }
    }
}
