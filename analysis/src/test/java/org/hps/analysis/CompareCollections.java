package org.hps.analysis;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import junit.framework.TestCase;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.recon.tracking.TrackType;
import org.hps.recon.tracking.TrackUtils;
import org.hps.record.StandardCuts;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Track;
import org.lcsim.event.Vertex;
import org.lcsim.geometry.Detector;
import org.lcsim.util.aida.AIDA;
import org.lcsim.util.loop.LCSimLoop;

public class CompareCollections extends TestCase {
    String oldFileName = "/nfs/slac/g/hps3/data/objectStandardization/engrun2015/5772_pass8_goldenV0.slcio";
    //String oldFileName = "/nfs/slac/work/mdiamond/hps-java/outputPass8.slcio";
    String[] newFileName = {"/nfs/slac/work/mdiamond/hps-java/outputMouse.slcio"};
    int nEvents = 1000;

    public void testClear() throws Exception {
        File oldFile = new File(oldFileName);
        LCSimLoop loop = new LCSimLoop();
        
        loop.setLCIORecordSource(oldFile);
        CompareCollectionsDriver ccd = new CompareCollectionsDriver();
        ccd.setOverlayFiles(newFileName);
        loop.add(ccd);

        try {
            loop.loop(nEvents);
        } catch (Exception e) {
            System.out.println(e.toString());
        }
        loop.dispose();
    }
    
    protected class CompareCollectionsDriver extends org.lcsim.util.OverlayDriver {
        String unconstrainedV0CandidatesColName = "TargetConstrainedV0Candidates";
        private String outputPlots = "CompareCollPlots.aida";
        public AIDA aida = null;
        private double timeOffset = 43;
        StandardCuts cuts;
        
        @Override
        public void endOfData() {
            if (outputPlots != null) {
                try {
                    aida.saveAs(outputPlots);
                } catch (IOException ex) {
                    System.out.println("aida write error");
                    Logger.getLogger(CompareCollectionsDriver.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        
        @Override 
        protected void detectorChanged(Detector detector) {
            if (aida == null)
                aida = AIDA.defaultInstance();
            aida.tree().cd("/");
            setupPlots();

	    System.out.printf("Detector name %s \n", detector.getName());
            
            cuts = new StandardCuts();
        }
        
        protected void setV0CandidatesColName(String input) {
            unconstrainedV0CandidatesColName = input;
        }
        
        protected void process(EventHeader event) {
            EventHeader extraEvent = getNextEvent(overlayEvents);
            List<ReconstructedParticle> V0sNew = null;
            List<ReconstructedParticle> V0sOld = null;
            
            
            if (extraEvent == null) {
                System.out.printf("error: null extraEvent %d \n", event.getEventNumber());
                return;
            }
            //System.out.println("Got extra event");
            if (!extraEvent.hasCollection(ReconstructedParticle.class, unconstrainedV0CandidatesColName)) {
                System.out.printf("error: extraEvent has no V0 collection %d \n", event.getEventNumber());
                return;
            }
            V0sNew = extraEvent.get(ReconstructedParticle.class, unconstrainedV0CandidatesColName);

            if (!event.hasCollection(ReconstructedParticle.class, unconstrainedV0CandidatesColName)) {
                System.out.printf("error: event has no V0 collection %d \n", event.getEventNumber());
                return;
            }
            V0sOld = event.get(ReconstructedParticle.class, unconstrainedV0CandidatesColName);
            
            
            if (!V0sOld.isEmpty() && V0sNew.isEmpty()) {
                RelationalTable hitToRotatedTable = TrackUtils.getHitToRotatedTable(event);
                RelationalTable hitToStripsTable = TrackUtils.getHitToStripsTable(event);
                System.out.printf("new V0 empty: event num %d \n", event.getEventNumber());
                doRecoParticles(V0sOld, hitToStripsTable, hitToRotatedTable);
            }

        }
        
        private void setupPlots() {
            aida.histogram1D("Match GoodnessOfPID", 100, 0, 25);
            aida.histogram1D("Track-Cluster dt", 50, -10, 10);
            aida.histogram1D("V0 Cluster dt", 50, 0, 10);
            aida.histogram1D("Electron P", 100, 0, 1.5);
            aida.histogram1D("V0 P", 100, 0, 3.0);
            aida.histogram1D("V0 Chi2", 500, 0, 500);
            aida.histogram1D("Track Chi2", 500, 0, 500);
            aida.histogram1D("same half", 2, -0.5, 1.5);
            aida.histogram1D("cut flow", 10, -0.5, 9.5);
        }
        
        private void doRecoParticles(List<ReconstructedParticle> V0s, RelationalTable hitToStrips, RelationalTable hitToRotated) {

            
            if (V0s == null || V0s.isEmpty())
                return;
            
            for (ReconstructedParticle V0 : V0s) {
            // vertex
                boolean passes = true;
                
                Vertex v = V0.getStartVertex();
                if (v == null)
                    continue;
                aida.histogram1D("cut flow").fill(0);
                
                if (V0.getParticles().size() < 2)
                    continue;
                ReconstructedParticle p1 = V0.getParticles().get(0);
                ReconstructedParticle p2 = V0.getParticles().get(1);
                if (!TrackType.isGBL(p1.getType()))
                    continue;
                if (!TrackType.isGBL(p2.getType()))
                    continue;

                boolean eleIsTop = (p1.getTracks().get(0).getTrackerHits().get(0).getPosition()[2] > 0);
                boolean posIsTop = (p2.getTracks().get(0).getTrackerHits().get(0).getPosition()[2] > 0);
                boolean sameTop = (eleIsTop == posIsTop);
                if (sameTop) {
                    aida.histogram1D("same half").fill(1);
                    passes = false;
                }
                else {
                    aida.histogram1D("same half").fill(0);
                    aida.histogram1D("cut flow").fill(1);
                }

                
                aida.histogram1D("V0 P").fill(V0.getMomentum().magnitude());
                if (V0.getMomentum().magnitude() < cuts.getMaxVertexP()) {
                    if (passes)
                        aida.histogram1D("cut flow").fill(2);
                }
                else {
                    passes = false;
		    
		}

                // electron/positron components
                List<ReconstructedParticle> parts = V0.getParticles();
                //boolean passPID = true;
                boolean passChi2 = true;
                boolean passTiming = true;
		boolean passElectronP = true;
		double maxChi2 = 0;
		double electronP = 0;
		double maxDt = 0;
		double maxPID = 0;
		boolean passPID = true;
                if (parts != null && !parts.isEmpty()) {
                    for (ReconstructedParticle part : parts) {
                        if (part.getCharge() == -1) {
			    electronP = part.getMomentum().magnitude();
                            if (part.getMomentum().magnitude() > cuts.getMaxElectronP()) 
                                passElectronP = false;
                        }
                        if (part.getCharge() != 0) {
			    if (part.getGoodnessOfPID() > maxPID)
				maxPID = part.getGoodnessOfPID();
                            if (part.getGoodnessOfPID() > cuts.getMaxMatchChisq())
                                passPID = false;
                        }
                        if (!part.getTracks().isEmpty()) {
                            Track trk = part.getTracks().get(0);
                            if (trk != null) {

                                if (trk.getChi2() > cuts.getMaxTrackChisq(trk.getTrackerHits().size()))
                                    passChi2 = false;
			        if (trk.getChi2() > maxChi2)
				    maxChi2 = trk.getChi2();
                                if (!part.getClusters().isEmpty()) {
                                    Cluster clus = part.getClusters().get(0);
                                    double clusTime = ClusterUtilities.getSeedHitTime(clus);
                                    double trkT = TrackUtils.getTrackTime(trk, hitToStrips, hitToRotated);
                                    double temp = Math.abs(clusTime - trkT - timeOffset);
				    if (temp > maxDt)
					maxDt = temp;
                                    if (temp > cuts.getMaxMatchDt())
                                        passTiming = false;
                                }
				else
				    passes=false;
                            }
			    else
				passes=false;
                        }
			else
			    passes=false;
			     
		    }

		    boolean passVtiming = false;
		    double diff = 0;
		    if (!p1.getClusters().isEmpty() && !p2.getClusters().isEmpty()) {
			double clus1 = ClusterUtilities.getSeedHitTime(p1.getClusters().get(0));
			double clus2 = ClusterUtilities.getSeedHitTime(p2.getClusters().get(0));
			diff = Math.abs(clus1 - clus2);
                        if (diff < cuts.getMaxVertexClusterDt()) {
			    passVtiming = true;
                        }

		    }
		    else
			passes=false;

		    if (!passes)
			continue;

		    double vProb = 1.0 - new ChiSquaredDistribution(4).cumulativeProbability(v.getChi2());

		    aida.histogram1D("Track-Cluster dt").fill(maxDt);
		    if (passTiming) {
			aida.histogram1D("cut flow").fill(3);
			aida.histogram1D("Match GoodnessOfPID").fill(maxPID);
		    
			if (passPID) {
			    aida.histogram1D("cut flow").fill(4);
			    aida.histogram1D("Electron P").fill(electronP);
			    if (passElectronP) {
				aida.histogram1D("V0 Cluster dt").fill(diff);
				aida.histogram1D("cut flow").fill(5);
				if (passVtiming) {
				    aida.histogram1D("cut flow").fill(6);
				    aida.histogram1D("V0 Chi2").fill(v.getChi2());
				    if (vProb > cuts.getMinVertexChisqProb()) {
					aida.histogram1D("cut flow").fill(7);
					aida.histogram1D("Track Chi2").fill(maxChi2);
					if (passChi2)
					    aida.histogram1D("cut flow").fill(8);
				    }
				}
			    }
			}
		    }
                }
            }
        }
    }    

}
