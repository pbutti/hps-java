package org.hps.analysis.MC;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.detector.svt.SvtDetectorSetup;
import org.hps.recon.tracking.CoordinateTransformations;
import org.hps.recon.tracking.TrackUtils;
import org.lcsim.event.EventHeader;
import org.lcsim.event.MCParticle;
//import org.lcsim.event.SimCalorimeterHit;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.FieldMap;
import org.lcsim.recon.tracking.digitization.sisim.config.RawTrackerHitSensorSetup;
import org.lcsim.recon.tracking.digitization.sisim.config.ReadoutCleanupDriver;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.lcsim.util.loop.LCSimLoop;
import org.lcsim.util.swim.Trajectory;

import junit.framework.TestCase;

public class TrackExtrapTest extends TestCase {

    public void testIt() throws Exception
    {
        int nEvents = 10;
        String fileName = "RecoCopy_tritrig-wab-beam.slcio";
        File inputFile = new File(fileName);
        String aidaOutput = fileName.replaceAll("slcio", "aida");
        LCSimLoop loop2 = new LCSimLoop();
        
        final DatabaseConditionsManager manager = new DatabaseConditionsManager();
        manager.addConditionsListener(new SvtDetectorSetup());

        loop2.setLCIORecordSource(inputFile);

        RawTrackerHitSensorSetup rthss = new RawTrackerHitSensorSetup();
        String[] readoutColl = { "SVTRawTrackerHits" };
        rthss.setReadoutCollections(readoutColl);
        loop2.add(rthss);

        TrackExtrapTestDriver trp = new TrackExtrapTestDriver();
        trp.setOutputPlots(aidaOutput);
        loop2.add(trp);

        ReadoutCleanupDriver rcd = new ReadoutCleanupDriver();
        loop2.add(rcd);
        
        loop2.loop(nEvents);
    }
    
    protected class TrackExtrapTestDriver extends Driver {
        private static final String ECAL_POSITION_CONSTANT_NAME = "ecal_dface";
        FieldMap bFieldMap = null;
        private double ecalPosition = 1338.0;
        private double stepSize = 5.0;
        private double epsilon = 0.05;
        public AIDA aida = null;
        private String outputPlots = null;
        private double lowMomThresh = 0.3;
        
        @Override
        public void endOfData() {
            if (outputPlots != null) {
                try {
                    aida.saveAs(outputPlots);
                } catch (IOException ex) {
                    Logger.getLogger(TrackExtrapTest.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        
        public void setOutputPlots(String input) {
            outputPlots = input;
        }
        
        public void setStepSize(double input) {
            stepSize=input;
        }
        public void setEpsilon(double input) {
            epsilon=input;
        }
        
        @Override
        protected void detectorChanged(Detector detector) {            
            bFieldMap = detector.getFieldMap();
            if (ecalPosition == 0)
                ecalPosition = detector.getConstants().get(ECAL_POSITION_CONSTANT_NAME).getValue();
            if (aida == null)
                aida = AIDA.defaultInstance();
            aida.tree().cd("/");
            setupPlots();
        }

        
        private void setupPlots() {
            aida.histogram2D("[Extrapolated SimTrackerHit - SimCalorimeterHit] Y vs SimCalorimeterHit Y", 100, -100, 100, 100, -100, 100);
            aida.histogram2D("[Extrapolated SimTrackerHit - SimCalorimeterHit] X vs SimCalorimeterHit Y", 100, -100, 100, 300, -300, 300);
            aida.histogram2D("[Extrapolated SimTrackerHit - SimCalorimeterHit] Y vs SimCalorimeterHit X", 100, -300, 300, 100, -100, 100);
            aida.histogram2D("[Extrapolated SimTrackerHit - SimCalorimeterHit] X vs SimCalorimeterHit X", 100, -300, 300, 300, -300, 300);
            
            aida.histogram2D("Extrapolated SimTrackerHit - SimCalorimeterHit: Y vs X", 300, -300, 300, 100, -100, 100);
            aida.histogram2D("Extrapolated SimTrackerHit - SimCalorimeterHit Y vs MCParticle Pz", 100, 0, 1.5, 100, -100, 100);
            aida.histogram2D("Extrapolated SimTrackerHit - SimCalorimeterHit X vs MCParticle Pz", 100, 0, 1.5, 300, -300, 300);
            
            aida.histogram2D("Top: Extrapolated SimTrackerHit - SimCalorimeterHit: Y vs X", 300, -300, 300, 100, -100, 100);
            aida.histogram2D("Top: Extrapolated SimTrackerHit - SimCalorimeterHit Y vs MCParticle Pz", 100, 0, 1.5, 100, -100, 100);
            aida.histogram2D("Top: Extrapolated SimTrackerHit - SimCalorimeterHit X vs MCParticle Pz", 100, 0, 1.5, 300, -300, 300);
            
            aida.histogram2D("Bot: Extrapolated SimTrackerHit - SimCalorimeterHit: Y vs X", 300, -300, 300, 100, -100, 100);
            aida.histogram2D("Bot: Extrapolated SimTrackerHit - SimCalorimeterHit Y vs MCParticle Pz", 100, 0, 1.5, 100, -100, 100);
            aida.histogram2D("Bot: Extrapolated SimTrackerHit - SimCalorimeterHit X vs MCParticle Pz", 100, 0, 1.5, 300, -300, 300);
            
            aida.histogram2D("e-: Extrapolated SimTrackerHit - SimCalorimeterHit: Y vs X", 300, -300, 300, 100, -100, 100);
            aida.histogram2D("e-: Extrapolated SimTrackerHit - SimCalorimeterHit Y vs MCParticle Pz", 100, 0, 1.5, 100, -100, 100);
            aida.histogram2D("e-: Extrapolated SimTrackerHit - SimCalorimeterHit X vs MCParticle Pz", 100, 0, 1.5, 300, -300, 300);
            
            aida.histogram2D("e+: Extrapolated SimTrackerHit - SimCalorimeterHit: Y vs X", 300, -300, 300, 100, -100, 100);
            aida.histogram2D("e+: Extrapolated SimTrackerHit - SimCalorimeterHit Y vs MCParticle Pz", 100, 0, 1.5, 100, -100, 100);
            aida.histogram2D("e+: Extrapolated SimTrackerHit - SimCalorimeterHit X vs MCParticle Pz", 100, 0, 1.5, 300, -300, 300);
        }
        
        @Override
        protected void process(EventHeader event) {
            List<SimTrackerHit> trackerHits = null;
            if (event.hasCollection(SimTrackerHit.class, "TrackerHits"))
                trackerHits = event.get(SimTrackerHit.class, "TrackerHits");
            else
                return;
//            List<SimCalorimeterHit> calHits = null;
//            if (event.hasCollection(SimCalorimeterHit.class, "EcalHits"))
//                calHits = event.get(SimCalorimeterHit.class, "EcalHits");
//            else
//                return;
            List<MCParticle> particles = null;
            if (event.hasCollection(MCParticle.class, "MCParticle"))
                particles = event.get(MCParticle.class, "MCParticle");  
            else
                return;
            List<SimTrackerHit> scoringHits = null;
            if (event.hasCollection(SimTrackerHit.class, "TrackerHitsECal"))
                scoringHits = event.get(SimTrackerHit.class, "TrackerHitsECal");
            else
                return;


            //System.out.println("analyzing event");
            
            Map<MCParticle, SimTrackerHit> scoringHitMap = new HashMap<MCParticle, SimTrackerHit>();
            for (SimTrackerHit scoringHit : scoringHits) {
                SimTrackerHit keyHit = scoringHitMap.get(scoringHit.getMCParticle());

                if (keyHit == null) {
                    scoringHitMap.put(scoringHit.getMCParticle(), scoringHit);
                } else if (scoringHit.getTime() < keyHit.getTime()) { //keep only the earliest hit from each particle
                    //System.out.println("Multiple scoring hits from same particle");
                    scoringHitMap.put(scoringHit.getMCParticle(), scoringHit);
                }
            }
            
            
            Map<MCParticle, List<SimTrackerHit>> trackerHitMap = MCFullDetectorTruth.BuildTrackerHitMap(trackerHits);
            //Map<MCParticle, List<SimCalorimeterHit>> calHitMap = MCFullDetectorTruth.BuildCalHitMap(calHits);
            
            for (MCParticle part : particles) {
                //System.out.printf("MCParticle charge %f Pz %f \n", part.getCharge(), Math.abs(part.getPZ()));
                if (part.getCharge() == 0)
                    continue;
                if (!trackerHitMap.containsKey(part) || (!scoringHitMap.containsKey(part)))
                    continue;
                List<SimTrackerHit> hits = trackerHitMap.get(part);
                //List<SimCalorimeterHit> caloHits = calHitMap.get(part);
                SimTrackerHit caloHit = scoringHitMap.get(part); 
                if ((hits == null) || (caloHit == null))
                    continue;
                
                //Hep3Vector startPosition = part.getOrigin();
                //Hep3Vector startMomentum = part.getMomentum();

                SimTrackerHit lastHit = null;
                int lay=0;
                for (SimTrackerHit hit : hits) {
                    if (hit.getLayer() > lay) {
                        lay=hit.getLayer();
                        lastHit = hit;
                    }
                }
                if (lastHit == null)
                    continue;
                    
                Hep3Vector hitPosition = lastHit.getPositionVec();
                Hep3Vector hitMomentum = new BasicHep3Vector(lastHit.getMomentum());
                Hep3Vector extrapPos = doTrackExtrap(CoordinateTransformations.transformVectorToTracking(hitPosition), CoordinateTransformations.transformVectorToTracking(hitMomentum), part.getCharge(), false);
                if (extrapPos.y() == Double.NaN) 
                    doTrackExtrap(CoordinateTransformations.transformVectorToTracking(hitPosition), CoordinateTransformations.transformVectorToTracking(hitMomentum), part.getCharge(), true);
                extrapPos=CoordinateTransformations.transformVectorToDetector(extrapPos);
                //System.out.printf("extrapPos %s \n", extrapPos.toString());
                
                //double residual=10000;
                //Hep3Vector residualVec = null;
                //for (SimCalorimeterHit caloHit : caloHits) {
                Hep3Vector caloHitPos = caloHit.getPositionVec();
                //System.out.printf("      caloHitPos %s \n", caloHitPos.toString());
                Hep3Vector residualVec = VecOp.sub(extrapPos, caloHitPos);
                    //double dist = Math.sqrt(diff.x() * diff.x() + diff.y() * diff.y());
//                    if (dist < residual) {
//                        residual =dist;
//                        residualVec = diff;
//                    }
                //}
                
                //System.out.printf("Pz %f  CaloHitPos %s  extrapPos %s \n", Math.abs(part.getPZ()), caloHitPos.toString(), extrapPos.toString());
                
                // filling plots
                if (Math.abs(part.getPZ()) < lowMomThresh) {
                    aida.histogram2D("Extrapolated SimTrackerHit - SimCalorimeterHit: Y vs X").fill(residualVec.x(), residualVec.y());
                    aida.histogram2D("[Extrapolated SimTrackerHit - SimCalorimeterHit] X vs SimCalorimeterHit X").fill(caloHitPos.x(), residualVec.x());
                    aida.histogram2D("[Extrapolated SimTrackerHit - SimCalorimeterHit] X vs SimCalorimeterHit Y").fill(caloHitPos.y(), residualVec.x());
                    aida.histogram2D("[Extrapolated SimTrackerHit - SimCalorimeterHit] Y vs SimCalorimeterHit X").fill(caloHitPos.x(), residualVec.y());
                    aida.histogram2D("[Extrapolated SimTrackerHit - SimCalorimeterHit] Y vs SimCalorimeterHit Y").fill(caloHitPos.y(), residualVec.y());
                }
                aida.histogram2D("Extrapolated SimTrackerHit - SimCalorimeterHit Y vs MCParticle Pz").fill(Math.abs(part.getPZ()), residualVec.y());
                aida.histogram2D("Extrapolated SimTrackerHit - SimCalorimeterHit X vs MCParticle Pz").fill(Math.abs(part.getPZ()), residualVec.x());
                System.out.printf("residual %s  caloHitPos %s  partPZ %f \n", residualVec.toString(), caloHitPos.toString(), part.getPZ());
              
                if (lastHit.getPosition()[1] > 0) {
                    if (Math.abs(part.getPZ()) < lowMomThresh)
                        aida.histogram2D("Top: Extrapolated SimTrackerHit - SimCalorimeterHit: Y vs X").fill(residualVec.x(), residualVec.y());
                    aida.histogram2D("Top: Extrapolated SimTrackerHit - SimCalorimeterHit Y vs MCParticle Pz").fill(Math.abs(part.getPZ()), residualVec.y());
                    aida.histogram2D("Top: Extrapolated SimTrackerHit - SimCalorimeterHit X vs MCParticle Pz").fill(Math.abs(part.getPZ()), residualVec.x());
                }
                else {
                    if (Math.abs(part.getPZ()) < lowMomThresh)
                        aida.histogram2D("Bot: Extrapolated SimTrackerHit - SimCalorimeterHit: Y vs X").fill(residualVec.x(), residualVec.y());
                    aida.histogram2D("Bot: Extrapolated SimTrackerHit - SimCalorimeterHit Y vs MCParticle Pz").fill(Math.abs(part.getPZ()), residualVec.y());
                    aida.histogram2D("Bot: Extrapolated SimTrackerHit - SimCalorimeterHit X vs MCParticle Pz").fill(Math.abs(part.getPZ()), residualVec.x());
                }
                
                if (part.getCharge() > 0) {
                    if (Math.abs(part.getPZ()) < lowMomThresh)
                        aida.histogram2D("e+: Extrapolated SimTrackerHit - SimCalorimeterHit: Y vs X").fill(residualVec.x(), residualVec.y());
                    aida.histogram2D("e+: Extrapolated SimTrackerHit - SimCalorimeterHit Y vs MCParticle Pz").fill(Math.abs(part.getPZ()), residualVec.y());
                    aida.histogram2D("e+: Extrapolated SimTrackerHit - SimCalorimeterHit X vs MCParticle Pz").fill(Math.abs(part.getPZ()), residualVec.x());
                }
                else {
                    if (Math.abs(part.getPZ()) < lowMomThresh)
                        aida.histogram2D("e-: Extrapolated SimTrackerHit - SimCalorimeterHit: Y vs X").fill(residualVec.x(), residualVec.y());
                    aida.histogram2D("e-: Extrapolated SimTrackerHit - SimCalorimeterHit Y vs MCParticle Pz").fill(Math.abs(part.getPZ()), residualVec.y());
                    aida.histogram2D("e-: Extrapolated SimTrackerHit - SimCalorimeterHit X vs MCParticle Pz").fill(Math.abs(part.getPZ()), residualVec.x());
                }
                
                    //Hep3Vector endPositionrot = VecOp.mult(beamAxisRotation, endPosition);
            }
        }
        
        // everything in tracking frame
        public Hep3Vector doTrackExtrap(Hep3Vector currentPosition, Hep3Vector currentMomentum, double q, boolean debug) {

            
            double bFieldY = bFieldMap.getField(new BasicHep3Vector(0, 0, 500)).y();
            if (bFieldY < 0)
                q = q * (-1);
            
            Hep3Vector currentPositionDet = null;

            double distance = ecalPosition - currentPosition.x();
            
            if (debug)
                System.out.printf("Track extrap, field %f distance %f pos %s mom %s \n", bFieldY, distance, currentPosition.toString(), currentMomentum.toString());
            
            
            if (stepSize == 0)
                stepSize = distance / 100.0;
            double sign = Math.signum(distance);
            distance = Math.abs(distance);

            while (distance > epsilon) {
                // The field map coordinates are in the detector frame so the
                // extrapolated track position needs to be transformed from the
                // track frame to detector.
                currentPositionDet = CoordinateTransformations.transformVectorToDetector(currentPosition);

                // Get the field at the current position along the track.
                bFieldY = bFieldMap.getField(currentPositionDet).y();

                // Get a trajectory (Helix or Line objects) created with the
                // track parameters at the current position.
                Trajectory trajectory = TrackUtils.getTrajectory(currentMomentum, new org.lcsim.spacegeom.SpacePoint(currentPosition), q, bFieldY);

                // Using the new trajectory, extrapolated the track by a step and
                // update the extrapolated position.
                Hep3Vector currentPositionTry = trajectory.getPointAtDistance(stepSize);
                if (debug)
                    System.out.printf("  > currentPositionTry %s \n", currentPositionTry.toString());

                if ((Math.abs(ecalPosition - currentPositionTry.x()) > epsilon) && (Math.signum(ecalPosition - currentPositionTry.x()) != sign)) {
                    // went too far, try again with smaller step-size
                    if (Math.abs(stepSize) > 0.01) {
                        stepSize /= 2.0;
                        if (debug)
                            System.out.println(" >>went too far");
                        continue;
                    } else {
                        if (debug)
                            System.out.println(" >>had to break");
                        break;
                    }
                }
                currentPosition = currentPositionTry;

                distance = Math.abs(ecalPosition - currentPosition.x());
                // Calculate the momentum vector at the new position. This will
                // be used when creating the trajectory that will be used to
                // extrapolate the track in the next iteration.
                currentMomentum = VecOp.mult(currentMomentum.magnitude(), trajectory.getUnitTangentAtLength(stepSize));
                if (debug)
                    System.out.printf("  > new distance %f mom %s \n", distance, currentMomentum.toString());
            }
            return currentPosition;
        }
    }
    
}
