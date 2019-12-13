package org.hps.recon.tracking;

//Configuration parser
import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.Path;
import java.nio.charset.StandardCharsets;

import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
//import hep.aida.IProfile;
import hep.physics.vec.Hep3Vector;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.LCIOParameters.ParameterName;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.event.GenericObject;
import org.lcsim.event.LCRelation;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.IDDecoder;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.lcsim.fit.helicaltrack.HelicalTrackFit;

//MCTruth
import org.lcsim.event.MCParticle;

/**
 * Analysis class to check recon.
 * 
 * @author phansson
 * @author mdiamond <mdiamond@slac.stanford.edu>
 * @author pf       <pbutti@slac.stanford.edu>
 */
public class TrackingReconstructionPlots extends Driver {

    //static {
    //    hep.aida.jfree.AnalysisFactory.register();
    //}

    public AIDA aida;
    private String helicalTrackHitCollectionName = "HelicalTrackHits";
    private String stripClusterCollectionName = "StripClusterer_SiTrackerHitStrip1D";
    private String plotConfig = "plotconfigs/trackingRecoPlots.xml";
    private boolean doAmplitudePlots = false;
    private boolean doECalClusterPlots = false;
    private boolean doHitsOnTrackPlots = false;
    private boolean doResidualPlots = false;
    private boolean doMatchedEcalClusterPlots = false;
    private boolean doElectronPositronPlots = false;
    private boolean doStripHitPlots = false;
    private boolean doMCTruthPlots  = true;
    private List<String> volumes = new ArrayList<String>();
    
    private String trackCollectionName = "GBLTracks";
    private String MCParticleCollectionName = "MCParticle";
    private String trackRelationCollectionName = "MatchedToGBLTrackRelations";
    
    
    String ecalSubdetectorName = "Ecal";
    String ecalCollectionName = "EcalClusters";
    IDDecoder dec;
    private Map<Track, Cluster> eCanditates;
    private Map<Track, Cluster> pCanditates;
    
    private String outputPlots = "TrackingRecoPlots.aida";
    private String trackSelection = "";
    
    
    ShaperFitAlgorithm _shaper = new DumbShaperFit();
    HelixConverter converter = new HelixConverter(0);
    private static Logger LOGGER = Logger.getLogger(TrackingReconstructionPlots.class.getName());
    private List<HpsSiSensor> sensors = new ArrayList<HpsSiSensor>();
    private double bfield;
    private double bfield_y;
    private int event_nr = 0;
    
    //Scaling factor from curvature to momentum
    private double momentum_param = 2.99792458e-04;

    @Override
    protected void detectorChanged(Detector detector) {
        if (aida == null)
            aida = AIDA.defaultInstance();

        aida.tree().cd("/");
        for (HpsSiSensor s : detector.getDetectorElement().findDescendants(HpsSiSensor.class)) {
            if (s.getName().startsWith("module_") && s.getName().endsWith("sensor0")) {
                sensors.add(s);
            }
        }
        LOGGER.info("Found " + sensors.size() + " SiSensors.");

        Hep3Vector fieldInTracker = TrackUtils.getBField(detector);
        this.bfield = Math.abs(fieldInTracker.y());
        this.bfield_y = fieldInTracker.y();
        volumes.add("");
        volumes.add("top_");
        volumes.add("bot_");
        setupPlots();
    }
    
    public TrackingReconstructionPlots() {
        LOGGER.setLevel(Level.WARNING);
    }
    
    public void setTrackSelection(String val) {
        trackSelection = val;
    }
    
    public void setPlotConfig(String val) {
        plotConfig = val;
    }

    public void setOutputPlots(String output) {
        this.outputPlots = output;
    }

    public void setHelicalTrackHitCollectionName(String helicalTrackHitCollectionName) {
        this.helicalTrackHitCollectionName = helicalTrackHitCollectionName;
    }

    public void setTrackCollectionName(String trackCollectionName) {
        this.trackCollectionName = trackCollectionName;
    }

    public void setDoAmplitudePlots(boolean value) {
        this.doAmplitudePlots = value;
    }

    public void setDoECalClusterPlots(boolean value) {
        this.doECalClusterPlots = value;
    }

    public void setDoHitsOnTrackPlots(boolean value) {
        this.doHitsOnTrackPlots = value;
    }

    public void setDoResidualPlots(boolean value) {
        this.doResidualPlots = value;
    }

    public void setDoMatchedEcalClusterPlots(boolean value) {
        this.doMatchedEcalClusterPlots = value;
    }

    public void setDoElectronPositronPlots(boolean value) {
        this.doElectronPositronPlots = value;
    }

    private void doStripHits(List<TrackerHit> stripClusters, Track trk, RelationalTable trackDataTable) {
        Map<HpsSiSensor, Integer> stripHits = new HashMap<HpsSiSensor, Integer>();

        for (TrackerHit stripHit : stripClusters) {
            HpsSiSensor sensor = ((HpsSiSensor) ((RawTrackerHit) stripHit.getRawHits().get(0)).getDetectorElement());

            int n;
            if (stripHits.containsKey(sensor)) {
                n = stripHits.get(sensor);
            } else {
                n = 0;
            }
            n++;
            stripHits.put(sensor, n);
        }

        for (Map.Entry<HpsSiSensor, Integer> sensor : stripHits.entrySet()) {
            aida.histogram1D(sensor.getKey().getName() + " strip hits").fill(stripHits.get(sensor.getKey()));
        }

        if (trackDataTable == null)
            return;
        GenericObject trackData = (GenericObject) trackDataTable.from(trk);
        if (trackData == null) {
            System.out.println("null TrackData for isolation");
            return;
        }

        int numIso = trackData.getNDouble();
        for (int i = 0; i < numIso; i++) {
            aida.histogram1D(String.format("Layer %d Isolation", i)).fill(trackData.getDoubleVal(i));
        }

    }

    private void doECalClusters(List<Cluster> clusters, boolean tracksPresent) {
        int nBotClusters = 0;
        int nTopClusters = 0;

        for (Cluster cluster : clusters) {
            // Get the ix and iy indices for the seed.
            //                final int ix = cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
            //                final int iy = cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");

            //System.out.println("cluser position = ("+cluster.getPosition()[0]+","+cluster.getPosition()[1]+") with energy = "+cluster.getEnergy());
            if (cluster.getPosition()[1] > 0) {
                nTopClusters++;
                //System.out.println("cl " + cluster.getPosition()[0] + " " + cluster.getPosition()[1] + "  ix  " + ix + " iy " + iy);
                aida.histogram2D("Top ECal Cluster Position").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                aida.histogram1D("Top ECal Cluster Energy").fill(cluster.getEnergy());
            }
            if (cluster.getPosition()[1] < 0) {
                nBotClusters++;
                aida.histogram2D("Bottom ECal Cluster Position").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                aida.histogram1D("Bottom ECal Cluster Energy").fill(cluster.getEnergy());
            }

            if (tracksPresent) {
                if (cluster.getPosition()[1] > 0) {
                    aida.histogram2D("Top ECal Cluster Position (>0 tracks)").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                }
                if (cluster.getPosition()[1] < 0) {
                    aida.histogram2D("Bottom ECal Cluster Position (>0 tracks)").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                }

                if (cluster.getEnergy() > 0.1) {
                    if (cluster.getPosition()[1] > 0) {
                        aida.histogram2D("Top ECal Cluster Position (E>0.1,>0 tracks)").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                        aida.histogram2D("Top ECal Cluster Position w_E (E>0.1,>0 tracks)").fill(cluster.getPosition()[0], cluster.getPosition()[1], cluster.getEnergy());
                    }
                    if (cluster.getPosition()[1] < 0) {
                        aida.histogram2D("Bottom ECal Cluster Position (E>0.1,>0 tracks)").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                        aida.histogram2D("Bottom ECal Cluster Position w_E (E>0.1,>0 tracks)").fill(cluster.getPosition()[0], cluster.getPosition()[1], cluster.getEnergy());
                    }
                }
            }

        }

        aida.histogram1D("Number of Clusters Top").fill(nTopClusters);
        aida.histogram1D("Number of Clusters Bot").fill(nBotClusters);
    }
    
    
    private void doMCTruth(List<Track> tracks,List<MCParticle> mcParticles){
        doMCTruth(tracks,mcParticles,"");
    }
    
    private void doMCTruth(List<Track> tracks,List<MCParticle> mcParticles,String name) { 
        
        Map <Track, MCParticle> trackTruthMatch  = new HashMap<Track, MCParticle>();
        Map <MCParticle,Track> duplicatedTracks = new HashMap<MCParticle,Track>();
        List<Track> fakeTracks       = new ArrayList<Track>();
        //List<Track> duplicatedTracks = new ArrayList<Track>();
        
        aida.histogram1D("n_mcParticles").fill(mcParticles.size());
        //Do the track matching
        int track_index   = 0;
        int matched_index = 0;
        int fake_index    = 0;
        for (Track track : tracks) {
            MCParticle part = TrackUtils.getMatchedTruthParticle(track);
            if (part != null) {
                //Check for duplicates
                if (trackTruthMatch.containsValue(part)) {
                    duplicatedTracks.put(part,track); // this holds the duplicated tracks only.
                }//Check for duplicates
                else {
                    trackTruthMatch.put(track,part); // uniquely matched tracks.
                }
            }
            else {
                System.out.println("fake track");
                fakeTracks.add(track);        
            }       
        }

        //Plot for the truth matched tracks
        //Plot for the fake  tracks 
        //NO plots for fake tracks at the moment
        //doBasicTracks(fakeTracks,"fakes");
        
        
        //Only on "truth matched tracks" (makes sense??)
        for (Track track : trackTruthMatch.keySet()) {
            
            MCParticle part = trackTruthMatch.get(track);
            //Get the MCparticle HTF
            HelicalTrackFit pHTF = null;
            double pTruth = -1.;
            double pTrackTruth = -1.;
            double pcharge = 0;
            pTruth = part.getMomentum().magnitude();
            pcharge = part.getCharge();
            pHTF = TrackUtils.getHTF(part, Math.abs(bfield_y));
            if (pHTF == null) 
                continue;
            
            pTrackTruth = pHTF.p(Math.abs(bfield_y));
            
            //Corrected for the first GBL Point
            double d0Gbl = track.getTrackStates().get(0).getD0();
            double z0Gbl = track.getTrackStates().get(0).getZ0();
            double CGbl = track.getTrackStates().get(0).getOmega();
            double phiGbl = track.getTrackStates().get(0).getPhi();
            double slopeGbl = track.getTrackStates().get(0).getTanLambda();
            double pGbl = getMag(track.getTrackStates().get(0).getMomentum());
            double chargeGbl = track.getCharge();
            
            //get the track errors
            
            //      d0  phi0 C    z0   tanL
            //d0    0  
            //phi0  1   2
            //C     3   4    5
            //z0    6   7    8    9
            //tanL  10  11   12   13   14 
            double[] covMatrix = track.getTrackStates().get(1).getCovMatrix();
            
            double d0Gbl_err = Math.sqrt(covMatrix[0]);
            double phiGbl_err = Math.sqrt(covMatrix[2]);
            double CGbl_err = Math.sqrt(covMatrix[5]);
            double z0Gbl_err = Math.sqrt(covMatrix[9]);
            double slopeGbl_err = Math.sqrt(covMatrix[14]);
            double pGbl_err = (Math.abs(Math.pow(1/CGbl,2)) * CGbl_err)*Math.abs(bfield)*momentum_param;
            
            
            
            double d0Truth = pHTF.dca();
            double z0Truth = pHTF.z0();
            double CTruth = pHTF.curvature();
            double phiTruth = pHTF.phi0();
            double slopeTruth = pHTF.slope();
            //double chargeTruth = 0. //pHTF.charge();
            
            //Fill truth info
            //aida.histogram1D(name+"truth_part_p").fill(pTruth);
            //aida.histogram1D(name+"truth_part_q").fill(pcharge);
            //aida.histogram1D(name+"truth_trk_p").fill(pTrackTruth);
            ////aida.histogram1D("truth_trk_q").fill(chargeTruth);
            //aida.histogram1D(name+"truth_trk_d0").fill(d0Truth);
            //aida.histogram1D(name+"truth_trk_phi").fill(phiTruth);
            //aida.histogram1D(name+"truth_trk_omega").fill(CTruth);
            //aida.histogram1D(name+"truth_trk_tanLambda").fill(slopeTruth);
            //aida.histogram1D("truth_trk_omega").fill(chargeTruth);
        
            
            for (String volume : volumes) {
                
                aida.histogram1D(name+volume+"reco_truth_p").fill(pTrackTruth - pGbl);
                aida.histogram1D(name+volume+"reco_truth_d0").fill(d0Truth - d0Gbl);
                aida.histogram1D(name+volume+"reco_truth_z0").fill(z0Truth - z0Gbl);
                aida.histogram1D(name+volume+"reco_truth_phi").fill(phiTruth - phiGbl);
                aida.histogram1D(name+volume+"reco_truth_omega").fill(CTruth - CGbl);
                aida.histogram1D(name+volume+"reco_truth_tanLambda").fill(slopeTruth - slopeGbl);
                //aida.histogram1D(volume+"reco_truth_charge").fill(chargeTruth - chargeGbl);
                
                aida.histogram1D(name+volume+"pull_p").fill((pTrackTruth - pGbl) / pGbl_err);
                aida.histogram1D(name+volume+"pull_d0").fill((d0Truth - d0Gbl)/d0Gbl_err);
                aida.histogram1D(name+volume+"pull_z0").fill((z0Truth - z0Gbl)/z0Gbl_err);
                aida.histogram1D(name+volume+"pull_phi").fill((phiTruth - phiGbl)/phiGbl_err);
                aida.histogram1D(name+volume+"pull_omega").fill((CTruth - CGbl)/CGbl_err);
                aida.histogram1D(name+volume+"pull_tanLambda").fill((slopeTruth - slopeGbl)/slopeGbl_err);

                aida.histogram1D(name+volume+"err_p").fill(pGbl_err);
                aida.histogram1D(name+volume+"err_d0").fill(d0Gbl_err);
                aida.histogram1D(name+volume+"err_z0").fill(z0Gbl_err);
                aida.histogram1D(name+volume+"err_phi").fill(phiGbl_err);
                aida.histogram1D(name+volume+"err_omega").fill(CGbl_err);
                aida.histogram1D(name+volume+"err_tanLambda").fill(slopeGbl_err);
            }
        }
    }
    
    //Argh - terrible
    private double getMag(double p[]) {
        return Math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    }
    
    private void doBasicTracks(List<Track> tracks) {
        doBasicTracks(tracks,"");
    }
    
        
    private void doBasicTracks(List<Track> tracks, String name) {
        int ntracksTop = 0;
        int ntracksBot = 0;
        
        aida.histogram1D(name+"n_tracks").fill(tracks.size());
        
        
        for (Track trk : tracks) {
            
            boolean isTop = false;
            int volume_index = 0;
            
            if (trk.getTrackerHits().get(0).getPosition()[2] > 0) {
                isTop = true;
            }
            
            double pt = Math.abs((1 / trk.getTrackStates().get(0).getOmega()) * bfield * momentum_param);
            
            double pz = pt * Math.cos(trk.getTrackStates().get(0).getPhi());
            double px = pt * Math.sin(trk.getTrackStates().get(0).getPhi());
            double py = pt * trk.getTrackStates().get(0).getTanLambda();

            double d0      = trk.getTrackStates().get(0).getD0();
            double d0prime = trk.getTrackStates().get(0).getParameter(ParameterName.d0.ordinal());
            double z0      = trk.getTrackStates().get(0).getZ0();
            double C       = trk.getTrackStates().get(0).getOmega();
            double phi     = trk.getTrackStates().get(0).getPhi();
            double slope   = trk.getTrackStates().get(0).getTanLambda();
            
            //TODO use this instead the strings
            if (slope > 0)
                volume_index = 1;
            else 
                volume_index = 2;

            aida.histogram1D(name+"pz").fill(pz);
            aida.histogram1D(name+"py").fill(py);
            aida.histogram1D(name+"px").fill(px);
            aida.histogram1D(name+"p").fill(pt);
            aida.histogram1D(name+"chi2").fill(trk.getChi2());
            aida.histogram1D(name+"chi2ndof").fill(trk.getChi2() / trk.getNDF());

            double[] covMatrix = trk.getTrackStates().get(0).getCovMatrix();
            aida.histogram1D(name+"err_p").fill((Math.pow(1/C,2)*Math.sqrt(covMatrix[5]))*Math.abs(bfield)*momentum_param);
            aida.histogram1D(name+"err_d0").fill(Math.sqrt(covMatrix[0]));
            aida.histogram1D(name+"err_z0").fill(Math.sqrt(covMatrix[9]));
            aida.histogram1D(name+"err_phi").fill(Math.sqrt(covMatrix[2]));
            aida.histogram1D(name+"err_omega").fill(Math.sqrt(covMatrix[5]));
            aida.histogram1D(name+"err_tanLambda").fill(Math.sqrt(covMatrix[14]));
            
            /*
            aida.histogram1D(name+"d0").fill(trk.getTrackStates().get(0).getParameter(ParameterName.d0.ordinal()));
            aida.histogram1D(name+"phi").fill(trk.getTrackStates().get(0).getParameter(ParameterName.phi0.ordinal()));
            aida.histogram1D(name+"omega").fill(trk.getTrackStates().get(0).getParameter(ParameterName.omega.ordinal()));
            aida.histogram1D(name+"tanLambda").fill(trk.getTrackStates().get(0).getParameter(ParameterName.tanLambda.ordinal()));
            aida.histogram1D(name+"z0").fill(trk.getTrackStates().get(0).getParameter(ParameterName.z0.ordinal()));
            */
            
            aida.histogram1D(name+"d0").fill(d0);
            aida.histogram1D(name+"phi").fill(phi);
            aida.histogram1D(name+"omega").fill(C);
            aida.histogram1D(name+"tanLambda").fill(slope);
            aida.histogram1D(name+"z0").fill(z0);
            aida.histogram1D(name+"HitsOnTrack").fill(trk.getTrackerHits().size());
            
            if (isTop) {
                
                aida.histogram1D(name+"top_chi2").fill(trk.getChi2());
                aida.histogram1D(name+"top_chi2ndof").fill(trk.getChi2() / trk.getNDF());
                aida.histogram1D(name+"top_HitsOnTrack").fill(trk.getTrackerHits().size());
                aida.histogram1D(name+"top_px").fill(px);
                aida.histogram1D(name+"top_py").fill(py);
                aida.histogram1D(name+"top_pz").fill(pz);
                aida.histogram1D(name+"top_p").fill(pt);
                aida.histogram1D(name+"top_chi2").fill(trk.getChi2());
                aida.histogram1D(name+"top_chi2").fill(trk.getChi2() / (double)trk.getNDF());
                aida.histogram1D(name+"top_d0").fill(trk.getTrackStates().get(0).getParameter(ParameterName.d0.ordinal()));
                aida.histogram1D(name+"top_phi").fill(trk.getTrackStates().get(0).getParameter(ParameterName.phi0.ordinal()));
                aida.histogram1D(name+"top_omega").fill(trk.getTrackStates().get(0).getParameter(ParameterName.omega.ordinal()));
                aida.histogram1D(name+"top_tanLambda").fill(trk.getTrackStates().get(0).getParameter(ParameterName.tanLambda.ordinal()));
                aida.histogram1D(name+"top_z0").fill(trk.getTrackStates().get(0).getParameter(ParameterName.z0.ordinal()));
                ntracksTop++;
            } else {
                aida.histogram1D(name+"bot_chi2").fill(trk.getChi2());
                aida.histogram1D(name+"bot_chi2ndof").fill(trk.getChi2() / trk.getNDF());
                aida.histogram1D(name+"bot_HitsOnTrack").fill(trk.getTrackerHits().size());
                aida.histogram1D(name+"bot_px").fill(px);
                aida.histogram1D(name+"bot_py").fill(py);
                aida.histogram1D(name+"bot_pz").fill(pz);
                aida.histogram1D(name+"bot_p").fill(pt);
                aida.histogram1D(name+"bot_chi2").fill(trk.getChi2());
                aida.histogram1D(name+"bot_chi2").fill(trk.getChi2()/(double)trk.getNDF());
                aida.histogram1D(name+"bot_d0").fill(trk.getTrackStates().get(0).getParameter(ParameterName.d0.ordinal()));
                aida.histogram1D(name+"bot_phi").fill(trk.getTrackStates().get(0).getParameter(ParameterName.phi0.ordinal()));
                aida.histogram1D(name+"bot_omega").fill(trk.getTrackStates().get(0).getParameter(ParameterName.omega.ordinal()));
                aida.histogram1D(name+"bot_tanLambda").fill(trk.getTrackStates().get(0).getParameter(ParameterName.tanLambda.ordinal()));
                aida.histogram1D(name+"bot_z0").fill(trk.getTrackStates().get(0).getParameter(ParameterName.z0.ordinal()));
                ntracksBot++;
            }
        }

        aida.histogram1D(name+"bot_n_tracks").fill(ntracksBot);
        aida.histogram1D(name+"top_n_tracks").fill(ntracksTop);
    }
    
    private void doHitsOnTrack(Track trk) {
        doHitsOnTrack(trk,"");
    }
    
    private void doHitsOnTrack(Track trk, String name) {
        Map<HpsSiSensor, Integer> stripHitsOnTrack = new HashMap<HpsSiSensor, Integer>();
        List<TrackerHit> hitsOnTrack = trk.getTrackerHits();
        
        boolean isTop = true;
        if (trk.getTrackStates().get(0).getTanLambda() < 0) {
            isTop = false;
        }

        for (TrackerHit hit : hitsOnTrack) {
            
            HpsSiSensor sensor = ((HpsSiSensor) ((RawTrackerHit) hit.getRawHits().get(0)).getDetectorElement());
            //Ly0-Ly6
            int layer = (sensor.getLayerNumber()+1)/2;
            aida.histogram1D(name+"hitDistribution").fill(layer);
            if (isTop)
                aida.histogram1D(name+"top_hitDistribution").fill(layer);
            else
                aida.histogram1D(name+"bot_hitDistribution").fill(layer);
            
            if (stripHitsOnTrack.containsKey(sensor)) {
                stripHitsOnTrack.put(sensor, stripHitsOnTrack.get(sensor) + 1);
            } else {
                stripHitsOnTrack.put(sensor, 1);
            }
        }
        boolean doDetailed = false;
        
        if (doDetailed) {
            for (Map.Entry<HpsSiSensor, Integer> sensor : stripHitsOnTrack.entrySet()) {
                aida.histogram1D(sensor.getKey().getName() + " strip hits on track").fill(stripHitsOnTrack.get(sensor.getKey()));
            }
        }
    }

    private void doResiduals(List<LCRelation> fittedHits, Track trk, RelationalTable trackResTable) {
        GenericObject trackRes = (GenericObject) trackResTable.from(trk);
        if (trackRes == null) {
            //System.out.println("null TrackResidualsData");
            return;
        }

        int numX = trackRes.getNDouble();
        for (int i = 0; i < numX; i++) {
            int layer = trackRes.getIntVal(i);
            String modNum = "Layer Unknown ";
            if (layer % 2 == 1)
                modNum = String.format("Layer %d ", layer / 2 + 1);
            aida.histogram1D(modNum + "Residual Y(mm)").fill(trackRes.getFloatVal(i));
            aida.histogram1D(modNum + "Residual X(mm)").fill(trackRes.getDoubleVal(i));
        }
    }

    private void doAmplitude(List<LCRelation> fittedHits, Track trk) {
        List<TrackerHit> hitsOnTrack = trk.getTrackerHits();

        for (TrackerHit hit : hitsOnTrack) {
            double clusterSum = 0;
            double clusterT0 = 0;
            int nHitsCluster = 0;

            for (RawTrackerHit rawHit : (List<RawTrackerHit>) hit.getRawHits()) {

                for (LCRelation fittedHit : fittedHits) {
                    if (rawHit.equals((RawTrackerHit) fittedHit.getFrom())) {
                        double amp = FittedRawTrackerHit.getAmp(fittedHit);
                        double t0 = FittedRawTrackerHit.getT0(fittedHit);
                        //System.out.println("to="+t0 + " amp=" + amp);
                        aida.histogram1D("Amp (HitOnTrack)").fill(amp);
                        if (trk.getTrackStates().get(0).getMomentum()[0] > 0.8) {
                            aida.histogram1D("Amp Pz>0.8 (HitOnTrack)").fill(amp);
                        }
                        aida.histogram1D("t0 (HitOnTrack)").fill(t0);
                        if (trk.getTrackStates().get(0).getMomentum()[0] > 0.8) {
                            aida.histogram1D("t0 Pz>0.8 (HitOnTrack)").fill(t0);
                        }
                        clusterSum += amp;
                        clusterT0 += t0;
                        nHitsCluster++;
                    }
                }
            }

            aida.histogram1D("Hits in Cluster (HitOnTrack)").fill(nHitsCluster);
            aida.histogram1D("Cluster Amp (HitOnTrack)").fill(clusterSum);
            if (trk.getTrackStates().get(0).getMomentum()[0] > 0.8) {
                aida.histogram1D("Cluster Amp Pz>0.8 (HitOnTrack)").fill(clusterSum);
            }
            if (nHitsCluster > 0) {
                aida.histogram1D("Cluster t0 (HitOnTrack)").fill(clusterT0 / nHitsCluster);
                if (trk.getTrackStates().get(0).getMomentum()[0] > 0.8) {
                    aida.histogram1D("Cluster t0 Pz>0.8 (HitOnTrack)").fill(clusterT0 / nHitsCluster);
                }
            }

        }
    }

    private void doEcalClustersOnTrack(Track trk, List<Cluster> clusters) {
        Hep3Vector posAtEcal = TrackUtils.getTrackPositionAtEcal(trk);
        Cluster clust = findClosestCluster(posAtEcal, clusters);
        if (clust == null)
            return;

        boolean isTop = false;
        if (trk.getTrackerHits().get(0).getPosition()[2] > 0) {
            isTop = true;
        }

        // track matching requirement
        if (Math.abs(posAtEcal.x() - clust.getPosition()[0]) < 30.0 && Math.abs(posAtEcal.y() - clust.getPosition()[1]) < 30.0) {

            if (doElectronPositronPlots) {
                if (trk.getCharge() < 0)
                    pCanditates.put(trk, clust);
                else
                    eCanditates.put(trk, clust);
            }

            posAtEcal = TrackUtils.extrapolateTrack(trk, clust.getPosition()[2]);//.positionAtEcal();

            aida.histogram2D("Energy Vs Momentum").fill(clust.getEnergy(), trk.getTrackStates().get(0).getMomentum()[0]);
            aida.histogram1D("Energy Over Momentum").fill(clust.getEnergy() / (trk.getTrackStates().get(0).getMomentum()[0]));
            aida.histogram1D("deltaX").fill(clust.getPosition()[0] - posAtEcal.x());
            aida.histogram1D("deltaY").fill(clust.getPosition()[1] - posAtEcal.y());
            aida.histogram2D("X ECal Vs Track").fill(clust.getPosition()[0], posAtEcal.x());
            aida.histogram2D("Y ECal Vs Track").fill(clust.getPosition()[1], posAtEcal.y());

            if (isTop) {
                aida.histogram2D("Top Energy Vs Momentum").fill(clust.getEnergy(), trk.getTrackStates().get(0).getMomentum()[0]);
                //                    aida.histogram2D("Top Energy Vs Momentum").fill(posAtEcal.y(), trk.getTrackStates().get(0).getMomentum()[0]);
                aida.histogram1D("Top Energy Over Momentum").fill(clust.getEnergy() / (trk.getTrackStates().get(0).getMomentum()[0]));
                aida.histogram1D("Top deltaX").fill(clust.getPosition()[0] - posAtEcal.x());
                aida.histogram1D("Top deltaY").fill(clust.getPosition()[1] - posAtEcal.y());
                aida.histogram2D("Top deltaX vs X").fill(clust.getPosition()[0], clust.getPosition()[0] - posAtEcal.x());
                aida.histogram2D("Top deltaY vs Y").fill(clust.getPosition()[1], clust.getPosition()[1] - posAtEcal.y());
                aida.histogram2D("Top X ECal Vs Track").fill(clust.getPosition()[0], posAtEcal.x());
                aida.histogram2D("Top Y ECal Vs Track").fill(clust.getPosition()[1], posAtEcal.y());
            } else {
                aida.histogram2D("Bottom Energy Vs Momentum").fill(clust.getEnergy(), trk.getTrackStates().get(0).getMomentum()[0]);
                aida.histogram1D("Bottom Energy Over Momentum").fill(clust.getEnergy() / (trk.getTrackStates().get(0).getMomentum()[0]));
                aida.histogram1D("Bottom deltaX").fill(clust.getPosition()[0] - posAtEcal.x());
                aida.histogram1D("Bottom deltaY").fill(clust.getPosition()[1] - posAtEcal.y());
                aida.histogram2D("Bottom deltaX vs X").fill(clust.getPosition()[0], clust.getPosition()[0] - posAtEcal.x());
                aida.histogram2D("Bottom deltaY vs Y").fill(clust.getPosition()[1], clust.getPosition()[1] - posAtEcal.y());
                aida.histogram2D("Bottom X ECal Vs Track").fill(clust.getPosition()[0], posAtEcal.x());
                aida.histogram2D("Bottom Y ECal Vs Track").fill(clust.getPosition()[1], posAtEcal.y());
            }

            aida.histogram1D("Tracks matched").fill(0);
            if (trk.getTrackStates().get(0).getMomentum()[0] > 0.8) {
                aida.histogram1D("Tracks matched (Pz>0.8)").fill(0);
            }
            if (isTop) {
                aida.histogram1D("Tracks matched Top").fill(0);
                if (trk.getTrackStates().get(0).getMomentum()[0] > 0.8) {
                    aida.histogram1D("Tracks matched Top (Pz>0.8)").fill(0);
                }
            } else {
                aida.histogram1D("Tracks matched Bottom").fill(0);
                if (trk.getTrackStates().get(0).getMomentum()[0] > 0.8) {
                    aida.histogram1D("Tracks matched Bottom (Pz>0.8)").fill(0);
                }
            }
        }

        else {
            aida.histogram1D("Tracks matched").fill(1);
            if (trk.getTrackStates().get(0).getMomentum()[0] > 0.8) {
                aida.histogram1D("Tracks matched (Pz>0.8)").fill(1);
            }

            if (isTop) {
                aida.histogram1D("Tracks matched Top").fill(1);
                if (trk.getTrackStates().get(0).getMomentum()[0] > 0.8) {
                    aida.histogram1D("Tracks matched Top (Pz>0.8)").fill(1);
                }
            } else {
                aida.histogram1D("Tracks matched Bottom").fill(1);
                if (trk.getTrackStates().get(0).getMomentum()[0] > 0.8) {
                    aida.histogram1D("Tracks matched Bottom (Pz>0.8)").fill(1);
                }
            }
        }
    }

    @Override
    public void process(EventHeader event) {
        if (event_nr % 100 == 0)
            System.out.println("Processed " + event_nr + " events");
        
            
        aida.tree().cd("/");
        if (!event.hasCollection(TrackerHit.class, helicalTrackHitCollectionName)) {
            System.out.println(helicalTrackHitCollectionName + " does not exist; skipping event");
            return;
        }

        if (!event.hasCollection(Track.class, trackCollectionName)) {
            System.out.println(trackCollectionName + " does not exist; skipping event");
            aida.histogram1D("Number Tracks/Event").fill(0);
            return;
        }
        List<Track> all_tracks = event.get(Track.class, trackCollectionName);      
        
        RelationalTable trackMatchTable = null;
        trackMatchTable = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);
        List<LCRelation> trackMatchRelation = event.get(LCRelation.class, trackRelationCollectionName);
        for (LCRelation relation : trackMatchRelation) {
            if (relation != null && relation.getFrom() != null && relation.getTo() != null) {
                trackMatchTable.add(relation.getFrom(), relation.getTo());
            }
        }

        //Select tracks => This should be improved
        List<Track> tracks = new ArrayList<Track>();
        TrackSelectorTool trackSelector = new TrackSelectorTool();
        
        //TODO do a proper mapping
        //Exactly 5 Hits and missing hit not in ly6, enhances the MS effects
        if (trackSelection.equals("hits5_L6")) {  
            trackSelector.setMinL6hits(1);
            trackSelector.setNHits(5);
        }
        
        for (Track track : all_tracks) { 
            if (trackSelector.pass(track)) {
                tracks.add(track);
            }
        }
        
        //Get the MCParticles
        List<MCParticle> mc_Particles = null;
        
        if (event.hasCollection(MCParticle.class,MCParticleCollectionName)) {
            
            mc_Particles = event.get(MCParticle.class,MCParticleCollectionName);
            if (doMCTruthPlots) {
                doMCTruth(tracks,mc_Particles);
            }
        }//Missing MCParticle collection
        
        
        List<Cluster> clusters = null;
        if (event.hasCollection(Cluster.class, ecalCollectionName)) {
            clusters = event.get(Cluster.class, ecalCollectionName);
        } else {
            doECalClusterPlots = false;
            doMatchedEcalClusterPlots = false;
            doElectronPositronPlots = false;
        }

        List<LCRelation> fittedHits = null;
        if (event.hasCollection(LCRelation.class, "SVTFittedRawTrackerHits")) {
            fittedHits = event.get(LCRelation.class, "SVTFittedRawTrackerHits");
        } else {
            doAmplitudePlots = false;
            doResidualPlots = false;
        }

        List<TrackerHit> stripClusters = null;
        if (event.hasCollection(TrackerHit.class, stripClusterCollectionName)) {
            stripClusters = event.get(TrackerHit.class, stripClusterCollectionName);
        } else {
            doStripHitPlots = false;
        }

        //RelationalTable hitToRotatedTable = TrackUtils.getHitToRotatedTable(event);
        //RelationalTable hitToStripsTable = TrackUtils.getHitToStripsTable(event);

        RelationalTable trackResidualsTable = null;
        if (event.hasCollection(LCRelation.class, "TrackResidualsRelations")) {
            trackResidualsTable = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);
            List<LCRelation> trackresRelation = event.get(LCRelation.class, "TrackResidualsRelations");
            for (LCRelation relation : trackresRelation) {
                if (relation != null && relation.getFrom() != null && relation.getTo() != null) {
                    trackResidualsTable.add(relation.getFrom(), relation.getTo());
                }
            }
        } else {
            doResidualPlots = false;
        }
        RelationalTable trackDataTable = null;
        if (event.hasCollection(LCRelation.class, "TrackDataRelations")) {
            trackDataTable = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);
            List<LCRelation> trackdataRelation = event.get(LCRelation.class, "TrackDataRelations");
            for (LCRelation relation : trackdataRelation) {
                if (relation != null && relation.getFrom() != null && relation.getTo() != null) {
                    trackDataTable.add(relation.getFrom(), relation.getTo());
                }
            }
        }
        
        doBasicTracks(tracks);

        List<Track> MatchedTracks = new ArrayList<Track>();
      
        for (Track trk : tracks)  {
            MatchedTracks.add((Track) trackMatchTable.from(trk));
        }
        
        doBasicTracks(MatchedTracks,"matchedTrack");
        //Switch Matched Tracks off
        //doMCTruth(MatchedTracks,mc_Particles,"matchedTrack");
        

        if (doECalClusterPlots)
            doECalClusters(clusters, tracks.size() > 0);

        if (doElectronPositronPlots) {
            eCanditates = new HashMap<Track, Cluster>();
            pCanditates = new HashMap<Track, Cluster>();
        }

        for (Track trk : tracks) {
            if (doStripHitPlots)
                doStripHits(stripClusters, trk, trackDataTable);

            if (doHitsOnTrackPlots)
                doHitsOnTrack(trk);

            if (doResidualPlots)
                doResiduals(fittedHits, trk, trackResidualsTable);

            if (doAmplitudePlots)
                doAmplitude(fittedHits, trk);

            if (doMatchedEcalClusterPlots)
                doEcalClustersOnTrack(trk, clusters);
        }

        if (doElectronPositronPlots)
            doElectronPositron();
        
        event_nr++;
        
    }

    private void doElectronPositron() {

        Map.Entry<Track, Cluster> ecand_highestP = null;
        double e_pmax = -1;
        Map.Entry<Track, Cluster> pcand_highestP = null;
        double p_pmax = -1;
        for (Map.Entry<Track, Cluster> ecand : eCanditates.entrySet()) {
            double p = getMomentum(ecand.getKey());
            aida.histogram1D("p(e-)").fill(p);
            if (ecand_highestP == null) {
                ecand_highestP = ecand;
                e_pmax = getMomentum(ecand_highestP.getKey());
            } else {
                if (p > e_pmax) {
                    ecand_highestP = ecand;
                    e_pmax = getMomentum(ecand_highestP.getKey());
                }
            }
        }

        for (Map.Entry<Track, Cluster> pcand : pCanditates.entrySet()) {
            double p = getMomentum(pcand.getKey());
            aida.histogram1D("p(e+)").fill(p);
            if (pcand_highestP == null) {
                pcand_highestP = pcand;
                p_pmax = getMomentum(pcand_highestP.getKey());
            } else {
                if (p > p_pmax) {
                    pcand_highestP = pcand;
                    p_pmax = getMomentum(pcand_highestP.getKey());
                }
            }
        }

        aida.histogram1D("n(e-)").fill(eCanditates.size());
        aida.histogram1D("n(e+)").fill(pCanditates.size());
        if (ecand_highestP != null) {
            aida.histogram1D("p(e-) max").fill(e_pmax);
        }
        if (pcand_highestP != null) {
            aida.histogram1D("p(e+) max").fill(p_pmax);
        }
        if (ecand_highestP != null && pcand_highestP != null) {
            aida.histogram2D("p(e-) vs p(e+) max").fill(e_pmax, p_pmax);
        }
    }

    private double getMomentum(Track trk) {
        double p = Math.sqrt(trk.getTrackStates().get(0).getMomentum()[0] * trk.getTrackStates().get(0).getMomentum()[0] + trk.getTrackStates().get(0).getMomentum()[1] * trk.getTrackStates().get(0).getMomentum()[1] + trk.getTrackStates().get(0).getMomentum()[2] * trk.getTrackStates().get(0).getMomentum()[2]);
        return p;
    }

    private Cluster findClosestCluster(Hep3Vector posonhelix, List<Cluster> clusters) {
        Cluster closest = null;
        double minDist = 9999;
        for (Cluster cluster : clusters) {
            double[] clPos = cluster.getPosition();
            double clEne = cluster.getEnergy();
            double dist = Math.sqrt(Math.pow(clPos[0] - posonhelix.x(), 2) + Math.pow(clPos[1] - posonhelix.y(), 2)); //coordinates!!!
            if (dist < minDist && clEne > 0.4) {
                closest = cluster;
                minDist = dist;
            }
        }
        return closest;
    }

    @Override
    public void endOfData() {
        if (outputPlots != null) {
            try {
                aida.saveAs(outputPlots);
            } catch (IOException ex) {
                Logger.getLogger(TrackingReconstructionPlots.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
    private void setupPlotsFromCfgFile() {
        Path path = Paths.get(plotConfig);
        try {
            BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8);
            while (true) {
                String line = reader.readLine();
                if (line == null) {
                    break;
                }
                String[] settings = line.split(" ");
                //skip "#" commented lines
                if (line.indexOf('#') > -1 )
                    continue;
                System.out.println("settings[0]=" + settings[0]);
                if (settings[0].equals("h1d")) {
                    System.out.println("doing histo " + settings[0]);
                    int    xbins= new Integer(settings[2]);
                    double xmin = new Double(settings[3]);
                    double xmax = new Double(settings[4]);
                    aida.histogram1D(settings[1],xbins,xmin,xmax);
                }
                if (settings[0].equals("h2d")) {
                    int    xbins= new Integer(settings[2]);
                    double xmin = new Double(settings[3]);
                    double xmax = new Double(settings[4]);
                    int    ybins= new Integer(settings[6]);
                    double ymin = new Double(settings[7]);
                    double ymax = new Double(settings[8]);
                    aida.histogram2D(settings[1],xbins,xmin,xmax,ybins,ymin,ymax);
                }
                System.out.println(line);
            }
        }
        
        catch(IOException ex) {
            System.out.println (ex.toString());
            System.out.println ("TrackingReconstructionPlots::Couldn't find file " + plotConfig);
        }
    }
    
    
    private void setupTruthPlots() {
        
        aida.histogram1D("n_mcParticles",10,0,10);
        aida.histogram1D("truth_part_p",100,0,4);
        aida.histogram1D("truth_part_q",4,-2,2);
        aida.histogram1D("truth_trk_p",100,0,4);
        aida.histogram1D("truth_trk_d0",100,-10.0,10.0);
        aida.histogram1D("truth_trk_phi",100,-0.2,0.2);
        aida.histogram1D("truth_trk_omega",100,-0.2,0.2);
        aida.histogram1D("truth_trk_tanLambda",100,-0.001,0.001);
        aida.histogram1D("truth_trk_charge",4,-2,2);
        
        for (String volume : volumes) {
            //Truth Residuals
            aida.histogram1D(volume+"reco_truth_p",100,-0.5,0.5);
            aida.histogram1D(volume+"reco_truth_d0",100,-0.5,0.5);
            aida.histogram1D(volume+"reco_truth_z0",100,-0.5,0.5);
            aida.histogram1D(volume+"reco_truth_phi",100,-0.05,0.05);
            aida.histogram1D(volume+"reco_truth_omega",100,-0.0001,0.0001);
            aida.histogram1D(volume+"reco_truth_tanLambda",100,-0.01,0.01);
            
            //Pulls
            aida.histogram1D(volume+"pull_p",100,-5,5);
            aida.histogram1D(volume+"pull_d0",100,-5,5);
            aida.histogram1D(volume+"pull_z0",100,-5,5);
            aida.histogram1D(volume+"pull_phi",100,-5,5);
            aida.histogram1D(volume+"pull_omega",100,-5.,5);
            aida.histogram1D(volume+"pull_tanLambda",100,-5,5);
        }
    }

    private void setupBasicTracksPlots(String name) {
        // Basic tracks
        
        for (String volume : volumes) {
            
            aida.histogram1D(name+volume+"px", 100, -0.15, 0.15);
            aida.histogram1D(name+volume+"py", 100, -0.15, 0.15);
            aida.histogram1D(name+volume+"pz", 100, 0, 1.5);
            aida.histogram1D(name+volume+"p",  100, 0, 4);
            aida.histogram1D(name+volume+"chi2", 25, 0, 25.0);
            aida.histogram1D(name+volume+"chi2ndf", 50, 0, 10.0);
            aida.histogram1D(name+volume+"d0", 100, -10.0, 10.0);
            aida.histogram1D(name+volume+"phi", 100, -0.2, 0.2);
            aida.histogram1D(name+volume+"omega", 100, -0.001, 0.001);
            aida.histogram1D(name+volume+"tanLambda", 100, -0.1, 0.1);
            aida.histogram1D(name+volume+"z0", 100, -4.0, 4.0);
            aida.histogram1D(name+volume+"n_tracks", 10, 0, 10);
            aida.histogram1D(name+volume+"HitsOnTrack", 4, 3, 7);
            
            //Track errors
            aida.histogram1D(name+volume+"err_p",  100, 0, 4);;
            aida.histogram1D(name+volume+"err_d0", 100, -10.0, 10.0);
            aida.histogram1D(name+volume+"err_phi", 100, -0.2, 0.2);
            aida.histogram1D(name+volume+"err_omega", 100, -0.001, 0.001);
            aida.histogram1D(name+volume+"err_tanLambda", 100, -0.1, 0.1);
            aida.histogram1D(name+volume+"err_z0", 100, -4.0, 4.0);
        }
    }
    
    private void setupPlots() {
        
        setupPlotsFromCfgFile();
        
        //setupBasicTracksPlots("");
        
        // if (doMCTruthPlots) {
        //  setupTruthPlots();
        //  setupBasicTracksPlots("fakes");
            //duplicates!
        //}
        
        if (doStripHitPlots) {
            int i = 0;
            for (SiSensor sensor : sensors) {
                IHistogram1D resX = aida.histogram1D(sensor.getName() + " strip hits", 10, 0, 10);
                i++;
            }

            for (i = 0; i < 12; i++) {
                IHistogram1D resX = aida.histogram1D(String.format("Layer %d Isolation", i), 50, 0, 5);
            }
        }

        if (doAmplitudePlots) {
            IHistogram1D nHitsCluster = aida.histogram1D("Hits in Cluster (HitOnTrack)", 4, 0, 4);

            IHistogram1D amp = aida.histogram1D("Amp (HitOnTrack)", 50, 0, 5000);
            IHistogram1D ampcl = aida.histogram1D("Cluster Amp (HitOnTrack)", 50, 0, 5000);
            IHistogram1D amp2 = aida.histogram1D("Amp Pz>0.8 (HitOnTrack)", 50, 0, 5000);
            IHistogram1D ampcl2 = aida.histogram1D("Cluster Amp Pz>0.8 (HitOnTrack)", 50, 0, 5000);

            IHistogram1D t0 = aida.histogram1D("t0 (HitOnTrack)", 50, -100, 100);
            IHistogram1D t0cl = aida.histogram1D("Cluster t0 (HitOnTrack)", 50, -100, 100);
            IHistogram1D t02 = aida.histogram1D("t0 Pz>0.8 (HitOnTrack)", 50, -100, 100);
            IHistogram1D t0cl2 = aida.histogram1D("Cluster t0 Pz>0.8 (HitOnTrack)", 50, -100, 100);
        }

        if (doResidualPlots) {

            IHistogram1D mod1ResX = aida.histogram1D("Layer 1 Residual X(mm)", 25, -1, 1);
            IHistogram1D mod1ResY = aida.histogram1D("Layer 1 Residual Y(mm)", 25, -0.04, 0.04);

            IHistogram1D mod2ResX = aida.histogram1D("Layer 2 Residual X(mm)", 25, -2, 2);
            IHistogram1D mod2ResY = aida.histogram1D("Layer 2 Residual Y(mm)", 25, -1, 1);

            IHistogram1D mod3ResX = aida.histogram1D("Layer 3 Residual X(mm)", 25, -2.5, 2.5);
            IHistogram1D mod3ResY = aida.histogram1D("Layer 3 Residual Y(mm)", 25, -1.5, 1.5);

            IHistogram1D mod4ResX = aida.histogram1D("Layer 4 Residual X(mm)", 25, -3.0, 3.0);
            IHistogram1D mod4ResY = aida.histogram1D("Layer 4 Residual Y(mm)", 25, -2, 2);

            IHistogram1D mod5ResX = aida.histogram1D("Layer 5 Residual X(mm)", 25, -4, 4);
            IHistogram1D mod5ResY = aida.histogram1D("Layer 5 Residual Y(mm)", 25, -3, 3);

            IHistogram1D mod6ResX = aida.histogram1D("Layer 6 Residual X(mm)", 25, -5, 5);
            IHistogram1D mod6ResY = aida.histogram1D("Layer 6 Residual Y(mm)", 25, -3, 3);
        }

        if (doMatchedEcalClusterPlots) {

            IHistogram2D eVsP = aida.histogram2D("Energy Vs Momentum", 50, 0, 0.50, 50, 0, 1.5);
            IHistogram1D eOverP = aida.histogram1D("Energy Over Momentum", 50, 0, 2);

            IHistogram1D distX = aida.histogram1D("deltaX", 50, -100, 100);
            IHistogram1D distY = aida.histogram1D("deltaY", 50, -40, 40);

            IHistogram2D xEcalVsTrk = aida.histogram2D("X ECal Vs Track", 100, -400, 400, 100, -400, 400);
            IHistogram2D yEcalVsTrk = aida.histogram2D("Y ECal Vs Track", 100, -100, 100, 100, -100, 100);

            IHistogram2D topeVsP = aida.histogram2D("Top Energy Vs Momentum", 50, 0, 0.500, 50, 0, 1.5);
            IHistogram1D topeOverP = aida.histogram1D("Top Energy Over Momentum", 50, 0, 2);

            IHistogram1D topdistX = aida.histogram1D("Top deltaX", 50, -100, 100);
            IHistogram1D topdistY = aida.histogram1D("Top deltaY", 50, -40, 40);

            IHistogram2D topxEcalVsTrk = aida.histogram2D("Top X ECal Vs Track", 100, -400, 400, 100, -100, 100);
            IHistogram2D topyEcalVsTrk = aida.histogram2D("Top Y ECal Vs Track", 100, 0, 100, 100, 0, 100);

            IHistogram2D BottomeVsP = aida.histogram2D("Bottom Energy Vs Momentum", 50, 0, 0.500, 50, 0, 1.5);
            IHistogram1D BottomeOverP = aida.histogram1D("Bottom Energy Over Momentum", 50, 0, 2);

            IHistogram1D BottomdistX = aida.histogram1D("Bottom deltaX", 50, -100, 100);
            IHistogram1D BottomdistY = aida.histogram1D("Bottom deltaY", 50, -40, 40);

            IHistogram2D BottomxEcalVsTrk = aida.histogram2D("Bottom X ECal Vs Track", 100, -400, 400, 100, -400, 400);
            IHistogram2D BottomyEcalVsTrk = aida.histogram2D("Bottom Y ECal Vs Track", 100, -100, 0, 100, -100, 0);

            IHistogram2D topdistXvsX = aida.histogram2D("Top deltaX vs X", 51, -400, 400, 25, -100, 100);
            IHistogram2D topdistYvsY = aida.histogram2D("Top deltaY vs Y", 51, 0, 100, 25, -40, 40);

            IHistogram2D botdistXvsX = aida.histogram2D("Bottom deltaX vs X", 51, -400, 400, 25, -100, 100);
            IHistogram2D botdistYvsY = aida.histogram2D("Bottom deltaY vs Y", 51, -100, 0, 25, -40, 40);

            IHistogram1D trackmatchN = aida.histogram1D("Tracks matched", 3, 0, 3);
            IHistogram1D toptrackmatchN = aida.histogram1D("Tracks matched Top", 3, 0, 3);
            IHistogram1D bottrackmatchN = aida.histogram1D("Tracks matched Bottom", 3, 0, 3);
            IHistogram1D trackmatchN2 = aida.histogram1D("Tracks matched (Pz>0.8)", 3, 0, 3);
            IHistogram1D toptrackmatchN2 = aida.histogram1D("Tracks matched Top (Pz>0.8)", 3, 0, 3);
            IHistogram1D bottrackmatchN2 = aida.histogram1D("Tracks matched Bottom (Pz>0.8)", 3, 0, 3);
        }

        if (doElectronPositronPlots) {
            IHistogram2D trackPCorr = aida.histogram2D("p(e-) vs p(e+) max", 25, 0, 1.2, 25, 0, 1.2);
            IHistogram1D ne = aida.histogram1D("n(e-)", 3, 0, 3);
            IHistogram1D np = aida.histogram1D("n(e+)", 3, 0, 3);
            IHistogram1D pem = aida.histogram1D("p(e-) max", 25, 0, 1.5);
            IHistogram1D pe = aida.histogram1D("p(e-)", 25, 0, 1.5);
            IHistogram1D ppm = aida.histogram1D("p(e+) max", 25, 0, 1.5);
            IHistogram1D pp = aida.histogram1D("p(e+)", 25, 0, 1.5);
        }

        if (doECalClusterPlots) {
            IHistogram2D topECal = aida.histogram2D("Top ECal Cluster Position", 50, -400, 400, 10, 0, 100);
            IHistogram2D botECal = aida.histogram2D("Bottom ECal Cluster Position", 50, -400, 400, 10, -100, 0);
            IHistogram2D topECal1 = aida.histogram2D("Top ECal Cluster Position (>0 tracks)", 50, -400, 400, 10, 0, 100);
            IHistogram2D botECal1 = aida.histogram2D("Bottom ECal Cluster Position (>0 tracks)", 50, -400, 400, 10, -100, 0);
            IHistogram2D topECal2 = aida.histogram2D("Top ECal Cluster Position (E>0.1,>0 tracks)", 50, -400, 400, 10, 0, 100);
            IHistogram2D botECal2 = aida.histogram2D("Bottom ECal Cluster Position (E>0.1,>0 tracks)", 50, -400, 400, 10, -100, 0);
            IHistogram2D topECal3 = aida.histogram2D("Top ECal Cluster Position w_E (E>0.1,>0 tracks)", 50, -400, 400, 10, 0, 100);
            IHistogram2D botECal3 = aida.histogram2D("Bottom ECal Cluster Position w_E (E>0.1,>0 tracks)", 50, -400, 400, 10, -100, 0);

            IHistogram1D topECalE = aida.histogram1D("Top ECal Cluster Energy", 50, 0, 2);
            IHistogram1D botECalE = aida.histogram1D("Bottom ECal Cluster Energy", 50, 0, 2);
            IHistogram1D topECalN = aida.histogram1D("Number of Clusters Top", 6, 0, 6);
            IHistogram1D botECalN = aida.histogram1D("Number of Clusters Bot", 6, 0, 6);

        }

        if (doHitsOnTrackPlots) {
            int i = 0;
            for (SiSensor sensor : sensors) {
                IHistogram1D resX = aida.histogram1D(sensor.getName() + " strip hits on track", 50, 0, 5);
                i++;
            }
        }
    }
}



