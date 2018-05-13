package org.hps.recon.tracking;

import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.physics.vec.BasicHep3Vector;
//import hep.aida.IProfile;
import hep.physics.vec.Hep3Vector;
//import hep.physics.vec.VecOp;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.recon.tracking.SimpleAmbiguityResolver.AmbiMode;
import org.hps.record.triggerbank.AbstractIntData;
import org.hps.record.triggerbank.SSPCluster;
import org.hps.record.triggerbank.SSPData;
import org.hps.record.triggerbank.TIData;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.LCIOParameters.ParameterName;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.event.GenericObject;
import org.lcsim.event.LCRelation;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.event.TrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.IDDecoder;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 * Analysis class to check recon.
 * 
 * @author phansson
 * @author mdiamond <mdiamond@slac.stanford.edu>
 * @version $id: 2.0 06/04/17$
 */
public class TrackingReconstructionPlots extends Driver {

    //static {
    //    hep.aida.jfree.AnalysisFactory.register();
    //}

    public AIDA aida;
    private String helicalTrackHitCollectionName = "HelicalTrackHits";
    private String stripClusterCollectionName = "StripClusterer_SiTrackerHitStrip1D";
    private boolean doAmplitudePlots = false;
    private boolean doECalClusterPlots = false;
    private boolean doHitsOnTrackPlots = false;
    private boolean doResidualPlots = false;
    private boolean doMatchedClusterPlots = false;
    private boolean doElectronPositronPlots = false;
    private boolean doStripHitPlots = false;
    private boolean doReconParticlePlots = false;
    private boolean doOccupancyPlots = false;
    private boolean doComparisonPlots = true;

    private String trackCollectionName = "GBLTracks";
    String ecalSubdetectorName = "Ecal";
    String ecalCollectionName = "EcalClusters";
    IDDecoder dec;
    private Map<Track, Cluster> eCanditates;
    private Map<Track, Cluster> pCanditates;

    private String outputPlots = "TrackingRecoPlots.aida";

    ShaperFitAlgorithm _shaper = new DumbShaperFit();
    HelixConverter converter = new HelixConverter(0);
    private static Logger LOGGER = Logger.getLogger(TrackingReconstructionPlots.class.getName());
    private List<HpsSiSensor> sensors = new ArrayList<HpsSiSensor>();
    private double bfield;

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

        setupPlots();
    }

    public TrackingReconstructionPlots() {
        LOGGER.setLevel(Level.WARNING);
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

    public void setDoMatchedClusterPlots(boolean value) {
        this.doMatchedClusterPlots = value;
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

    private void doECalClusters(List<Cluster> clusters, boolean isHighEsum, boolean isLowEsum, boolean hasTopTrack, boolean hasBottomTrack) {
        //int nBotClusters = 0;
        //int nTopClusters = 0;

        for (Cluster cluster : clusters) {
            // Get the ix and iy indices for the seed.
            //                final int ix = cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix");
            //                final int iy = cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy");

            //System.out.println("cluser position = ("+cluster.getPosition()[0]+","+cluster.getPosition()[1]+") with energy = "+cluster.getEnergy());
            if (isHighEsum) {
                //nTopClusters++;
                //System.out.println("cl " + cluster.getPosition()[0] + " " + cluster.getPosition()[1] + "  ix  " + ix + " iy " + iy);
                if (hasTopTrack && !hasBottomTrack)
                    aida.histogram2D("HighEsum TopTrack ECal Cluster Position").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                else if (!hasTopTrack && hasBottomTrack)
                    aida.histogram2D("HighEsum BottomTrack ECal Cluster Position").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                //aida.histogram1D("Top ECal Cluster Energy").fill(cluster.getEnergy());
            } else if (isLowEsum) {
                //nBotClusters++;
                if (hasTopTrack && !hasBottomTrack)
                    aida.histogram2D("LowEsum TopTrack ECal Cluster Position").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                else if (!hasTopTrack && hasBottomTrack)
                    aida.histogram2D("LowEsum BottomTrack ECal Cluster Position").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
                //aida.histogram1D("Bottom ECal Cluster Energy").fill(cluster.getEnergy());
            }

            //            if (tracksPresent) {
            //                if (cluster.getPosition()[1] > 0) {
            //                    aida.histogram2D("Top ECal Cluster Position (>0 tracks)").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
            //                }
            //                if (cluster.getPosition()[1] < 0) {
            //                    aida.histogram2D("Bottom ECal Cluster Position (>0 tracks)").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
            //                }
            //
            //                if (cluster.getEnergy() > 0.1) {
            //                    if (cluster.getPosition()[1] > 0) {
            //                        aida.histogram2D("Top ECal Cluster Position (E>0.1,>0 tracks)").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
            //                        aida.histogram2D("Top ECal Cluster Position w_E (E>0.1,>0 tracks)").fill(cluster.getPosition()[0], cluster.getPosition()[1], cluster.getEnergy());
            //                    }
            //                    if (cluster.getPosition()[1] < 0) {
            //                        aida.histogram2D("Bottom ECal Cluster Position (E>0.1,>0 tracks)").fill(cluster.getPosition()[0], cluster.getPosition()[1]);
            //                        aida.histogram2D("Bottom ECal Cluster Position w_E (E>0.1,>0 tracks)").fill(cluster.getPosition()[0], cluster.getPosition()[1], cluster.getEnergy());
            //                    }
            //                }
            //            }

        }
        //
        //        aida.histogram1D("Number of Clusters Top").fill(nTopClusters);
        //        aida.histogram1D("Number of Clusters Bot").fill(nBotClusters);
    }

    private void doBasicTracks(List<Track> tracks) {
        int ntracksTop = 0;
        int ntracksBot = 0;
        double momentum_param = 2.99792458e-04;

        aida.histogram1D("Tracks per Event", 3, 0, 3).fill(tracks.size());

        for (Track trk : tracks) {

            boolean isTop = false;
            if (trk.getTrackerHits().get(0).getPosition()[2] > 0) {
                isTop = true;
            }

            double pt = Math.abs((1 / trk.getTrackStates().get(0).getOmega()) * bfield * momentum_param);
            double pz = pt * Math.cos(trk.getTrackStates().get(0).getPhi());
            double px = pt * Math.sin(trk.getTrackStates().get(0).getPhi());
            double py = pt * trk.getTrackStates().get(0).getTanLambda();
            aida.histogram1D("Track Momentum (Pz)").fill(pz);
            aida.histogram1D("Track Momentum (Py)").fill(py);
            aida.histogram1D("Track Momentum (Px)").fill(px);
            aida.histogram1D("Track Chi2").fill(trk.getChi2());

            aida.histogram1D("Hits per Track").fill(trk.getTrackerHits().size());
            if (isTop)
                aida.histogram1D("Hits per Track Top").fill(trk.getTrackerHits().size());
            else
                aida.histogram1D("Hits per Track Bottom").fill(trk.getTrackerHits().size());

            aida.histogram1D("d0 ").fill(trk.getTrackStates().get(0).getParameter(ParameterName.d0.ordinal()));
            aida.histogram1D("sinphi ").fill(Math.sin(trk.getTrackStates().get(0).getParameter(ParameterName.phi0.ordinal())));
            aida.histogram1D("omega ").fill(trk.getTrackStates().get(0).getParameter(ParameterName.omega.ordinal()));
            aida.histogram1D("tan(lambda) ").fill(trk.getTrackStates().get(0).getParameter(ParameterName.tanLambda.ordinal()));
            aida.histogram1D("z0 ").fill(trk.getTrackStates().get(0).getParameter(ParameterName.z0.ordinal()));

            if (isTop) {
                aida.histogram1D("Top Track Momentum (Px)").fill(px);
                aida.histogram1D("Top Track Momentum (Py)").fill(py);
                aida.histogram1D("Top Track Momentum (Pz)").fill(pz);
                aida.histogram1D("Top Track Chi2").fill(trk.getChi2());

                aida.histogram1D("d0 Top").fill(trk.getTrackStates().get(0).getParameter(ParameterName.d0.ordinal()));
                aida.histogram1D("sinphi Top").fill(Math.sin(trk.getTrackStates().get(0).getParameter(ParameterName.phi0.ordinal())));
                aida.histogram1D("omega Top").fill(trk.getTrackStates().get(0).getParameter(ParameterName.omega.ordinal()));
                aida.histogram1D("tan(lambda) Top").fill(trk.getTrackStates().get(0).getParameter(ParameterName.tanLambda.ordinal()));
                aida.histogram1D("z0 Top").fill(trk.getTrackStates().get(0).getParameter(ParameterName.z0.ordinal()));
                ntracksTop++;
            } else {
                aida.histogram1D("Bottom Track Momentum (Px)").fill(px);
                aida.histogram1D("Bottom Track Momentum (Py)").fill(py);
                aida.histogram1D("Bottom Track Momentum (Pz)").fill(pz);
                aida.histogram1D("Bottom Track Chi2").fill(trk.getChi2());

                aida.histogram1D("d0 Bottom").fill(trk.getTrackStates().get(0).getParameter(ParameterName.d0.ordinal()));
                aida.histogram1D("sinphi Bottom").fill(Math.sin(trk.getTrackStates().get(0).getParameter(ParameterName.phi0.ordinal())));
                aida.histogram1D("omega Bottom").fill(trk.getTrackStates().get(0).getParameter(ParameterName.omega.ordinal()));
                aida.histogram1D("tan(lambda) Bottom").fill(-1.0 * trk.getTrackStates().get(0).getParameter(ParameterName.tanLambda.ordinal()));
                aida.histogram1D("z0 Bottom").fill(trk.getTrackStates().get(0).getParameter(ParameterName.z0.ordinal()));
                ntracksBot++;
            }
        }

        aida.histogram1D("Tracks per Event Bot").fill(ntracksBot);
        aida.histogram1D("Tracks per Event Top").fill(ntracksTop);
    }

    private void doHitsOnTrack(Track trk) {
        Map<HpsSiSensor, Integer> stripHitsOnTrack = new HashMap<HpsSiSensor, Integer>();
        List<TrackerHit> hitsOnTrack = trk.getTrackerHits();

        for (TrackerHit hit : hitsOnTrack) {

            HpsSiSensor sensor = ((HpsSiSensor) ((RawTrackerHit) hit.getRawHits().get(0)).getDetectorElement());

            if (stripHitsOnTrack.containsKey(sensor)) {
                stripHitsOnTrack.put(sensor, stripHitsOnTrack.get(sensor) + 1);
            } else {
                stripHitsOnTrack.put(sensor, 1);
            }
        }

        for (Map.Entry<HpsSiSensor, Integer> sensor : stripHitsOnTrack.entrySet()) {
            aida.histogram1D(sensor.getKey().getName() + " strip hits on track").fill(stripHitsOnTrack.get(sensor.getKey()));
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

    private void doMissingHits(List<Track> GBLTrackList, List<TrackerHit> hthList, List<TrackerHit> sthList, List<Cluster> clusterList, int trigPhase) {
        int ntracksTop = 0;
        int ntracksBot = 0;
        int ntracksPos = 0;
        int ntracksEle = 0;

        double maxClusE = 0;
        Cluster topClus = null;
        Cluster botClus = null;
        for (Cluster clus : clusterList) {
            if (clus.getPosition()[1] > 0)
                topClus = clus;
            else
                botClus = clus;
            if (clus.getEnergy() > maxClusE)
                maxClusE = clus.getEnergy();
        }

        for (Track trk : GBLTrackList) {
            boolean isTop = trk.getTrackerHits().get(0).getPosition()[2] > 0;
            int charge = -(int) Math.signum(TrackUtils.getR(trk.getTrackStates().get(0)));
            if (isTop)
                ntracksTop++;
            else
                ntracksBot++;
            if (charge == -1)
                ntracksEle++;
            else if (charge == 1)
                ntracksPos++;
        }
        aida.histogram1D("Tracks per Event Bot").fill(ntracksBot);
        aida.histogram1D("Tracks per Event Top").fill(ntracksTop);

        Map<Integer, Integer> HitsInLayer = new HashMap<Integer, Integer>();
        //System.out.println("Starting 3D hits in event:");
        for (TrackerHit hth : hthList) {
            HpsSiSensor sensor = ((HpsSiSensor) ((RawTrackerHit) hth.getRawHits().get(0)).getDetectorElement());
            int lay = sensor.getLayerNumber() / 2 + 1;
            int n;
            if (HitsInLayer.containsKey(lay)) {
                n = HitsInLayer.get(lay);
            } else {
                n = 0;
            }
            //System.out.printf("lay %d pos %f %f %f \n", lay, hth.getPosition()[0], hth.getPosition()[1], hth.getPosition()[2]);
            if ((hth.getPosition()[1] > 0 && ntracksTop == 0) || (hth.getPosition()[1] < 0 && ntracksBot == 0)) {
                //System.out.printf("lay %d time %f \n", lay, hth.getTime());
                if (ntracksBot == 0) {
                    aida.histogram1D("Times of Bottom 3D-Hits in Events with No Bottom Track").fill(hth.getTime());
                }
                if (ntracksTop == 0) {
                    aida.histogram1D("Times of Top 3D-Hits in Events with No Top Track").fill(hth.getTime());
                }
                n++;
            }
            HitsInLayer.put(lay, n);
        }
        int numLayHit = 0;
        for (int lay = 1; lay <= 6; lay++) {
            if (HitsInLayer.containsKey(lay) && HitsInLayer.get(lay) > 0)
                numLayHit++;
        }
        if (ntracksTop == 0) {
            aida.histogram1D("Top Layers with 3D-Hits in Events with No Top Track").fill(numLayHit);
        }
        if (ntracksBot == 0) {
            aida.histogram1D("Bottom Layers with 3D-Hits in Events with No Bottom Track").fill(numLayHit);
        }

        //        System.out.println("Starting 2D hits in event:");
        HitsInLayer.clear();
        for (TrackerHit stripHit : sthList) {

            HpsSiSensor sensor = ((HpsSiSensor) ((RawTrackerHit) stripHit.getRawHits().get(0)).getDetectorElement());
            int lay = sensor.getLayerNumber();
            int n;
            if (HitsInLayer.containsKey(lay)) {
                n = HitsInLayer.get(lay);
            } else {
                n = 0;
            }
            //System.out.printf("lay %d pos %f %f %f \n", lay, stripHit.getPosition()[0], stripHit.getPosition()[1], stripHit.getPosition()[2]);
            if (lay >= 7) {
                if (stripHit.getPosition()[1] < 0 && ntracksBot == 0) {
                    aida.histogram1D("Times of Bottom 2D-Hits in Events with No Bottom Track").fill(stripHit.getTime());
                    String temp = String.format("Times of Bottom 2D-Hits in Events with No Bottom Track - Layer %d", lay);
                    aida.histogram1D(temp).fill(stripHit.getTime());
                    temp = String.format("Times of Bottom 2D-Hits in Events with No Bottom Track - Trigger Phase %d", trigPhase);
                    aida.histogram1D(temp).fill(stripHit.getTime());
                    if (ntracksEle > 0)
                        aida.histogram1D("Times of Bottom 2D-Hits in Electron Events with No Bottom Track").fill(stripHit.getTime());
                    else if (ntracksPos > 0)
                        aida.histogram1D("Times of Bottom 2D-Hits in Positron Events with No Bottom Track").fill(stripHit.getTime());
                    if (maxClusE > 0.3)
                        aida.histogram1D("Times of Bottom 2D-Hits in HighMaxClusE Events with No Bottom Track").fill(stripHit.getTime());
                    else
                        aida.histogram1D("Times of Bottom 2D-Hits in LowMaxClusE Events with No Bottom Track").fill(stripHit.getTime());
                    if (botClus != null && topClus != null) {
                        if ((botClus.getEnergy() < 0.3) && (topClus.getEnergy() > 0.3))
                            aida.histogram1D("Times of Bottom 2D-Hits in LowMaxClusEBot Events with No Bottom Track").fill(stripHit.getTime());
                        else if ((botClus.getEnergy() > 0.3) && (topClus.getEnergy() < 0.3))
                            aida.histogram1D("Times of Bottom 2D-Hits in LowMaxClusETop Events with No Bottom Track").fill(stripHit.getTime());
                    }
                }
                if (stripHit.getPosition()[1] > 0 && ntracksTop == 0) {
                    aida.histogram1D("Times of Top 2D-Hits in Events with No Top Track").fill(stripHit.getTime());
                    String temp = String.format("Times of Top 2D-Hits in Events with No Top Track - Layer %d", lay);
                    aida.histogram1D(temp).fill(stripHit.getTime());
                    temp = String.format("Times of Top 2D-Hits in Events with No Top Track - Trigger Phase %d", trigPhase);
                    aida.histogram1D(temp).fill(stripHit.getTime());
                    if (ntracksEle > 0)
                        aida.histogram1D("Times of Top 2D-Hits in Electron Events with No Top Track").fill(stripHit.getTime());
                    else if (ntracksPos > 0)
                        aida.histogram1D("Times of Top 2D-Hits in Positron Events with No Top Track").fill(stripHit.getTime());
                    if (maxClusE > 0.3)
                        aida.histogram1D("Times of Top 2D-Hits in HighMaxClusE Events with No Top Track").fill(stripHit.getTime());
                    else
                        aida.histogram1D("Times of Top 2D-Hits in LowMaxClusE Events with No Top Track").fill(stripHit.getTime());
                    if (botClus != null && topClus != null) {
                        if ((botClus.getEnergy() < 0.3) && (topClus.getEnergy() > 0.3))
                            aida.histogram1D("Times of Top 2D-Hits in LowMaxClusEBot Events with No Top Track").fill(stripHit.getTime());
                        if ((botClus.getEnergy() > 0.3) && (topClus.getEnergy() < 0.3))
                            aida.histogram1D("Times of Top 2D-Hits in LowMaxClusETop Events with No Top Track").fill(stripHit.getTime());
                    }
                }

                //System.out.printf("lay %d time %f \n", lay, stripHit.getTime());
            }
            if (Math.abs(stripHit.getTime()) <= 12.0)
                n++;
            HitsInLayer.put(lay, n);
        }
        numLayHit = 0;
        for (int lay = 1; lay <= 12; lay++) {
            if (HitsInLayer.containsKey(lay) && HitsInLayer.get(lay) > 0)
                numLayHit++;
        }
        if (ntracksTop == 0) {
            aida.histogram1D("Top Layers with 2D-Hits in Events with No Top Track").fill(numLayHit);
        }
        if (ntracksBot == 0) {
            aida.histogram1D("Bottom Layers with 2D-Hits in Events with No Bottom Track").fill(numLayHit);
        }
    }

    private boolean doRecoParticles(EventHeader event, List<ReconstructedParticle> fspList, List<Track> tracks, List<Cluster> clusterList, List<GenericObject> TriggerBank) {

        //boolean hasClusters = true;
        //        boolean passesEsum = false;
        boolean highEclusTop = false;
        boolean highEclusBot = false;
        boolean lowEclusBot = false;
        boolean lowEclusTop = false;
        List<SSPCluster> sspClusters = null;
        boolean isOK = false;
        boolean hasTop = false;
        boolean hasBot = false;
        boolean isMatched = false;
        Cluster topClus = null;
        Cluster botClus = null;
        double totE = 0;
        double trigE = 0;
        //
        if (clusterList != null) {
            for (Cluster clus : clusterList) {
                totE += clus.getEnergy();
                if (clus.getPosition()[1] > 0) {
                    if (clus.getEnergy() > 0.15)
                        highEclusTop = true;
                    else
                        lowEclusTop = true;
                    hasTop = true;
                    topClus = clus;
                } else {
                    if (clus.getEnergy() > 0.15)
                        highEclusBot = true;
                    else
                        lowEclusBot = true;
                    hasBot = true;
                    botClus = clus;
                }
            }
            isOK = (hasTop && hasBot);

        }

        if (!isOK)
            return false;

        if ((ClusterUtilities.getSeedHitTime(topClus) - ClusterUtilities.getSeedHitTime(botClus)) > 2.0)
            return false;

        //if (!(totE > 0.8 && totE < 1.2))
        //    return;

        //if (fspList.isEmpty())
        //    return;

        //        ReconstructedParticle mostEnergeticEle = null;
        //        ReconstructedParticle mostEnergeticPos = null;
        //        for (ReconstructedParticle fsp : fspList) {
        //            if (fsp.getCharge() < 0) {
        //                if (mostEnergeticEle == null)
        //                    mostEnergeticEle = fsp;
        //                else if (fsp.getEnergy() > mostEnergeticEle.getEnergy())
        //                    mostEnergeticEle = fsp;
        //            } else if (fsp.getCharge() > 0) {
        //                if (mostEnergeticPos == null)
        //                    mostEnergeticPos = fsp;
        //                else if (fsp.getEnergy() > mostEnergeticPos.getEnergy())
        //                    mostEnergeticPos = fsp;
        //            }
        //        }
        //        if (mostEnergeticPos != null && mostEnergeticEle != null) {
        //            Hep3Vector temp = VecOp.add(mostEnergeticPos.getMomentum(), mostEnergeticEle.getMomentum());
        //            if (temp.magnitude() > 0.8)
        //                passesEsum = true;
        //        }

        //        double fsppt = 0;
        //        double fsppz = 0;
        //        if (mostEnergeticFsp != null) {
        //            fsppt = Math.abs((1.0 / mostEnergeticFsp.getTracks().get(0).getTrackStates().get(0).getOmega()) * bfield * 2.99792458e-04);
        //            fsppz = fsppt * Math.cos(mostEnergeticFsp.getTracks().get(0).getTrackStates().get(0).getPhi());
        //            fsppz *= mostEnergeticFsp.getCharge();
        //        }

        //        for (ReconstructedParticle fsp : fspList) {
        //            if (fsp.getTracks().isEmpty())
        //                continue;
        //
        //            Track trk = fsp.getTracks().get(0);
        //            boolean isEle = false;
        //            boolean isPos = false;
        //            boolean isTop = false;
        //            if (fsp.getCharge() < 0)
        //                isEle = true;
        //            else if (fsp.getCharge() > 0)
        //                isPos = true;
        //            if (trk.getTrackerHits().get(0).getPosition()[2] > 0)
        //                isTop = true;
        //
        //            double pt = Math.abs((1.0 / trk.getTrackStates().get(0).getOmega()) * bfield * 2.99792458e-04);
        //            double pz = pt * Math.cos(trk.getTrackStates().get(0).getPhi());
        //
        //            if (isTop) {
        //                aida.histogram1D("Reco Particles Top Track Chi2").fill(trk.getChi2());
        //                aida.histogram1D("Reco Particles Top Track Pz").fill(pz);
        //                if (isEle)
        //                    aida.histogram1D("Reco Electrons Top Track Chi2").fill(trk.getChi2());
        //                else if (isPos)
        //                    aida.histogram1D("Reco Positrons Top Track Chi2").fill(trk.getChi2());
        //            } else {
        //                aida.histogram1D("Reco Particles Bottom Track Chi2").fill(trk.getChi2());
        //                aida.histogram1D("Reco Particles Bottom Track Pz").fill(pz);
        //                if (isEle)
        //                    aida.histogram1D("Reco Electrons Bottom Track Chi2").fill(trk.getChi2());
        //                else if (isPos)
        //                    aida.histogram1D("Reco Positrons Bottom Track Chi2").fill(trk.getChi2());
        //            }
        //
        //            if (!hasClusters)
        //                continue;
        //        }

        //            boolean isSingles = false;
        //            boolean isSingles0 = false;
        //            boolean isSingles1 = false;
        //            boolean isPairs = false;
        for (GenericObject gob : TriggerBank) {
            //                if (AbstractIntData.getTag(gob) == TIData.BANK_TAG) {
            //
            //                    TIData tid = new TIData(gob);
            //                    if (tid.isSingle0Trigger()) {
            //                        isSingles = true;
            //                        isSingles0 = true;
            //
            //                    }
            //                    if (tid.isSingle1Trigger()) {
            //                        isSingles = true;
            //                        isSingles1 = true;
            //
            //                    }
            //                    if (tid.isPair0Trigger() || tid.isPair1Trigger()) {
            //                        isPairs = true;
            //
            //                    }
            //                }
            if (AbstractIntData.getTag(gob) == SSPData.BANK_TAG) {
                SSPData sspBank = new SSPData(gob);
                sspClusters = sspBank.getClusters();
            }
        }

        //            if (sspClusters != null) {
        //                double totE = 0;
        //                for (SSPCluster clus : sspClusters) {
        //                    totE += clus.getEnergy();
        //                    if (clus.getYIndex() < 0) {
        //                        if (clus.getEnergy() < 0.1)
        //                            lowEclusBot = true;
        //                        else if (clus.getEnergy() > 0.15)
        //                            highEclusBot = true;
        //                    }
        //                                        if (clus.getEnergy() > 0.5) {
        //                                            if (clus.getYIndex() > 0)
        //                                                highEclusTop = true;
        //                                            else
        //                                                highEclusBot = true;
        //                                        }
        //                }
        //                if (totE > 0.8)
        //                    passesEsum = true;
        //            }

        //            List<Cluster> matchedClusters = fsp.getClusters();

        if (!clusterList.isEmpty() && sspClusters != null && !sspClusters.isEmpty()) {
            isMatched = true;
            for (SSPCluster sspclus : sspClusters) {
                boolean m = false;
                for (Cluster matchedCluster : clusterList) {
                    if (isMatchedCluster(matchedCluster, sspclus)) {
                        m = true;
                        break;
                    }
                }
                if (m == false) {
                    isMatched = false;
                    //break;
                }
                trigE += sspclus.getEnergy();
            }
            //                    if (matchedCluster.getEnergy() > 0.06 && isSingles0) {
            //                        isMatched = true;
            //                        break;
            //                    }
            //                    if (matchedCluster.getEnergy() > 0.4 && isSingles1) {
            //                        isMatched = true;
            //                        break;
            //                    }
            //                    if (matchedCluster.getEnergy() > 0.054 && isPairs) {
            //                        isMatched = true;
            //                        break;
            //                    }

        }

        if (!isMatched)
            return false;

        int nTracksTop = 0;
        int nTracksBot = 0;
        for (Track trk : tracks) {
            boolean isTop = trk.getTrackerHits().get(0).getPosition()[2] > 0;
            int charge = -(int) Math.signum(TrackUtils.getR(trk.getTrackStates().get(0)));
            if (isTop)
                nTracksTop++;
            else
                nTracksBot++;
            double pt = Math.abs((1 / trk.getTrackStates().get(0).getOmega()) * bfield * 2.99792458e-04);
            double pz = pt * Math.cos(trk.getTrackStates().get(0).getPhi());

            //            if (isTop) {
            //                if (isEle) {
            //                    aida.histogram1D("Reco Pairs Electrons Top Track Chi2").fill(trk.getChi2());
            //                    aida.histogram1D("Reco Pairs Electrons Top Track Pz").fill(pz);
            //                } else if (isPos) {
            //                    aida.histogram1D("Reco Pairs Positrons Top Track Chi2").fill(trk.getChi2());
            //                    aida.histogram1D("Reco Pairs Positrons Top Track Pz").fill(pz);
            //                }

            //                    if (isPos) {
            //                        aida.histogram1D("Reco Pairs Positrons Top Track Chi2").fill(trk.getChi2());
            //                    } else if (isEle) {
            //                        aida.histogram1D("Reco Pairs Electrons Top Track Chi2").fill(trk.getChi2());
            //                    }
            //if (isMatched)
            //    aida.histogram1D("Reco MatchedPairs Particles Top Track Pz").fill(pz);
            //                    if (highEclusTop) {
            //                        aida.histogram1D("Reco Particles Top Track Pz HighEclusTopEvent").fill(pz);
            //                        if (tracks.size() == 1)
            //                            aida.histogram1D("Reco Particles Top Track Pz HighEclusTopEvent Tracks=1").fill(pz);
            //                        else if (tracks.size() == 2)
            //                            aida.histogram1D("Reco Particles Top Track Pz HighEclusTopEvent Tracks=2").fill(pz);
            //                        else if (tracks.size() >= 3)
            //                            aida.histogram1D("Reco Particles Top Track Pz HighEclusTopEvent minTracks=3").fill(pz);
            //                    }
            //                    if (lowEclusBot && !highEclusBot) {
            //                        aida.histogram1D("Reco Particles Top Track Pz LowEclusBotEvent").fill(pz);
            //                        if (tracks.size() == 1)
            //                            aida.histogram1D("Reco Particles Top Track Pz LowEclusBotEvent Tracks=1").fill(pz);
            //                        else if (tracks.size() == 2)
            //                            aida.histogram1D("Reco Particles Top Track Pz LowEclusBotEvent Tracks=2").fill(pz);
            //                        else if (tracks.size() >= 3)
            //                            aida.histogram1D("Reco Particles Top Track Pz LowEclusBotEvent minTracks=3").fill(pz);
            //                    }
            //                    if (highEclusBot) {
            //                        aida.histogram1D("Reco Particles Top Track Pz HighEclusBotEvent").fill(pz);
            //                        if (tracks.size() == 1)
            //                            aida.histogram1D("Reco Particles Top Track Pz HighEclusBotEvent Tracks=1").fill(pz);
            //                        else if (tracks.size() == 2)
            //                            aida.histogram1D("Reco Particles Top Track Pz HighEclusBotEvent Tracks=2").fill(pz);
            //                        else if (tracks.size() >= 3)
            //                            aida.histogram1D("Reco Particles Top Track Pz HighEclusBotEvent minTracks=3").fill(pz);
            //                    }
            //if (passesEsum)
            //  aida.histogram1D("Reco Particles Top Track Pz HighEsumEvent").fill(pz);
            //            } else {
            //                if (isEle) {
            //                    aida.histogram1D("Reco Pairs Electrons Bottom Track Chi2").fill(trk.getChi2());
            //                    aida.histogram1D("Reco Pairs Electrons Bottom Track Pz").fill(pz);
            //                } else if (isPos) {
            //                    aida.histogram1D("Reco Pairs Positrons Bottom Track Chi2").fill(trk.getChi2());
            //                    aida.histogram1D("Reco Pairs Positrons Bottom Track Pz").fill(pz);
            //                }
            //            }

            if (highEclusBot) {
                if (isTop) {
                    aida.histogram1D("Reco Top Track Pz HighEclusBotEvent").fill(pz);
                    aida.histogram1D("Reco Top Track Chi2 HighEclusBotEvent").fill(trk.getChi2());
                } else {
                    aida.histogram1D("Reco Bottom Track Pz HighEclusBotEvent").fill(pz);
                    aida.histogram1D("Reco Bottom Track Chi2 HighEclusBotEvent").fill(trk.getChi2());
                }
            }
            if (highEclusTop) {
                if (isTop) {
                    aida.histogram1D("Reco Top Track Pz HighEclusTopEvent").fill(pz);
                    aida.histogram1D("Reco Top Track Chi2 HighEclusTopEvent").fill(trk.getChi2());
                } else {
                    aida.histogram1D("Reco Bottom Track Chi2 HighEclusTopEvent").fill(trk.getChi2());
                    aida.histogram1D("Reco Bottom Track Pz HighEclusTopEvent").fill(pz);
                }
            }
            if (lowEclusTop) {
                if (isTop) {
                    aida.histogram1D("Reco Top Track Pz LowEclusTopEvent").fill(pz);
                    aida.histogram1D("Reco Top Track Chi2 LowEclusTopEvent").fill(trk.getChi2());
                } else {
                    aida.histogram1D("Reco Bottom Track Chi2 LowEclusTopEvent").fill(trk.getChi2());
                    aida.histogram1D("Reco Bottom Track Pz LowEclusTopEvent").fill(pz);
                }
            }
            if (lowEclusBot) {
                if (isTop) {
                    aida.histogram1D("Reco Top Track Pz LowEclusBotEvent").fill(pz);
                    aida.histogram1D("Reco Top Track Chi2 LowEclusBotEvent").fill(trk.getChi2());
                } else {
                    aida.histogram1D("Reco Bottom Track Chi2 LowEclusBotEvent").fill(trk.getChi2());
                    aida.histogram1D("Reco Bottom Track Pz LowEclusBotEvent").fill(pz);
                }
            }

            if (isTop) {
                if (charge == -1)
                    aida.histogram1D("Reco Pairs Electrons Top Track Pz").fill(pz);
                else if (charge == 1)
                    aida.histogram1D("Reco Pairs Positrons Top Track Pz").fill(pz);
            } else {
                if (charge == -1)
                    aida.histogram1D("Reco Pairs Electrons Bottom Track Pz").fill(pz);
                else if (charge == 1)
                    aida.histogram1D("Reco Pairs Positrons Bottom Track Pz").fill(pz);
            }

        }

        for (Track trk : tracks) {
            TrackState tsAtEcal = TrackStateUtils.getTrackStateAtECal(trk);
            Hep3Vector atEcal = null;
            if (tsAtEcal == null)
                continue;
            atEcal = new BasicHep3Vector(tsAtEcal.getReferencePoint());
            atEcal = CoordinateTransformations.transformVectorToDetector(atEcal);
            if (totE > 0.7) {
                if (nTracksTop > 0 && nTracksBot == 0)
                    aida.histogram2D("HighEsum TopTrack ECal Extrap Position").fill(atEcal.x(), atEcal.y());
                else if (nTracksTop == 0 && nTracksBot > 0)
                    aida.histogram2D("HighEsum BottomTrack ECal Extrap Position").fill(atEcal.x(), atEcal.y());
            } else if (totE < 0.55) {
                double pt = Math.abs((1 / trk.getTrackStates().get(0).getOmega()) * bfield * 2.99792458e-04);
                double pz = pt * Math.cos(trk.getTrackStates().get(0).getPhi());
                if (nTracksTop > 0 && nTracksBot == 0) {
                    aida.histogram2D("LowEsum TopTrack ECal Extrap Position").fill(atEcal.x(), atEcal.y());
                    aida.histogram2D("LowEsum TopTrack ClusterE vs TrackPz").fill(pz, topClus.getEnergy());
                } else if (nTracksTop == 0 && nTracksBot > 0) {
                    aida.histogram2D("LowEsum BottomTrack ECal Extrap Position").fill(atEcal.x(), atEcal.y());
                    aida.histogram2D("LowEsum BottomTrack ClusterE vs TrackPz").fill(pz, botClus.getEnergy());
                }
            }
        }

        if (nTracksTop > 0) {
            if (nTracksBot > 0) {
                aida.histogram2D("Top vs Bottom ClusE with tracks in both").fill(topClus.getEnergy(), botClus.getEnergy());
                aida.histogram1D("Reco ESum with tracks in both").fill(totE);
                aida.histogram1D("Trigger ESum with tracks in both").fill(trigE);
            } else {
                aida.histogram2D("Top vs Bottom ClusE with track in top").fill(topClus.getEnergy(), botClus.getEnergy());
                aida.histogram1D("Reco ESum with track in top").fill(totE);
                aida.histogram1D("Trigger ESum with track in top").fill(trigE);
            }
        } else {
            if (nTracksBot > 0) {
                aida.histogram2D("Top vs Bottom ClusE with track in bottom").fill(topClus.getEnergy(), botClus.getEnergy());
                aida.histogram1D("Reco ESum with track in bottom").fill(totE);
                aida.histogram1D("Trigger ESum with track in bottom").fill(trigE);
            } else {
                aida.histogram2D("Top vs Bottom ClusE with no tracks").fill(topClus.getEnergy(), botClus.getEnergy());
                aida.histogram1D("Reco ESum with no tracks").fill(totE);
                aida.histogram1D("Trigger ESum with no tracks").fill(trigE);
            }
        }

        doECalClusters(clusterList, totE > 0.7, totE < 0.55, nTracksTop > 0, nTracksBot > 0);

        //                    if (isPos) {
        //                        aida.histogram1D("Reco Pairs Positrons Bottom Track Chi2").fill(trk.getChi2());
        //                    } else if (isEle) {
        //                        aida.histogram1D("Reco Pairs Electrons Bottom Track Chi2").fill(trk.getChi2());
        //                    }
        //                    if (isMatched)
        //                        aida.histogram1D("Reco MatchedPairs Particles Bottom Track Pz").fill(pz);
        //                    if (highEclusTop) {
        //                        aida.histogram1D("Reco Particles Bottom Track Pz HighEclusTopEvent").fill(pz);
        //                        if (tracks.size() == 1)
        //                            aida.histogram1D("Reco Particles Bottom Track Pz HighEclusTopEvent Tracks=1").fill(pz);
        //                        else if (tracks.size() == 2)
        //                            aida.histogram1D("Reco Particles Bottom Track Pz HighEclusTopEvent Tracks=2").fill(pz);
        //                        else if (tracks.size() >= 3)
        //                            aida.histogram1D("Reco Particles Bottom Track Pz HighEclusTopEvent minTracks=3").fill(pz);
        //                    }
        //                    if (lowEclusBot && !highEclusBot) {
        //                        aida.histogram1D("Reco Particles Bottom Track Pz LowEclusBotEvent").fill(pz);
        //                        if (tracks.size() == 1)
        //                            aida.histogram1D("Reco Particles Bottom Track Pz LowEclusBotEvent Tracks=1").fill(pz);
        //                        else if (tracks.size() == 2)
        //                            aida.histogram1D("Reco Particles Bottom Track Pz LowEclusBotEvent Tracks=2").fill(pz);
        //                        else if (tracks.size() >= 3)
        //                            aida.histogram1D("Reco Particles Bottom Track Pz LowEclusBotEvent minTracks=3").fill(pz);
        //                    }
        //                    if (highEclusBot) {
        //                        aida.histogram1D("Reco Particles Bottom Track Pz HighEclusBotEvent").fill(pz);
        //                        if (tracks.size() == 1)
        //                            aida.histogram1D("Reco Particles Bottom Track Pz HighEclusBotEvent Tracks=1").fill(pz);
        //                        else if (tracks.size() == 2)
        //                            aida.histogram1D("Reco Particles Bottom Track Pz HighEclusBotEvent Tracks=2").fill(pz);
        //                        else if (tracks.size() >= 3)
        //                            aida.histogram1D("Reco Particles Bottom Track Pz HighEclusBotEvent minTracks=3").fill(pz);
        //                    }
        //                    if (passesEsum)
        //                        aida.histogram1D("Reco Particles Bottom Track Pz HighEsumEvent").fill(pz);
        //                }

        //                if (pz < 0.8) {
        //                    if (isTop) {
        //                        aida.histogram1D("Reco PairsLowPz Particles Top Track Chi2").fill(trk.getChi2());
        //                        if (isPos) {
        //                            aida.histogram1D("Reco PairsLowPz Positrons Top Track Chi2").fill(trk.getChi2());
        //                        } else if (isEle) {
        //                            aida.histogram1D("Reco PairsLowPz Electrons Top Track Chi2").fill(trk.getChi2());
        //                        }
        //                    } else {
        //                        aida.histogram1D("Reco PairsLowPz Particles Bottom Track Chi2").fill(trk.getChi2());
        //                        if (isPos) {
        //                            aida.histogram1D("Reco PairsLowPz Positrons Bottom Track Chi2").fill(trk.getChi2());
        //                        } else if (isEle) {
        //                            aida.histogram1D("Reco PairsLowPz Electrons Bottom Track Chi2").fill(trk.getChi2());
        //                        }
        //                    }
        //
        //                    if (isMatched) {
        //                        if (isTop) {
        //
        //                            aida.histogram1D("Reco PairsLowPzMatched Particles Top Track Chi2").fill(trk.getChi2());
        //                            aida.histogram1D("Reco PairsLowPzMatched Particles Top Track d0").fill(trk.getTrackStates().get(0).getD0());
        //                            aida.histogram1D("Reco PairsLowPzMatched Particles Top Track z0").fill(trk.getTrackStates().get(0).getZ0());
        //                            if (isPos) {
        //                                aida.histogram1D("Reco PairsLowPzMatched Positrons Top Track Chi2").fill(trk.getChi2());
        //                            } else if (isEle) {
        //                                aida.histogram1D("Reco PairsLowPzMatched Electrons Top Track Chi2").fill(trk.getChi2());
        //                            }
        //                        } else {
        //                            aida.histogram1D("Reco PairsLowPzMatched Particles Bottom Track Chi2").fill(trk.getChi2());
        //                            aida.histogram1D("Reco PairsLowPzMatched Particles Bottom Track d0").fill(trk.getTrackStates().get(0).getD0());
        //                            aida.histogram1D("Reco PairsLowPzMatched Particles Bottom Track z0").fill(trk.getTrackStates().get(0).getZ0());
        //                            if (isPos) {
        //                                aida.histogram1D("Reco PairsLowPzMatched Positrons Bottom Track Chi2").fill(trk.getChi2());
        //                            } else if (isEle) {
        //                                aida.histogram1D("Reco PairsLowPzMatched Electrons Bottom Track Chi2").fill(trk.getChi2());
        //                            }
        //                        }
        //                    }
        //
        //                }
        //            }

        //            if (isOK) {
        //                if (isTop) {
        //                    aida.histogram1D("Reco Energetic Particles Top Track Chi2").fill(trk.getChi2());
        //                    aida.histogram1D("Reco Energetic Particles Top Track Pz").fill(pz);
        //                } else {
        //                    aida.histogram1D("Reco Energetic Particles Bottom Track Chi2").fill(trk.getChi2());
        //                    aida.histogram1D("Reco Energetic Particles Bottom Track Pz").fill(pz);
        //                }

        //                if (isSingles) {
        //                    if (isTop) {
        //                        aida.histogram1D("Reco FEE Particles Top Track Chi2").fill(trk.getChi2());
        //                        aida.histogram1D("Reco FEE Particles Top Track Pz").fill(pz);
        //                    } else {
        //                        aida.histogram1D("Reco FEE Particles Bottom Track Chi2").fill(trk.getChi2());
        //                        aida.histogram1D("Reco FEE Particles Bottom Track Pz").fill(pz);
        //                    }

        //                    if (isMatched) {
        //                                                if (isTop) {
        //                                                    aida.histogram1D("Reco FEEMatched Particles Top Track Chi2").fill(trk.getChi2());
        //                                                    aida.histogram1D("Reco FEEMatched Particles Top Track Pz").fill(pz);
        //                                                } else {
        //                                                    aida.histogram1D("Reco FEEMatched Particles Bottom Track Chi2").fill(trk.getChi2());
        //                                                    aida.histogram1D("Reco FEEMatched Particles Bottom Track Pz").fill(pz);
        //                                                }
        //                    }

        //                    if (pz < 0.8) {
        //                        if (mostEnergeticFsp == null)
        //                            mostEnergeticFsp = fsp;
        //                        // 4 plots: top/top, bottom/bottom, bottom/top, top/bottom
        //                        // 4-quadrant 2d plot, Pz (negative=electron, positive=positron)
        //                        double newPz = pz * fsp.getCharge();
        //                        if (isTop) {
        //                            aida.histogram1D("Reco FEELowPz Particles Top Track Chi2").fill(trk.getChi2());
        //                            aida.histogram1D("Reco FEELowPz Particles Top Track d0").fill(trk.getTrackStates().get(0).getD0());
        //                            aida.histogram1D("Reco FEELowPz Particles Top Track z0").fill(trk.getTrackStates().get(0).getZ0());
        //                            if (isPos) {
        //                                aida.histogram1D("Reco FEELowPz Positrons Top Track Chi2").fill(trk.getChi2());
        //                            } else if (isEle) {
        //                                aida.histogram1D("Reco FEELowPz Electrons Top Track Chi2").fill(trk.getChi2());
        //                            }
        //
        //                            if (mostEnergeticFsp.getTracks().get(0).getTrackerHits().get(0).getPosition()[2] > 0) {
        //                                aida.histogram2D("Reco FEELowPz vs MostEnergetic top-top").fill(newPz, fsppz);
        //                            } else {
        //                                aida.histogram2D("Reco FEELowPz vs MostEnergetic top-bottom").fill(newPz, fsppz);
        //                            }
        //                        } else {
        //                            aida.histogram1D("Reco FEELowPz Particles Bottom Track Chi2").fill(trk.getChi2());
        //                            aida.histogram1D("Reco FEELowPz Particles Bottom Track d0").fill(trk.getTrackStates().get(0).getD0());
        //                            aida.histogram1D("Reco FEELowPz Particles Bottom Track z0").fill(trk.getTrackStates().get(0).getZ0());
        //                            if (isPos) {
        //                                aida.histogram1D("Reco FEELowPz Positrons Bottom Track Chi2").fill(trk.getChi2());
        //                            } else if (isEle) {
        //                                aida.histogram1D("Reco FEELowPz Electrons Bottom Track Chi2").fill(trk.getChi2());
        //                            }
        //                            if (mostEnergeticFsp.getTracks().get(0).getTrackerHits().get(0).getPosition()[2] > 0) {
        //                                aida.histogram2D("Reco FEELowPz vs MostEnergetic bottom-top").fill(newPz, fsppz);
        //                            } else {
        //                                aida.histogram2D("Reco FEELowPz vs MostEnergetic bottom-bottom").fill(newPz, fsppz);
        //                            }
        //                        }
        //
        //                    }
        //                }
        //}
        return true;
    }

    private boolean isMatchedCluster(Cluster matchedCluster, SSPCluster sspclus) {
        if (Math.abs(matchedCluster.getEnergy() - sspclus.getEnergy()) > 0.05)
            return false;

        if (matchedCluster.getPosition()[1] > 0 && sspclus.getYIndex() < 0)
            return false;

        if (matchedCluster.getPosition()[1] < 0 && sspclus.getYIndex() > 0)
            return false;

        if (matchedCluster.getPosition()[0] < 0 && sspclus.getXIndex() > 0)
            return false;

        if (matchedCluster.getPosition()[0] > 0 && sspclus.getXIndex() < 0)
            return false;

        return true;
    }

    private void doOccupancy(List<RawTrackerHit> rawHits, List<GenericObject> TriggerBank) {

        boolean isSingles1 = false;
        boolean isPairs1 = false;
        boolean isPulser = false;

        for (GenericObject gob : TriggerBank) {
            if (AbstractIntData.getTag(gob) == TIData.BANK_TAG) {
                TIData tid = new TIData(gob);
                if (tid.isSingle1Trigger())
                    isSingles1 = true;
                if (tid.isPair1Trigger())
                    isPairs1 = true;
                if (tid.isPulserTrigger())
                    isPulser = true;
            }
        }

        for (RawTrackerHit rawHit : rawHits) {

            // Obtain the raw ADC samples for each of the six samples readout
            short[] adcValues = rawHit.getADCValues();

            // Find the sample that has the largest amplitude. This should
            // correspond to the peak of the shaper signal if the SVT is timed
            // in correctly. Otherwise, the maximum sample value will default
            // to 0.
            int maxAmplitude = 0;
            int maxSamplePositionFound = -1;
            for (int sampleN = 0; sampleN < 6; sampleN++) {
                if (adcValues[sampleN] > maxAmplitude) {
                    maxAmplitude = adcValues[sampleN];
                    maxSamplePositionFound = sampleN;
                }
            }

            String tempName = ((HpsSiSensor) (rawHit.getDetectorElement())).getName() + " - Occupancy";
            aida.histogram1D(tempName).fill(rawHit.getIdentifierFieldValue("strip"));
            String tempName2 = ((HpsSiSensor) (rawHit.getDetectorElement())).getName() + " - Max Sample Number";
            aida.histogram1D(tempName2).fill(maxSamplePositionFound);
            if (isPairs1) {
                aida.histogram1D(tempName + " Pairs1").fill(rawHit.getIdentifierFieldValue("strip"));
                aida.histogram1D(tempName2 + " Pairs1").fill(maxSamplePositionFound);
            }
            if (isSingles1) {
                aida.histogram1D(tempName + " Singles1").fill(rawHit.getIdentifierFieldValue("strip"));
                aida.histogram1D(tempName2 + " Singles1").fill(maxSamplePositionFound);
            }
            if (isPulser) {
                aida.histogram1D(tempName + " Pulser").fill(rawHit.getIdentifierFieldValue("strip"));
                aida.histogram1D(tempName2 + " Pulser").fill(maxSamplePositionFound);
            }

        }
    }

    private void doClustersOnTrack(Track trk, List<Cluster> clusters) {
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

    private void doComparison(List<Track> tracks, List<TrackerHit> hthList, List<Track> extraTracks, List<TrackerHit> extraHits) {
        //        System.out.println("Printing first track set");
        //        for (Track trk : tracks) {
        //            System.out.printf("%s \n", trk.toString());
        //        }
        //        System.out.println("Printing second track set");
        //        for (Track trk : extraTracks) {
        //            System.out.printf("%s \n", trk.toString());
        //        }

        //doBasicTracks(extraTracks);
        aida.histogram2D("HelicalTrackHits Per Event New vs Old").fill(hthList.size(), extraHits.size());

        List<List<Track>> temp = new ArrayList<List<Track>>();
        ArrayList<Track> temp1 = new ArrayList<Track>();
        temp1.addAll(extraTracks);
        temp.add(temp1);

        // new tracks shared amongst themselves
        DualAmbiguityResolver dar = new DualAmbiguityResolver(temp, AmbiMode.SHARED, 4, 0);
        dar.resolve();
        List<Track> sharedNew = new ArrayList<Track>();
        sharedNew.addAll(dar.getSharedTracks());
        dar.addToTrackList(tracks);
        dar.setMode(AmbiMode.DUPS);
        dar.resolve();
        List<Track> dupsNew = dar.getDuplicateTracks();
        aida.histogram1D("New Duplicate Tracks Per Event").fill(dupsNew.size());

        // with duplicates removed... find partials
        temp.get(0).removeAll(dupsNew);
        List<Track> temp2 = new ArrayList<Track>();
        temp2.addAll(tracks);
        temp2.removeAll(dupsNew);
        dar.resetResolver();
        dar.setMode(AmbiMode.PARTIALS);
        dar.initializeFromCollection(temp);
        List<Track> partialsNew = new ArrayList<Track>();
        partialsNew.addAll(dar.getPartialTracks());
        dar.addToTrackList(temp2);
        dar.resolve();
        List<Track> partialsDual = new ArrayList<Track>();
        partialsDual.addAll(dar.getPartialTracks());

        // with shared-new removed... find shared-old
        temp.get(0).removeAll(sharedNew);
        dar.resetResolver();
        dar.setMode(AmbiMode.SHARED);
        dar.initializeFromCollection(temp);
        dar.addToTrackList(temp2);
        dar.resolve();
        List<Track> sharedDual = new ArrayList<Track>();
        sharedDual.addAll(dar.getSharedTracks());

        // System.out.printf("shared with new %d old %d \n", sharedDual.size(), sharedNew.size());
        aida.histogram2D("New Shared Tracks with New vs Old - Per Event").fill(sharedDual.size(), sharedNew.size());
        aida.histogram2D("New Partial Tracks with New vs Old - Per Event").fill(partialsDual.size(), partialsNew.size());

        //temp.add(tracks);
        //        dar.resetResolver();
        //        dar.setMode(AmbiMode.DUPS);
        //        dar.initializeFromCollection(temp);
        //        dar.addToTrackList(tracks);
        //        dar.resolve();

        extraTracks.removeAll(dupsNew);
        doBasicTracks(extraTracks);
    }

    @Override
    public void process(EventHeader event) {
        aida.tree().cd("/");
        if (!event.hasCollection(TrackerHit.class, helicalTrackHitCollectionName)) {
            System.out.println(helicalTrackHitCollectionName + " does not exist; skipping event");
            return;
        }
        List<TrackerHit> hthList = event.get(TrackerHit.class, helicalTrackHitCollectionName);

        if (!event.hasCollection(Track.class, trackCollectionName)) {
            System.out.println(trackCollectionName + " does not exist; skipping event");
            aida.histogram1D("Number Tracks/Event").fill(0);
            return;
        }
        List<Track> tracks = event.get(Track.class, trackCollectionName);

        List<GenericObject> tb = event.get(GenericObject.class, "TriggerBank");

        List<Cluster> clusters = null;
        if (event.hasCollection(Cluster.class, ecalCollectionName)) {
            clusters = event.get(Cluster.class, ecalCollectionName);
        } else {
            doECalClusterPlots = false;
            doMatchedClusterPlots = false;
            doElectronPositronPlots = false;
        }

        List<ReconstructedParticle> fspList = null;
        if (event.hasCollection(ReconstructedParticle.class, "FinalStateParticles")) {
            fspList = event.get(ReconstructedParticle.class, "FinalStateParticles");
        } else
            doReconParticlePlots = false;

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

        List<RawTrackerHit> rawHits = null;
        if (!event.hasCollection(RawTrackerHit.class, "SVTRawTrackerHits")) {
            doOccupancyPlots = false;
        }
        // Get RawTrackerHit collection from event.
        rawHits = event.get(RawTrackerHit.class, "SVTRawTrackerHits");

        //        doBasicTracks(tracks);

        if (doOccupancyPlots)
            doOccupancy(rawHits, tb);

        // if (doECalClusterPlots)
        //     doECalClusters(clusters, tracks.size() > 0);

        if (doElectronPositronPlots) {
            eCanditates = new HashMap<Track, Cluster>();
            pCanditates = new HashMap<Track, Cluster>();
        }

        String extraHitCollection = "RotatedHelicalTrackHits-20s";
        List<TrackerHit> extraHits = null;
        if (event.hasCollection(TrackerHit.class, extraHitCollection)) {
            extraHits = event.get(TrackerHit.class, extraHitCollection);
        } else {
            doComparisonPlots = false;
        }
        String extraTracksCollection = "GBLTracks-20s";
        List<Track> extraTracks = null;
        if (event.hasCollection(Track.class, extraTracksCollection)) {
            extraTracks = event.get(Track.class, extraTracksCollection);
            //            System.out.println("Printing first track set");
            //            for (Track trk : tracks) {
            //                System.out.printf("%s \n", trk.toString());
            //            }
            //            //          
            //            System.out.printf("tracks size %d extratracks size %d \n", tracks.size(), extraTracks.size());
            //            System.out.println("Printing second track set");
            //            for (Track trk : extraTracks) {
            //                System.out.printf("%s \n", trk.toString());
            //            }
        } else {
            doComparisonPlots = false;
        }
        if (doComparisonPlots) {
            // partials test
            //if (!extraTracks.isEmpty())
            //    extraTracks.get(0).getTrackerHits().remove(0);

            doComparison(tracks, hthList, extraTracks, extraHits);
        }

        for (Track trk : tracks) {
            //            List<TrackerHit> temp = trk.getTrackerHits();
            //            for (TrackerHit hit : temp) {
            //                System.out.printf("%f \n", hit.getPosition()[0]);
            //            }

            if (doStripHitPlots)
                doStripHits(stripClusters, trk, trackDataTable);

            if (doHitsOnTrackPlots)
                doHitsOnTrack(trk);

            if (doResidualPlots)
                doResiduals(fittedHits, trk, trackResidualsTable);

            if (doAmplitudePlots)
                doAmplitude(fittedHits, trk);

            if (doMatchedClusterPlots)
                doClustersOnTrack(trk, clusters);
        }

        if (doElectronPositronPlots)
            doElectronPositron();

        if (doReconParticlePlots)
            if (doRecoParticles(event, fspList, tracks, clusters, tb))
                doMissingHits(tracks, hthList, stripClusters, clusters, (int) ((event.getTimeStamp() / 4) % 6));

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

    private void setupPlots() {

        // Basic tracks
        IHistogram1D trkPx = aida.histogram1D("Track Momentum (Px)", 100, -0.15, 0.15);
        IHistogram1D trkPy = aida.histogram1D("Track Momentum (Py)", 100, -0.15, 0.15);
        IHistogram1D trkPz = aida.histogram1D("Track Momentum (Pz)", 100, 0, 1.5);
        IHistogram1D trkChi2 = aida.histogram1D("Track Chi2", 100, 0, 100.0);

        IHistogram1D toptrkPx = aida.histogram1D("Top Track Momentum (Px)", 100, -0.15, 0.15);
        IHistogram1D toptrkPy = aida.histogram1D("Top Track Momentum (Py)", 100, -0.15, 0.15);
        IHistogram1D toptrkPz = aida.histogram1D("Top Track Momentum (Pz)", 100, 0, 1.5);
        IHistogram1D toptrkChi2 = aida.histogram1D("Top Track Chi2", 100, 0, 100.0);

        IHistogram1D bottrkPx = aida.histogram1D("Bottom Track Momentum (Px)", 100, -0.15, 0.15);
        IHistogram1D bottrkPy = aida.histogram1D("Bottom Track Momentum (Py)", 100, -0.15, 0.15);
        IHistogram1D bottrkPz = aida.histogram1D("Bottom Track Momentum (Pz)", 100, 0, 1.5);
        IHistogram1D bottrkChi2 = aida.histogram1D("Bottom Track Chi2", 100, 0, 100.0);

        IHistogram1D trkd0 = aida.histogram1D("d0 ", 100, -10.0, 10.0);
        IHistogram1D trkphi = aida.histogram1D("sinphi ", 100, -0.2, 0.2);
        IHistogram1D trkomega = aida.histogram1D("omega ", 100, -0.001, 0.001);
        IHistogram1D trklam = aida.histogram1D("tan(lambda) ", 100, -0.1, 0.1);
        IHistogram1D trkz0 = aida.histogram1D("z0 ", 100, -4.0, 4.0);

        IHistogram1D toptrkd0 = aida.histogram1D("d0 Top", 100, -10.0, 10.0);
        IHistogram1D toptrkphi = aida.histogram1D("sinphi Top", 100, -0.2, 0.2);
        IHistogram1D toptrkomega = aida.histogram1D("omega Top", 100, -0.001, 0.001);
        IHistogram1D toptrklam = aida.histogram1D("tan(lambda) Top", 100, 0, 0.1);
        IHistogram1D toptrkz0 = aida.histogram1D("z0 Top", 100, -4.0, 4.0);

        IHistogram1D bottrkd0 = aida.histogram1D("d0 Bottom", 100, -10.0, 10.0);
        IHistogram1D bottrkphi = aida.histogram1D("sinphi Bottom", 100, -0.2, 0.2);
        IHistogram1D bottrkomega = aida.histogram1D("omega Bottom", 100, -0.001, 0.001);
        IHistogram1D bottrklam = aida.histogram1D("tan(lambda) Bottom", 100, 0, 0.1);
        IHistogram1D bottrkz0 = aida.histogram1D("z0 Bottom", 100, -4.0, 4.0);

        IHistogram1D nTracksBot = aida.histogram1D("Tracks per Event Bot", 10, 0, 10);
        IHistogram1D nTracksTop = aida.histogram1D("Tracks per Event Top", 10, 0, 10);
        IHistogram1D nHitsTop = aida.histogram1D("Hits per Track Top", 4, 3, 7);
        IHistogram1D nHitsBot = aida.histogram1D("Hits per Track Bottom", 4, 3, 7);
        IHistogram1D nHits = aida.histogram1D("Hits per Track", 4, 3, 7);
        IHistogram1D nTracks = aida.histogram1D("Tracks per Event", 10, 0, 10);

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

        if (doMatchedClusterPlots) {

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
            IHistogram2D topECal_HighEsum = aida.histogram2D("HighEsum TopTrack ECal Cluster Position", 40, -400, 400, 20, -100, 100);
            IHistogram2D botECal_HighEsum = aida.histogram2D("HighEsum BottomTrack ECal Cluster Position", 40, -400, 400, 20, -100, 100);
            IHistogram2D topECal_LowEsum = aida.histogram2D("LowEsum TopTrack ECal Cluster Position", 40, -400, 400, 20, -100, 100);
            IHistogram2D botECal_LowEsum = aida.histogram2D("LowEsum BottomTrack ECal Cluster Position", 40, -400, 400, 20, -100, 100);
            //            IHistogram2D botECal = aida.histogram2D("Bottom ECal Cluster Position", 50, -400, 400, 10, -100, 0);
            //            IHistogram2D topECal1 = aida.histogram2D("Top ECal Cluster Position (>0 tracks)", 50, -400, 400, 10, 0, 100);
            //            IHistogram2D botECal1 = aida.histogram2D("Bottom ECal Cluster Position (>0 tracks)", 50, -400, 400, 10, -100, 0);
            //            IHistogram2D topECal2 = aida.histogram2D("Top ECal Cluster Position (E>0.1,>0 tracks)", 50, -400, 400, 10, 0, 100);
            //            IHistogram2D botECal2 = aida.histogram2D("Bottom ECal Cluster Position (E>0.1,>0 tracks)", 50, -400, 400, 10, -100, 0);
            //            IHistogram2D topECal3 = aida.histogram2D("Top ECal Cluster Position w_E (E>0.1,>0 tracks)", 50, -400, 400, 10, 0, 100);
            //            IHistogram2D botECal3 = aida.histogram2D("Bottom ECal Cluster Position w_E (E>0.1,>0 tracks)", 50, -400, 400, 10, -100, 0);
            //
            //            IHistogram1D topECalE = aida.histogram1D("Top ECal Cluster Energy", 50, 0, 2);
            //            IHistogram1D botECalE = aida.histogram1D("Bottom ECal Cluster Energy", 50, 0, 2);
            //            IHistogram1D topECalN = aida.histogram1D("Number of Clusters Top", 6, 0, 6);
            //            IHistogram1D botECalN = aida.histogram1D("Number of Clusters Bot", 6, 0, 6);

        }

        if (doOccupancyPlots) {
            for (HpsSiSensor sensor : sensors) {
                aida.histogram1D(sensor.getName() + " - Occupancy", 640, 0, 640);
                aida.histogram1D(sensor.getName() + " - Max Sample Number", 6, -0.5, 5.5);
                aida.histogram1D(sensor.getName() + " - Occupancy Pairs1", 640, 0, 640);
                aida.histogram1D(sensor.getName() + " - Max Sample Number Pairs1", 6, -0.5, 5.5);
                aida.histogram1D(sensor.getName() + " - Occupancy Singles1", 640, 0, 640);
                aida.histogram1D(sensor.getName() + " - Max Sample Number Singles1", 6, -0.5, 5.5);
                aida.histogram1D(sensor.getName() + " - Occupancy Pulser", 640, 0, 640);
                aida.histogram1D(sensor.getName() + " - Max Sample Number Pulser", 6, -0.5, 5.5);
            }
        }

        if (doHitsOnTrackPlots) {
            int i = 0;
            for (SiSensor sensor : sensors) {
                IHistogram1D resX = aida.histogram1D(sensor.getName() + " strip hits on track", 50, 0, 5);
                i++;
            }
        }

        if (doComparisonPlots) {
            aida.histogram2D("HelicalTrackHits Per Event New vs Old", 50, 0, 50, 50, 0, 50);
            aida.histogram2D("New Shared Tracks with New vs Old - Per Event", 10, 0, 10, 10, 0, 10);
            aida.histogram2D("New Partial Tracks with New vs Old - Per Event", 10, 0, 10, 10, 0, 10);
            //            aida.histogram1D("New Shared Tracks Per Event", 5, 0, 5);
            //           aida.histogram1D("New Partial Tracks Per Event", 5, 0, 5);
            aida.histogram1D("New Duplicate Tracks Per Event", 10, 0, 10);
        }

        if (doReconParticlePlots) {
            IHistogram1D botRecoPartETrkPz = aida.histogram1D("Reco Pairs Electrons Bottom Track Pz", 100, 0, 1.5);
            IHistogram1D topRecoPartETrkPz = aida.histogram1D("Reco Pairs Electrons Top Track Pz", 100, 0, 1.5);
            IHistogram1D botRecoPartPTrkPz = aida.histogram1D("Reco Pairs Positrons Bottom Track Pz", 100, 0, 1.5);
            IHistogram1D topRecoPartPTrkPz = aida.histogram1D("Reco Pairs Positrons Top Track Pz", 100, 0, 1.5);

            //            IHistogram1D botRecoPartETrkChi2 = aida.histogram1D("Reco Pairs Electrons Bottom Track Chi2", 100, 0, 100);
            //            IHistogram1D topRecoPartETrkChi2 = aida.histogram1D("Reco Pairs Electrons Top Track Chi2", 100, 0, 100);
            //            IHistogram1D botRecoPartPTrkChi2 = aida.histogram1D("Reco Pairs Positrons Bottom Track Chi2", 100, 0, 100);
            //            IHistogram1D topRecoPartPTrkChi2 = aida.histogram1D("Reco Pairs Positrons Top Track Chi2", 100, 0, 100);

            //            IHistogram1D botRecoPartTrkChi2 = aida.histogram1D("Reco Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoPartTrkChi2 = aida.histogram1D("Reco Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoPartTrkPz = aida.histogram1D("Reco Bottom Track Pz", 100, 0, 1.5);
            //            IHistogram1D topRecoPartTrkPz = aida.histogram1D("Reco Top Track Pz", 100, 0, 1.5);

            //IHistogram1D botRecoPartTrkPz_Esum = aida.histogram1D("Reco Bottom Track Pz HighEsumEvent", 100, 0, 1.5);
            //IHistogram1D topRecoPartTrkPz_Esum = aida.histogram1D("Reco Top Track Pz HighEsumEvent", 100, 0, 1.5);

            IHistogram1D botRecoPartTrkPz_HighEclusTop = aida.histogram1D("Reco Bottom Track Pz HighEclusTopEvent", 50, 0, 1.5);
            IHistogram1D topRecoPartTrkPz_HighEclusTop = aida.histogram1D("Reco Top Track Pz HighEclusTopEvent", 50, 0, 1.5);
            IHistogram1D botRecoPartTrkPz_HighEclusBot = aida.histogram1D("Reco Bottom Track Pz HighEclusBotEvent", 50, 0, 1.5);
            IHistogram1D topRecoPartTrkPz_HighEclusBot = aida.histogram1D("Reco Top Track Pz HighEclusBotEvent", 50, 0, 1.5);
            IHistogram1D botRecoPartTrkChi2_HighEclusTop = aida.histogram1D("Reco Bottom Track Chi2 HighEclusTopEvent", 50, 0, 100);
            IHistogram1D topRecoPartTrkChi2_HighEclusTop = aida.histogram1D("Reco Top Track Chi2 HighEclusTopEvent", 50, 0, 100);
            IHistogram1D botRecoPartTrkChi2_HighEclusBot = aida.histogram1D("Reco Bottom Track Chi2 HighEclusBotEvent", 50, 0, 100);
            IHistogram1D topRecoPartTrkChi2_HighEclusBot = aida.histogram1D("Reco Top Track Chi2 HighEclusBotEvent", 50, 0, 100);
            IHistogram1D botRecoPartTrkPz_LowEclusTop = aida.histogram1D("Reco Bottom Track Pz LowEclusTopEvent", 50, 0, 1.5);
            IHistogram1D topRecoPartTrkPz_LowEclusTop = aida.histogram1D("Reco Top Track Pz LowEclusTopEvent", 50, 0, 1.5);
            IHistogram1D botRecoPartTrkPz_LowEclusBot = aida.histogram1D("Reco Bottom Track Pz LowEclusBotEvent", 50, 0, 1.5);
            IHistogram1D topRecoPartTrkPz_LowEclusBot = aida.histogram1D("Reco Top Track Pz LowEclusBotEvent", 50, 0, 1.5);
            IHistogram1D botRecoPartTrkChi2_LowEclusTop = aida.histogram1D("Reco Bottom Track Chi2 LowEclusTopEvent", 50, 0, 100);
            IHistogram1D topRecoPartTrkChi2_LowEclusTop = aida.histogram1D("Reco Top Track Chi2 LowEclusTopEvent", 50, 0, 100);
            IHistogram1D botRecoPartTrkChi2_LowEclusBot = aida.histogram1D("Reco Bottom Track Chi2 LowEclusBotEvent", 50, 0, 100);
            IHistogram1D topRecoPartTrkChi2_LowEclusBot = aida.histogram1D("Reco Top Track Chi2 LowEclusBotEvent", 50, 0, 100);

            IHistogram1D top3DHitSVTlayers = aida.histogram1D("Top Layers with 3D-Hits in Events with No Top Track", 7, -0.5, 6.5);
            IHistogram1D bot3DHitSVTlayers = aida.histogram1D("Bottom Layers with 3D-Hits in Events with No Bottom Track", 7, -0.5, 6.5);
            IHistogram1D top2DHitSVTlayers = aida.histogram1D("Top Layers with 2D-Hits in Events with No Top Track", 13, -0.5, 12.5);
            IHistogram1D bot2DHitSVTlayers = aida.histogram1D("Bottom Layers with 2D-Hits in Events with No Bottom Track", 13, -0.5, 12.5);

            IHistogram1D top3DHitSVTtimes = aida.histogram1D("Times of Top 3D-Hits in Events with No Top Track", 30, -15, 15);
            IHistogram1D bot3DHitSVTtimes = aida.histogram1D("Times of Bottom 3D-Hits in Events with No Bottom Track", 30, -15, 15);
            IHistogram1D top2DHitSVTtimes = aida.histogram1D("Times of Top 2D-Hits in Events with No Top Track", 100, -100, 100);
            IHistogram1D bot2DHitSVTtimes = aida.histogram1D("Times of Bottom 2D-Hits in Events with No Bottom Track", 100, -100, 100);
            IHistogram1D top2DHitSVTtimes_ele = aida.histogram1D("Times of Top 2D-Hits in Electron Events with No Top Track", 100, -100, 100);
            IHistogram1D bot2DHitSVTtimes_ele = aida.histogram1D("Times of Bottom 2D-Hits in Electron Events with No Bottom Track", 100, -100, 100);
            IHistogram1D top2DHitSVTtimes_pos = aida.histogram1D("Times of Top 2D-Hits in Positron Events with No Top Track", 100, -100, 100);
            IHistogram1D bot2DHitSVTtimes_pos = aida.histogram1D("Times of Bottom 2D-Hits in Positron Events with No Bottom Track", 100, -100, 100);
            IHistogram1D top2DHitSVTtimes_high = aida.histogram1D("Times of Top 2D-Hits in HighMaxClusE Events with No Top Track", 100, -100, 100);
            IHistogram1D bot2DHitSVTtimes_high = aida.histogram1D("Times of Bottom 2D-Hits in HighMaxClusE Events with No Bottom Track", 100, -100, 100);
            IHistogram1D top2DHitSVTtimes_low = aida.histogram1D("Times of Top 2D-Hits in LowMaxClusE Events with No Top Track", 100, -100, 100);
            IHistogram1D bot2DHitSVTtimes_low = aida.histogram1D("Times of Bottom 2D-Hits in LowMaxClusE Events with No Bottom Track", 100, -100, 100);
            IHistogram1D top2DHitSVTtimes_lowBot = aida.histogram1D("Times of Top 2D-Hits in LowMaxClusEBot Events with No Top Track", 100, -100, 100);
            IHistogram1D bot2DHitSVTtimes_lowBot = aida.histogram1D("Times of Bottom 2D-Hits in LowMaxClusEBot Events with No Bottom Track", 100, -100, 100);
            IHistogram1D top2DHitSVTtimes_lowTop = aida.histogram1D("Times of Top 2D-Hits in LowMaxClusETop Events with No Top Track", 100, -100, 100);
            IHistogram1D bot2DHitSVTtimes_lowTop = aida.histogram1D("Times of Bottom 2D-Hits in LowMaxClusETop Events with No Bottom Track", 100, -100, 100);
            for (int i = 1; i <= 12; i++) {
                String temp = String.format("Times of Top 2D-Hits in Events with No Top Track - Layer %d", i);
                aida.histogram1D(temp, 100, -100, 100);
                temp = String.format("Times of Bottom 2D-Hits in Events with No Bottom Track - Layer %d", i);
                aida.histogram1D(temp, 100, -100, 100);
            }
            for (int i = 0; i < 6; i++) {
                String temp = String.format("Times of Top 2D-Hits in Events with No Top Track - Trigger Phase %d", i);
                aida.histogram1D(temp, 100, -100, 100);
                temp = String.format("Times of Bottom 2D-Hits in Events with No Bottom Track - Trigger Phase %d", i);
                aida.histogram1D(temp, 100, -100, 100);
            }

            IHistogram2D clusE_2D_both = aida.histogram2D("Top vs Bottom ClusE with tracks in both", 25, 0, 0.75, 25, 0, 0.75);
            IHistogram2D clusE_2D_top = aida.histogram2D("Top vs Bottom ClusE with track in top", 25, 0, 0.75, 25, 0, 0.75);
            IHistogram2D clusE_2D_bot = aida.histogram2D("Top vs Bottom ClusE with track in bottom", 25, 0, 0.75, 25, 0, 0.75);
            IHistogram2D clusE_2D_none = aida.histogram2D("Top vs Bottom ClusE with no tracks", 25, 0, 0.75, 25, 0, 0.75);

            IHistogram1D ESum_Reco_both = aida.histogram1D("Reco ESum with tracks in both", 50, 0, 1.5);
            IHistogram1D ESum_Reco_top = aida.histogram1D("Reco ESum with track in top", 50, 0, 1.5);
            IHistogram1D ESum_Reco_bot = aida.histogram1D("Reco ESum with track in bottom", 50, 0, 1.5);
            IHistogram1D ESum_Reco_none = aida.histogram1D("Reco ESum with no tracks", 50, 0, 1.5);
            IHistogram1D ESum_Trigger_both = aida.histogram1D("Trigger ESum with tracks in both", 50, 0, 1.5);
            IHistogram1D ESum_Trigger_top = aida.histogram1D("Trigger ESum with track in top", 50, 0, 1.5);
            IHistogram1D ESum_Trigger_bot = aida.histogram1D("Trigger ESum with track in bottom", 50, 0, 1.5);
            IHistogram1D ESum_Trigger_none = aida.histogram1D("Trigger ESum with no tracks", 50, 0, 1.5);

            IHistogram2D topECalExtrap_HighEsum = aida.histogram2D("HighEsum TopTrack ECal Extrap Position", 40, -400, 400, 20, -100, 100);
            IHistogram2D botECalExtrap_HighEsum = aida.histogram2D("HighEsum BottomTrack ECal Extrap Position", 40, -400, 400, 20, -100, 100);
            IHistogram2D topECalExtrap_LowEsum = aida.histogram2D("LowEsum TopTrack ECal Extrap Position", 40, -400, 400, 20, -100, 100);
            IHistogram2D botECalExtrap_LowEsum = aida.histogram2D("LowEsum BottomTrack ECal Extrap Position", 40, -400, 400, 20, -100, 100);

            IHistogram2D botClusE_TrackPz_LowEsum = aida.histogram2D("LowEsum BottomTrack ClusterE vs TrackPz", 50, 0, 1.5, 25, 0, 0.75);
            IHistogram2D topClusE_TrackPz_LowEsum = aida.histogram2D("LowEsum TopTrack ClusterE vs TrackPz", 50, 0, 1.5, 25, 0, 0.75);
            //            IHistogram1D botRecoPartTrkPz_LowEclusBot = aida.histogram1D("Reco Particles Bottom Track Pz LowEclusBotEvent", 100, 0, 1.5);
            //            IHistogram1D topRecoPartTrkPz_LowEclusBot = aida.histogram1D("Reco Particles Top Track Pz LowEclusBotEvent", 100, 0, 1.5);

            //            IHistogram1D botRecoPartTrkPz_HighEclusTop_3 = aida.histogram1D("Reco Particles Bottom Track Pz HighEclusTopEvent minTracks=3", 100, 0, 1.5);
            //            IHistogram1D topRecoPartTrkPz_HighEclusTop_3 = aida.histogram1D("Reco Particles Top Track Pz HighEclusTopEvent minTracks=3", 100, 0, 1.5);
            //            IHistogram1D botRecoPartTrkPz_HighEclusBot_3 = aida.histogram1D("Reco Particles Bottom Track Pz HighEclusBotEvent minTracks=3", 100, 0, 1.5);
            //            IHistogram1D topRecoPartTrkPz_HighEclusBot_3 = aida.histogram1D("Reco Particles Top Track Pz HighEclusBotEvent minTracks=3", 100, 0, 1.5);
            //            IHistogram1D botRecoPartTrkPz_LowEclusBot_3 = aida.histogram1D("Reco Particles Bottom Track Pz LowEclusBotEvent minTracks=3", 100, 0, 1.5);
            //            IHistogram1D topRecoPartTrkPz_LowEclusBot_3 = aida.histogram1D("Reco Particles Top Track Pz LowEclusBotEvent minTracks=3", 100, 0, 1.5);
            //
            //            IHistogram1D botRecoPartTrkPz_HighEclusTop_2 = aida.histogram1D("Reco Particles Bottom Track Pz HighEclusTopEvent Tracks=2", 100, 0, 1.5);
            //            IHistogram1D topRecoPartTrkPz_HighEclusTop_2 = aida.histogram1D("Reco Particles Top Track Pz HighEclusTopEvent Tracks=2", 100, 0, 1.5);
            //            IHistogram1D botRecoPartTrkPz_HighEclusBot_2 = aida.histogram1D("Reco Particles Bottom Track Pz HighEclusBotEvent Tracks=2", 100, 0, 1.5);
            //            IHistogram1D topRecoPartTrkPz_HighEclusBot_2 = aida.histogram1D("Reco Particles Top Track Pz HighEclusBotEvent Tracks=2", 100, 0, 1.5);
            //            IHistogram1D botRecoPartTrkPz_LowEclusBot_2 = aida.histogram1D("Reco Particles Bottom Track Pz LowEclusBotEvent Tracks=2", 100, 0, 1.5);
            //            IHistogram1D topRecoPartTrkPz_LowEclusBot_2 = aida.histogram1D("Reco Particles Top Track Pz LowEclusBotEvent Tracks=2", 100, 0, 1.5);
            //
            //            IHistogram1D botRecoPartTrkPz_HighEclusTop_1 = aida.histogram1D("Reco Particles Bottom Track Pz HighEclusTopEvent Tracks=1", 100, 0, 1.5);
            //            IHistogram1D topRecoPartTrkPz_HighEclusTop_1 = aida.histogram1D("Reco Particles Top Track Pz HighEclusTopEvent Tracks=1", 100, 0, 1.5);
            //            IHistogram1D botRecoPartTrkPz_HighEclusBot_1 = aida.histogram1D("Reco Particles Bottom Track Pz HighEclusBotEvent Tracks=1", 100, 0, 1.5);
            //            IHistogram1D topRecoPartTrkPz_HighEclusBot_1 = aida.histogram1D("Reco Particles Top Track Pz HighEclusBotEvent Tracks=1", 100, 0, 1.5);
            //            IHistogram1D botRecoPartTrkPz_LowEclusBot_1 = aida.histogram1D("Reco Particles Bottom Track Pz LowEclusBotEvent Tracks=1", 100, 0, 1.5);
            //            IHistogram1D topRecoPartTrkPz_LowEclusBot_1 = aida.histogram1D("Reco Particles Top Track Pz LowEclusBotEvent Tracks=1", 100, 0, 1.5);
            //
            //            IHistogram1D botRecoPairsPartETrkChi2 = aida.histogram1D("Reco Pairs Electrons Bottom Track Pz", 100, 0, 1.5);
            //            IHistogram1D topRecoPairsPartETrkChi2 = aida.histogram1D("Reco Pairs Electrons Top Track Pz", 100, 0, 1.5);
            //            IHistogram1D botRecoPairsPartPTrkChi2 = aida.histogram1D("Reco Pairs Positrons Bottom Track Pz", 100, 0, 1.5);
            //            IHistogram1D topRecoPairsPartPTrkChi2 = aida.histogram1D("Reco Pairs Positrons Top Track Pz", 100, 0, 1.5);
            //            IHistogram1D botRecoPairsPartTrkChi2 = aida.histogram1D("Reco Pairs Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoPairsPartTrkChi2 = aida.histogram1D("Reco Pairs Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoPairsPartTrkPz = aida.histogram1D("Reco Pairs Top Track Pz", 100, 0, 1.5);
            //            IHistogram1D topRecoPairsPartTrkPz = aida.histogram1D("Reco Pairs Bottom Track Pz", 100, 0, 1.5);

            //            IHistogram1D botRecoTrkChi2_energetic = aida.histogram1D("Reco Energetic Particles Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoTrkChi2_energetic = aida.histogram1D("Reco Energetic Particles Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoTrkPz_energetic = aida.histogram1D("Reco Energetic Particles Bottom Track Pz", 100, 0, 1.5);
            //            IHistogram1D topRecoTrkPz_energetic = aida.histogram1D("Reco Energetic Particles Top Track Pz", 100, 0, 1.5);
            //
            //            IHistogram1D botRecoMatchedTrkPz = aida.histogram1D("Reco MatchedPairs Particles Bottom Track Pz", 100, 0, 1.5);
            //            IHistogram1D topRecoMatchedTrkPz = aida.histogram1D("Reco MatchedPairs Particles Top Track Pz", 100, 0, 1.5);

            //
            //            IHistogram1D botRecoFEETrigTrkChi2 = aida.histogram1D("Reco FEE Particles Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoFEETrigTrkChi2 = aida.histogram1D("Reco FEE Particles Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoFEETrigTrkPz = aida.histogram1D("Reco FEE Particles Bottom Track Pz", 100, 0, 1.5);
            //            IHistogram1D topRecoFEETrigTrkPz = aida.histogram1D("Reco FEE Particles Top Track Pz", 100, 0, 1.5);
            //            IHistogram1D botRecoFEEMatchedTrkChi2 = aida.histogram1D("Reco FEEMatched Particles Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoFEEMatchedTrkChi2 = aida.histogram1D("Reco FEEMatched Particles Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoFEEMatchedTrkPz = aida.histogram1D("Reco FEEMatched Particles Bottom Track Pz", 100, 0, 1.5);
            //            IHistogram1D topRecoFEEMatchedTrkPz = aida.histogram1D("Reco FEEMatched Particles Top Track Pz", 100, 0, 1.5);
            //
            //            IHistogram1D botRecoPartFEELowPzETrkChi2 = aida.histogram1D("Reco FEELowPz Electrons Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoPartFEELowPzETrkChi2 = aida.histogram1D("Reco FEELowPz Electrons Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoPartFEELowPzPTrkChi2 = aida.histogram1D("Reco FEELowPz Positrons Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoPartFEELowPzPTrkChi2 = aida.histogram1D("Reco FEELowPz Positrons Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoPartFEELowPzTrkChi2 = aida.histogram1D("Reco FEELowPz Particles Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoPartFEELowPzTrkChi2 = aida.histogram1D("Reco FEELowPz Particles Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoPartFEELowPzTrkd0 = aida.histogram1D("Reco FEELowPz Particles Bottom Track d0", 100, -10, 10);
            //            IHistogram1D topRecoPartFEELowPzTrkd0 = aida.histogram1D("Reco FEELowPz Particles Top Track d0", 100, -10, 10);
            //            IHistogram1D botRecoPartFEELowPzTrkz0 = aida.histogram1D("Reco FEELowPz Particles Bottom Track z0", 100, -5, 5);
            //            IHistogram1D topRecoPartFEELowPzTrkz0 = aida.histogram1D("Reco FEELowPz Particles Top Track z0", 100, -5, 5);
            //
            //            IHistogram1D botRecoPartPairsLowPzETrkChi2 = aida.histogram1D("Reco PairsLowPz Electrons Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoPartPairsLowPzETrkChi2 = aida.histogram1D("Reco PairsLowPz Electrons Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoPartPairsLowPzPTrkChi2 = aida.histogram1D("Reco PairsLowPz Positrons Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoPartPairsLowPzPTrkChi2 = aida.histogram1D("Reco PairsLowPz Positrons Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoPartPairsLowPzTrkChi2 = aida.histogram1D("Reco PairsLowPz Particles Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoPartPairsLowPzTrkChi2 = aida.histogram1D("Reco PairsLowPz Particles Top Track Chi2", 100, 0, 100.0);
            //
            //            IHistogram1D botRecoPartPairsLowPzMatchedETrkChi2 = aida.histogram1D("Reco PairsLowPzMatched Electrons Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoPartPairsLowPzMatchedETrkChi2 = aida.histogram1D("Reco PairsLowPzMatched Electrons Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoPartPairsLowPzMatchedPTrkChi2 = aida.histogram1D("Reco PairsLowPzMatched Positrons Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoPartPairsLowPzMatchedPTrkChi2 = aida.histogram1D("Reco PairsLowPzMatched Positrons Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoPartPairsLowPzMatchedTrkChi2 = aida.histogram1D("Reco PairsLowPzMatched Particles Bottom Track Chi2", 100, 0, 100.0);
            //            IHistogram1D topRecoPartPairsLowPzMatchedTrkChi2 = aida.histogram1D("Reco PairsLowPzMatched Particles Top Track Chi2", 100, 0, 100.0);
            //            IHistogram1D botRecoPartPairsLowPzMatchedTrkd0 = aida.histogram1D("Reco PairsLowPzMatched Particles Bottom Track d0", 100, -10, 10);
            //            IHistogram1D topRecoPartPairsLowPzMatchedTrkd0 = aida.histogram1D("Reco PairsLowPzMatched Particles Top Track d0", 100, -10, 10);
            //            IHistogram1D botRecoPartPairsLowPzMatchedTrkz0 = aida.histogram1D("Reco PairsLowPzMatched Particles Bottom Track z0", 100, -5, 5);
            //            IHistogram1D topRecoPartPairsLowPzMatchedTrkz0 = aida.histogram1D("Reco PairsLowPzMatched Particles Top Track z0", 100, -5, 5);
            //
            //            IHistogram2D botbotRecoPartFEELowPz2D = aida.histogram2D("Reco FEELowPz vs MostEnergetic bottom-bottom", 100, -1.5, 1.5, 100, -1.5, 1.5);
            //            IHistogram2D bottopRecoPartFEELowPz2D = aida.histogram2D("Reco FEELowPz vs MostEnergetic bottom-top", 100, -1.5, 1.5, 100, -1.5, 1.5);
            //            IHistogram2D topbotRecoPartFEELowPz2D = aida.histogram2D("Reco FEELowPz vs MostEnergetic top-bottom", 100, -1.5, 1.5, 100, -1.5, 1.5);
            //            IHistogram2D toptopRecoPartFEELowPz2D = aida.histogram2D("Reco FEELowPz vs MostEnergetic top-top", 100, -1.5, 1.5, 100, -1.5, 1.5);
        }
    }
}
