package org.hps.recon.tracking.gbl;

import hep.physics.matrix.SymmetricMatrix;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;
import org.hps.recon.tracking.gbl.matrix.SymMatrix;
import org.hps.recon.tracking.gbl.matrix.Vector;


import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.io.IOException;

import org.apache.commons.math3.util.Pair;
import org.hps.recon.tracking.CoordinateTransformations;
import org.hps.recon.tracking.MultipleScattering;
import org.hps.recon.tracking.TrackType;
import org.hps.recon.tracking.TrackUtils;
import org.lcsim.constants.Constants;
import org.lcsim.detector.ITransform3D;
import org.lcsim.detector.IDetectorElement;
import org.lcsim.detector.tracker.silicon.ChargeCarrier;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.detector.tracker.silicon.SiSensorElectrodes;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
//import org.hps.recon.tracking.TrackStateUtils;
import org.lcsim.event.TrackerHit;
import org.lcsim.event.base.BaseTrack;
import org.lcsim.event.base.BaseTrackState;
import org.lcsim.fit.helicaltrack.HelicalTrackFit;
import org.lcsim.fit.helicaltrack.HelicalTrackStrip;
import org.lcsim.fit.helicaltrack.HelixUtils;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;
import org.lcsim.recon.tracking.digitization.sisim.TrackerHitType;
import org.lcsim.recon.tracking.seedtracker.ScatterAngle;

import org.lcsim.util.aida.AIDA;

/**
 * Utilities that create track objects from fitted GBL trajectories.
 *
 * @author Per Hansson Adrian <phansson@slac.stanford.edu>
 * @author Miriam Diamond
 * @author PF
 *
 */
public class MakeGblTracks {

    private Logger LOGGER = Logger.getLogger(MakeGblTracks.class.getPackage().getName());
    public AIDA aida = null;
    private String outputPlots = "RefitterPlots.root";
    private String plotDir = "RefitterPlots";
    /*
    static {
        LOGGER.setLevel(Level.WARNING);
    }
    */
    public MakeGblTracks() {
        aida = AIDA.defaultInstance();
        //aida.tree().mkdir(plotDir);
        //aida.tree().cd(plotDir);
        aida.tree().cd("/");
        setupPlots();
    }

    public  void setDebug(boolean debug) {
        if (debug) {
            LOGGER.setLevel(Level.INFO);
        } else {
            LOGGER.setLevel(Level.OFF);
        }
    }

    private void setupPlots() {
        int niter = 20;
        aida.histogram1D("timing",100,0,100);
        for (int i = 0 ; i<15; i++) { 
            //For each label (trackState I have 5 corrections)
            aida.histogram1D("Correction_q_p_l"  + String.valueOf(i),niter,0,niter);
            aida.histogram1D("Correction_yTp_l"   + String.valueOf(i),niter,0,niter);
            aida.histogram1D("Correction_xTp_l" + String.valueOf(i),niter,0,niter);
            aida.histogram1D("Correction_xT_l" + String.valueOf(i),niter,0,niter);
            aida.histogram1D("Correction_yT_l"   + String.valueOf(i),niter,0,niter);

            for (int it = 0 ; it<niter; it++) {
                aida.histogram1D("Correction_q_p_l"   + String.valueOf(i)+"it_"+String.valueOf(it),100,-0.05,0.05);
                aida.histogram1D("Correction_yTp_l"    + String.valueOf(i)+"it_"+String.valueOf(it),100,-0.05,0.05);
                aida.histogram1D("Correction_xTp_l"  + String.valueOf(i)+"it_"+String.valueOf(it),100,-0.05,0.05);
                aida.histogram1D("Correction_xT_l"  + String.valueOf(i)+"it_"+String.valueOf(it),100,-3,3);
                aida.histogram1D("Correction_yT_l"    + String.valueOf(i)+"it_"+String.valueOf(it),100,-3,3);
            }
        }
    }
    
    private void FillIterationHistos(FittedGblTrajectory fit, int iter) { 
        Vector locPar = new Vector(5);
        SymMatrix locCov = new SymMatrix(5);
        
        Integer [] sensorsFromMapArray = fit.getSensorMap().keySet().toArray(new Integer[0]);
        int sensors_size = fit.getSensorMap().keySet().size();
        
        for (int i = 0; i<sensors_size; i++) { 
            int iLabel = sensorsFromMapArray[i];
            fit.getResults(iLabel,locPar,locCov);
            aida.histogram1D("Correction_q_p_l"  +String.valueOf(i)).fill(iter,Math.abs(locPar.get(0)));
            aida.histogram1D("Correction_yTp_l"   +String.valueOf(i)).fill(iter,Math.abs(locPar.get(1)));
            aida.histogram1D("Correction_xTp_l" +String.valueOf(i)).fill(iter,Math.abs(locPar.get(2)));
            aida.histogram1D("Correction_xT_l" +String.valueOf(i)).fill(iter,Math.abs(locPar.get(3)));
            aida.histogram1D("Correction_yT_l"   +String.valueOf(i)).fill(iter,Math.abs(locPar.get(4)));
            
            aida.histogram1D("Correction_q_p_l"  + String.valueOf(i)+"it_"+String.valueOf(iter)).fill(locPar.get(0));
            aida.histogram1D("Correction_yTp_l"   + String.valueOf(i)+"it_"+String.valueOf(iter)).fill(locPar.get(1));
            aida.histogram1D("Correction_xTp_l" + String.valueOf(i)+"it_"+String.valueOf(iter)).fill(locPar.get(2));
            aida.histogram1D("Correction_xT_l" + String.valueOf(i)+"it_"+String.valueOf(iter)).fill(locPar.get(3));
            aida.histogram1D("Correction_yT_l"   + String.valueOf(i)+"it_"+String.valueOf(iter)).fill(locPar.get(4));
        }
    }
    
    public void savePlots() {
        
        try { 
            aida.saveAs(outputPlots);
        } catch  (IOException ex) {
            
        }
        
    }
    
    

    public Pair<Track, GBLKinkData> makeCorrectedTrack(FittedGblTrajectory fittedGblTrajectory, HelicalTrackFit helicalTrackFit, List<TrackerHit> hitsOnTrack, int trackType, double bfield) {
        return makeCorrectedTrack(fittedGblTrajectory, helicalTrackFit, hitsOnTrack, trackType, bfield, false);
    }


    public Pair<Track,GBLKinkData> makeCorrectedTrack(FittedGblTrajectory fittedGblTrajectory, Track track, List<TrackerHit> hitsOnTrack, int trackType, double bfield) { 
        
        double[] ref = new double[] { 0., 0., 0. };
        BaseTrack trk = new BaseTrack();
        trk.getTrackStates().clear();
        
        for (TrackerHit hit : hitsOnTrack) { 
            trk.addHit(hit);
        }
        
        // Set state at IP
        TrackState ts_ip = track.getTrackStates().get(0);
        Pair<double[], SymmetricMatrix> correctedHelixParams = fittedGblTrajectory.getCorrectedPerigeeParameters(TrackUtils.getHTF(ts_ip), FittedGblTrajectory.GBLPOINT.IP, bfield);
        trk.setTrackParameters(correctedHelixParams.getFirst(), bfield);// hack to set the track charge
        trk.getTrackStates().clear();
        TrackState stateIP = new BaseTrackState(correctedHelixParams.getFirst(), ref, correctedHelixParams.getSecond().asPackedArray(true), TrackState.AtIP, bfield);
        trk.getTrackStates().add(stateIP);
        int prevID = 0;
        int dummyCounter = -1;
        
        Integer[] sensorsFromMapArray = fittedGblTrajectory.getSensorMap().keySet().toArray(new Integer[0]);
        
        //sensors FromMapArray length = gblStripClusterDataSize  = > easier to keep the TSOS on the gblStripClusterData
        int TScounter = 0;
        for (int i = 0; i < sensorsFromMapArray.length; i++) {
            //The 0 already has been used for the IP - so we start from 1
            TScounter+=1;
            while (track.getTrackStates().get(TScounter).getLocation() <  0 ) {
                //System.out.println("TSCounter = " + TScounter + " negative/dummy shit " + track.getTrackStates().get(TScounter).getLocation());
                TScounter+=1;
            }
            
            //Correct original helix - I need to get the right TSOS. They can be 0
            double[] TSOS =  track.getTrackStates().get(TScounter).getParameters();
            
            int ilabel = sensorsFromMapArray[i];
            int millepedeID = fittedGblTrajectory.getSensorMap().get(ilabel);
            
            // if sensors are missing from track, insert blank TrackState objects
            for (int k = 1; k < millepedeID - prevID; k++) {
                // uses new lcsim constructor
                BaseTrackState dummy = new BaseTrackState(dummyCounter);
                trk.getTrackStates().add(dummy);
                dummyCounter--;
            }
            prevID = millepedeID;
            Pair<double[], SymmetricMatrix> correctedHelixParamsSensor = fittedGblTrajectory.getCorrectedPerigeeParameters(TrackUtils.getHTF(TSOS), ilabel, bfield);
            // set TrackState location code
            int loc = TrackState.AtOther;
            if (i == 0) {
                loc = TrackState.AtFirstHit;
            } else if (i == sensorsFromMapArray.length - 1) {
                loc = TrackState.AtLastHit;
            }
            // insert TrackState at sensor
            TrackState stateSensor = new BaseTrackState(correctedHelixParamsSensor.getFirst(), ref, correctedHelixParamsSensor.getSecond().asPackedArray(true), loc, bfield);
            trk.getTrackStates().add(stateSensor);
        }
        GBLKinkData kinkData = fittedGblTrajectory.getKinks();
        return new Pair<Track, GBLKinkData>(trk, kinkData);
        
    }

        /**
         * Create a new {@link BaseTrack} from a {@link FittedGblTrajectory}.
         *
         * @param fittedGblTrajectory
         * @param helicalTrackFit
         * @param hitsOnTrack
         * @param trackType
         * @param bfield
         * @return the new {@link BaseTrack} and the kinks along the
         * {@link GblTrajectory} as a {@link Pair}.
     */
    public  Pair<Track, GBLKinkData> makeCorrectedTrack(FittedGblTrajectory fittedGblTrajectory, HelicalTrackFit helicalTrackFit, List<TrackerHit> hitsOnTrack, int trackType, double bfield, boolean storeTrackStates) {
        //  Initialize the reference point to the origin
        double[] ref = new double[] { 0., 0., 0. };

        // Create a new SeedTrack
        BaseTrack trk = new BaseTrack();

        // Add the hits to the track
        for (TrackerHit hit : hitsOnTrack) {
            trk.addHit(hit);
        }

        // Set state at IP
        Pair<double[], SymmetricMatrix> correctedHelixParams = fittedGblTrajectory.getCorrectedPerigeeParameters(helicalTrackFit, FittedGblTrajectory.GBLPOINT.IP, bfield);
        trk.setTrackParameters(correctedHelixParams.getFirst(), bfield);// hack to set the track charge
        trk.getTrackStates().clear();
        TrackState stateIP = new BaseTrackState(correctedHelixParams.getFirst(), ref, correctedHelixParams.getSecond().asPackedArray(true), TrackState.AtIP, bfield);
        trk.getTrackStates().add(stateIP);

        if (!storeTrackStates) {
            // just store last state
            Pair<double[], SymmetricMatrix> correctedHelixParamsLast = fittedGblTrajectory.getCorrectedPerigeeParameters(helicalTrackFit, FittedGblTrajectory.GBLPOINT.LAST, bfield);
            TrackState stateLast = new BaseTrackState(correctedHelixParamsLast.getFirst(), ref, correctedHelixParamsLast.getSecond().asPackedArray(true), TrackState.AtLastHit, bfield);
            trk.getTrackStates().add(stateLast);
        } else {

            // store states at all 18 sensors
            int prevID = 0;
            int dummyCounter = -1;
            // note: SensorMap doesn't include IP
            Integer[] sensorsFromMapArray = fittedGblTrajectory.getSensorMap().keySet().toArray(new Integer[0]);

            for (int i = 0; i < sensorsFromMapArray.length; i++) {
                int ilabel = sensorsFromMapArray[i];
                int millepedeID = fittedGblTrajectory.getSensorMap().get(ilabel);

                // if sensors are missing from track, insert blank TrackState objects
                for (int k = 1; k < millepedeID - prevID; k++) {
                    // uses new lcsim constructor
                    BaseTrackState dummy = new BaseTrackState(dummyCounter);
                    trk.getTrackStates().add(dummy);
                    dummyCounter--;
                }
                prevID = millepedeID;
                Pair<double[], SymmetricMatrix> correctedHelixParamsSensor = fittedGblTrajectory.getCorrectedPerigeeParameters(helicalTrackFit, ilabel, bfield);

                // set TrackState location code
                int loc = TrackState.AtOther;
                if (i == 0) {
                    loc = TrackState.AtFirstHit;
                } else if (i == sensorsFromMapArray.length - 1) {
                    loc = TrackState.AtLastHit;
                }
                // insert TrackState at sensor
                TrackState stateSensor = new BaseTrackState(correctedHelixParamsSensor.getFirst(), ref, correctedHelixParamsSensor.getSecond().asPackedArray(true), loc, bfield);
                trk.getTrackStates().add(stateSensor);
            }

        }

        // Extract kinks from trajectory
        GBLKinkData kinkData = fittedGblTrajectory.getKinks();

        // Set other info needed
        trk.setChisq(fittedGblTrajectory.get_chi2());
        trk.setNDF(fittedGblTrajectory.get_ndf());
        trk.setFitSuccess(true);
        trk.setRefPointIsDCA(true);
        trk.setTrackType(TrackType.setGBL(trackType, true));

        LOGGER.fine(String.format("helix chi2 %f ndf %d gbl chi2 %f ndf %d\n", helicalTrackFit.chisqtot(), helicalTrackFit.ndf()[0] + helicalTrackFit.ndf()[1], trk.getChi2(), trk.getNDF()));

        return new Pair<Track, GBLKinkData>(trk, kinkData);
    }

    public  Pair<Track, GBLKinkData> refitTrack(HelicalTrackFit helix, Collection<TrackerHit> stripHits, Collection<TrackerHit> hth, int nIterations, int trackType, MultipleScattering scattering, double bfield) {
        return refitTrack(helix, stripHits, hth, nIterations, trackType, scattering, bfield, false);
    }

    /**
     * Do a GBL fit to an arbitrary set of strip hits, with a starting value of
     * the helix parameters.
     *
     * @param helix Initial helix parameters. Only track parameters are used
     * (not covariance)
     * @param stripHits Strip hits to be used for the GBL fit. Does not need to
     * be in sorted order.
     * @param hth Stereo hits for the track's hit list (these are not used in
     * the GBL fit). Does not need to be in sorted order.
     * @param nIterations Number of times to iterate the GBL fit.
     * @param scattering Multiple scattering manager.
     * @param bfield B-field
     * @return The refitted track.
     */
    public  Pair<Track, GBLKinkData> refitTrack(HelicalTrackFit helix, Collection<TrackerHit> stripHits, Collection<TrackerHit> hth, int nIterations, int trackType, MultipleScattering scattering, double bfield, boolean storeTrackStates) {
        Pair<Pair<Track, GBLKinkData>, FittedGblTrajectory> refitTrack = refitTrackWithTraj(helix, stripHits, hth, nIterations, trackType, scattering, bfield, storeTrackStates,false);
        if(refitTrack == null)
            return null;
        else
            return refitTrack.getFirst();
    }
    
    /*
    public Pair<Pair<Track, GBLKinkData>, FittedGblTrajectory> refitTrackWithTraj_2(HelicalTrackFit helix, Collection<TrackerHit> stripHits, Collection<TrackerHit> hth, int nIterations, int trackType, MultipleScattering scattering, double bfield, boolean storeTrackStates,boolean includeMS) {
        
        double bfac = Constants.fieldConversion * bfield;
        //Initial fit on the seed
        List<TrackerHit> allHthList = TrackUtils.sortHits(hth);
        List<TrackerHit> sortedStripHits = TrackUtils.sortHits(stripHits);
        //List<HpsSensor> sortedSensors    = 
        FittedGblTrajectory fit = doGBLFit(helix, sortedStripHits, scattering, bfield, 0,includeMS);
        if(fit==null) return null;    
        //At this stage I have the FittedGblTrajectory with the GblStripClusterData and the stripHits.
        
        //Create corrected track storing all track states on surface
        Pair<Track, GBLKinkData> newTrack = null;
        
        //Get the sensorArray in the fittedGblTrajector
        Integer [] sensorsFromMapArray = fit.getSensorMap().keySet().toArray(new Integer[0]);
        //Make a cache of the GBLStripClusterData list for next iterations
        List<GBLStripClusterData> updated_gblStripData = new ArrayList<GBLStripClusterData>();
        
        for (int i =0; i<fit.getGBLStripClusterData().size(); i++) {
            updated_gblStripData.add(new GBLStripClusterData(fit.getGBLStripClusterData().get(i)));
        }
        
        for (int iter = 0; iter < nIterations; iter++) { 
            //Create corrected track storing all track states on surface
            long startTime = System.nanoTime();

            newTrack = makeCorrectedTrack(fit, helix, allHthList, trackType, bfield, false);
            helix    = TrackUtils.getHTF(newTrack.getFirst());
                        
            //for (int i =0; i<fit.getGBLStripClusterData().size(); i++) {
            for (int i =0; i<updated_gblStripData.size(); i++) {
                
                //GBLStripClusterData gblStripData = new GBLStripClusterData(fit.getGBLStripClusterData().get(i));
                GBLStripClusterData gblStripData = updated_gblStripData.get(i);
                int iLabel = sensorsFromMapArray[i];
                
                //measurement
                if (gblStripData.getScatterOnly() == 0) {
                    
                    SiTrackerHitStrip1D h = null;
                    //Could cache this too
                    if (sortedStripHits.get(i) instanceof SiTrackerHitStrip1D) {
                        h = (SiTrackerHitStrip1D) sortedStripHits.get(i);
                    }
                    else {
                        h = new SiTrackerHitStrip1D(sortedStripHits.get(i));
                    }
                    
                    
                    HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit) h.getRawHits().get(0)).getDetectorElement();
                                        
                    //Get global position
                    //Hep3Vector ts_extr_pos = TrackStateUtils.getLocationAtSensor(newTrack.getFirst(),sensor,bfield);
                    Hep3Vector ts_extr_pos = TrackStateUtils.getLocationAtSensor(helix,sensor,bfield);
                    
                    //Rotate to tracking
                    Hep3Vector ts_extr_pos_trk = CoordinateTransformations.transformVectorToTracking(ts_extr_pos);
                    //Get s_path and direction
                    double s_corr = HelixUtils.PathToXPlane(helix, ts_extr_pos_trk.x(),0.,0).get(0);
                    Hep3Vector dir_corr = HelixUtils.Direction(helix,s_corr);
                    
                    //Get local position on sensor
                    Hep3Vector ts_extr_pos_l = sensor.getGeometry().getGlobalToLocal().transformed(ts_extr_pos);
                    
                    gblStripData.setPath(s_corr);
                    gblStripData.setPath3D(s_corr / Math.cos(Math.atan(helix.slope())));
                    gblStripData.setTrackDir(dir_corr);
                    gblStripData.setTrackPhi(helix.phi0() - s_corr / helix.R());
                    gblStripData.setTrackLambda(Math.atan(helix.slope()));
                    gblStripData.setTrackPos(ts_extr_pos_l);
                    gblStripData.setTrackPosTrk(ts_extr_pos_trk);
                }
                //updated_gblStripData.add(gblStripData);
            }

            //update fit
            fit = HpsGblRefitter.fit(updated_gblStripData, bfac, false);
            //fit.addGBLStripClusterData(updated_gblStripData);
            //fit.addStripHits(sortedStripHits);
            newTrack = makeCorrectedTrack(fit, helix, allHthList, trackType, bfield, false);
            long endTime = System.nanoTime();
            long duration = (endTime - startTime);//divide by 1000000 to get milliseconds.
            aida.histogram1D("timing").fill((double)duration / 1000000.);
            //Fill the histograms
            FillIterationHistos(fit,iter);
        }
                
        Pair<Track, GBLKinkData> mergedTrack = makeCorrectedTrack(fit,TrackUtils.getHTF(newTrack.getFirst()),allHthList, trackType, bfield, storeTrackStates);
        return new Pair<Pair<Track, GBLKinkData>, FittedGblTrajectory>(mergedTrack, fit);
    }


    public Pair<Pair<Track, GBLKinkData>, FittedGblTrajectory> refitTrackWithTraj_3(HelicalTrackFit helix, Collection<TrackerHit> stripHits, Collection<TrackerHit> hth, int nIterations, int trackType, MultipleScattering scattering, double bfield, boolean storeTrackStates,boolean includeMS) {
        
        double bfac = Constants.fieldConversion * bfield;
        //Initial fit on the seed
        List<TrackerHit> allHthList = TrackUtils.sortHits(hth);
        List<TrackerHit> sortedStripHits = TrackUtils.sortHits(stripHits);
        FittedGblTrajectory fit = doGBLFit(helix, sortedStripHits, scattering, bfield, 0,includeMS);
        if(fit==null) return null;    
        //At this stage I have the FittedGblTrajectory with the GblStripClusterData and the stripHits.
        
        //Create corrected track storing all track states
        Pair<Track, GBLKinkData> newTrack = makeCorrectedTrack(fit, helix, allHthList, trackType, bfield, true);
        
        //Get the sensorArray in the fittedGblTrajector
        Integer [] sensorsFromMapArray = fit.getSensorMap().keySet().toArray(new Integer[0]);
        //System.out.println("SensorsFromMapArray lenght = " + sensorsFromMapArray.length + " gblStripClusterDataSize = " + fit.getGBLStripClusterData().size());
        //Make a cache of the GBLStripClusterData list for next iterations
        List<GBLStripClusterData> updated_gblStripData = new ArrayList<GBLStripClusterData>();
        
        for (int i =0; i<fit.getGBLStripClusterData().size(); i++) {
            updated_gblStripData.add(new GBLStripClusterData(fit.getGBLStripClusterData().get(i)));
        }
        
        FillIterationHistos(fit,0);
        for (int iter = 1; iter < nIterations; iter++) { 
            
            //Create corrected track storing all track states on surface
            long startTime = System.nanoTime();
            
            //System.out.println("SensorsFromMapArray lenght = " + sensorsFromMapArray.length + " gblStripClusterDataSize = " + updated_gblStripData.size() + "TrackStates size=" + newTrack.getFirst().getTrackStates().size());
            //for (int i =0; i<fit.getGBLStripClusterData().size(); i++) {
            int TScounter = -1;
            
            for (int i =0; i<updated_gblStripData.size(); i++) {
                
                TScounter++;
                //GBLStripClusterData gblStripData = new GBLStripClusterData(fit.getGBLStripClusterData().get(i));
                GBLStripClusterData gblStripData = updated_gblStripData.get(i);
                int iLabel = sensorsFromMapArray[i];
                //System.out.println("GblStripData MIP=" + gblStripData.getId() + " SensorsFromMapArray ID = " + fit.getSensorMap().get(iLabel));
                
                
                //measurement
                if (gblStripData.getScatterOnly() == 0) {
                    
                    SiTrackerHitStrip1D h = null;
                    //Could cache this too
                    if (sortedStripHits.get(i) instanceof SiTrackerHitStrip1D) {
                        h = (SiTrackerHitStrip1D) sortedStripHits.get(i);
                    }
                    else {
                        h = new SiTrackerHitStrip1D(sortedStripHits.get(i));
                    }
                    
                    HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit) h.getRawHits().get(0)).getDetectorElement();
                    
                    //Bit wasteful - find better way to make the htf out of track states.
                    //At each stripData, get the corrected track parameters
                    
                    
                    //Correct original helix - I need to get the right TSOS
                    double[] TSOS =  newTrack.getFirst().getTrackStates().get(TScounter).getParameters();
                    while (Math.abs(TSOS[0]) < 0.0000001) {
                        //System.out.println("Track Params tsos: " + TSOS[0] + "," +TSOS[1] + "," + TSOS[2] + "," + TSOS[3] + "," + TSOS[4]); 
                        //System.out.println("TScounter " + TScounter);
                        TScounter+=1;
                        TSOS = newTrack.getFirst().getTrackStates().get(TScounter).getParameters();
                    }
                    
                    //Now I have the TSOS at the right place (actually they are n-1 but let's see if that is OK)
                    HelicalTrackFit corr_htf = TrackUtils.getHTF(TSOS);
                    
                    //System.out.println("iLabel: " + iLabel + " i:"+i + " nIter:"+iter);
                    //System.out.println("Track Params tsos: " + TSOS[0] + "," +TSOS[1] + "," + TSOS[2] + "," + TSOS[3] + "," + TSOS[4]); 
                    
                    
                    //System.out.println("Track Params tsos: " + TSOS[0] + "," +TSOS[1] + "," + TSOS[2] + "," + TSOS[3] + "," + TSOS[4]); 
                    //System.out.println("before pos:" + gblStripData.getTrackPos().x() + "," + gblStripData.getTrackPos().y() + "," + gblStripData.getTrackPos().z());
                    //System.out.println("before dir:" + gblStripData.getTrackDirection().x() + "," + gblStripData.getTrackDirection().y() + "," + gblStripData.getTrackDirection().z());
                    //System.out.println("before s / s3D:" + gblStripData.getPath() + "/"+gblStripData.getPath3D());
                    
                    //Get global position
                    Hep3Vector ts_extr_pos_ts = TrackStateUtils.getLocationAtSensor(newTrack.getFirst(),sensor,bfield);
                    Hep3Vector ts_extr_pos    = TrackStateUtils.getLocationAtSensor(corr_htf,sensor,bfield);
                    
                    //Rotate to tracking
                    Hep3Vector ts_extr_pos_trk = CoordinateTransformations.transformVectorToTracking(ts_extr_pos_ts);
                                        
                    //Iterative in tracking
                    double s_origin = HelixUtils.PathToXPlane(corr_htf, gblStripData.getOrigin().x(), 0., 0).get(0);
                    //Hep3Vector extr_pos_it = TrackUtils.getHelixPlaneIntercept(corr_htf,gblStripData.getW(),gblStripData.getOrigin(),bfield,s_origin);
                                        
                    //Get s_path and direction
                    double s_corr = HelixUtils.PathToXPlane(corr_htf, ts_extr_pos_trk.x(),0.,0).get(0);
                    Hep3Vector dir_corr = HelixUtils.Direction(corr_htf,s_corr);

                    //Get local position on sensor
                    Hep3Vector ts_extr_pos_l = sensor.getGeometry().getGlobalToLocal().transformed(ts_extr_pos_ts);
                    Hep3Vector ts_extr_pos_l_htf = sensor.getGeometry().getGlobalToLocal().transformed(ts_extr_pos);
                    
                    //System.out.println("pos:            " + ts_extr_pos_l.x() + "," + ts_extr_pos_l.y() + "," + ts_extr_pos_l.z());
                    //System.out.println("pos corr_htf:   " + ts_extr_pos_l_htf.x() + "," + ts_extr_pos_l_htf.y() + "," + ts_extr_pos_l_htf.z());
                    
                    
                    //System.out.println("pos tr :" + ts_extr_pos_trk.x() + "," + ts_extr_pos_trk.y() + "," + ts_extr_pos_trk.z());
                    //System.out.println("pos it :" + extr_pos_it.x() + "," + extr_pos_it.y() + "," + extr_pos_it.z());
                    //System.out.println("dir tr :" + dir_corr.x() + "," + dir_corr.y() + "," +dir_corr.z());
                    //System.out.println("s / s3D:" + s_corr + "/"+(s_corr / Math.cos(Math.atan(corr_htf.slope()))));
                    
                    gblStripData.setPath(s_corr);
                    gblStripData.setPath3D(s_corr / Math.cos(Math.atan(corr_htf.slope())));
                    gblStripData.setTrackDir(dir_corr);
                    gblStripData.setTrackPhi(corr_htf.phi0() - s_corr / corr_htf.R());
                    gblStripData.setTrackLambda(Math.atan(corr_htf.slope()));
                    gblStripData.setTrackPos(ts_extr_pos_l);
                    gblStripData.setTrackPosTrk(ts_extr_pos_trk);
                }
                //updated_gblStripData.add(gblStripData);
            }
            
            //update fit
            fit = HpsGblRefitter.fit(updated_gblStripData, bfac, false);
            //fit.addGBLStripClusterData(updated_gblStripData);
            //fit.addStripHits(sortedStripHits);
            newTrack = makeCorrectedTrack(fit, newTrack.getFirst(), allHthList, trackType, bfield);
            
            long endTime = System.nanoTime();
            long duration = (endTime - startTime);//divide by 1000000 to get milliseconds.
            aida.histogram1D("timing").fill((double)duration / 1000000.);
            //Fill the histograms
            FillIterationHistos(fit,iter);
        }
        
        Pair<Track, GBLKinkData> mergedTrack = makeCorrectedTrack(fit,TrackUtils.getHTF(newTrack.getFirst()),allHthList, trackType, bfield, storeTrackStates);
        return new Pair<Pair<Track, GBLKinkData>, FittedGblTrajectory>(mergedTrack, fit);
    }
    
    */
    public  Pair<Pair<Track, GBLKinkData>, FittedGblTrajectory> refitTrackWithTraj(HelicalTrackFit helix, Collection<TrackerHit> stripHits, Collection<TrackerHit> hth, int nIterations, int trackType, MultipleScattering scattering, double bfield, boolean storeTrackStates,boolean includeMS) {
        
        List<TrackerHit> allHthList = TrackUtils.sortHits(hth);
        List<TrackerHit> sortedStripHits = TrackUtils.sortHits(stripHits);
        FittedGblTrajectory fit = doGBLFit(helix, sortedStripHits, scattering, bfield, 0,includeMS);
        FillIterationHistos(fit,0);
        if(fit==null) return null;
        for (int i = 0; i < nIterations; i++) {
            long startTime = System.nanoTime();
            Pair<Track, GBLKinkData> newTrack = makeCorrectedTrack(fit, helix, allHthList, trackType, bfield);
            helix = TrackUtils.getHTF(newTrack.getFirst());
            fit = doGBLFit(helix, sortedStripHits, scattering, bfield, 0,includeMS);
            if (fit == null)
                return null;
            long endTime = System.nanoTime();
            long duration = (endTime - startTime);//divide by 1000000 to get milliseconds.
            FillIterationHistos(fit,i+1);
            aida.histogram1D("timing").fill((double)duration / 1000000.);
        }
        
        Pair<Track, GBLKinkData> mergedTrack = makeCorrectedTrack(fit, helix, allHthList, trackType, bfield, storeTrackStates);
        return new Pair<Pair<Track, GBLKinkData>, FittedGblTrajectory>(mergedTrack, fit);
    }
    
    
    /**
     * Do a GBL fit to a list of {@link TrackerHit}.
     *
     * @param htf - seed fit
     * @param stripHits - list of {@link TrackerHit}.
     * @param _scattering - estimation of the multiple scattering
     * {@link MultipleScattering}.
     * @param bfield - magnitude of B-field.
     * @param debug - debug flag.
     * @return the fitted GBL trajectory
     */
    public  FittedGblTrajectory doGBLFit(HelicalTrackFit htf, List<TrackerHit> stripHits, MultipleScattering _scattering, double bfield, int debug,boolean includeMS) {
        List<GBLStripClusterData> stripData = makeStripData(htf, stripHits, _scattering, bfield, debug, includeMS);
        if (stripData == null)
            return null;
        double bfac = Constants.fieldConversion * bfield;
        FittedGblTrajectory fit = HpsGblRefitter.fit(stripData, bfac, debug > 0);
        

        //Get the sensorArray in the fittedGblTrajector
        Integer [] sensorsFromMapArray = fit.getSensorMap().keySet().toArray(new Integer[0]);
        
        //Store the GBL strip cluster data to fitted traj
        //fit.addGBLStripClusterData(stripData);
        //Store the strip Hits used to build the stripData
        //fit.addStripHits(stripHits);
        
        
        /*
        //make Corrected track with all track states on surface
        //empty list 
        List<TrackerHit> hth_list = new ArrayList<TrackerHit>();
        Pair<Track, GBLKinkData> newTrack = makeCorrectedTrack(fit, htf,hth_list,0, bfield, true);
        
       
        System.out.println("Some checks on the fit results...");
        //Start from 1 to skip ID
        for (int i = 0  ; i<stripData.size(); i++) {
            
            GBLStripClusterData strip = stripData.get(i);
            int iLabel = sensorsFromMapArray[i];                    
            System.out.println("iLabel = " + iLabel + " i=" +i);
            //At each stripData, get the corrected track parameters
            Pair<double[],SymmetricMatrix> corrHelixParams_surface = fit.getCorrectedPerigeeParameters(htf,iLabel,bfield);
            double [] chisq = {0,0};
            int [] ndf   = {0,0};
            HelicalTrackFit corr_htf = new HelicalTrackFit(corrHelixParams_surface.getFirst(),
                                                           corrHelixParams_surface.getSecond(),
                                                           chisq,ndf,null,null);
            
            //Get the measurements 
            if (stripData.get(i).getScatterOnly()==0) { 
                                
                SiTrackerHitStrip1D h = null;
                //Could cache this too
                if (stripHits.get(i) instanceof SiTrackerHitStrip1D) {
                    h = (SiTrackerHitStrip1D) stripHits.get(i);
                }
                else {
                    h = new SiTrackerHitStrip1D(stripHits.get(i));
                }
                

                SiTrackerHitStrip1D local  = h.getTransformedHit(TrackerHitType.CoordinateSystem.SENSOR);
                SiTrackerHitStrip1D global = h.getTransformedHit(TrackerHitType.CoordinateSystem.GLOBAL);
                
                ITransform3D trans = local.getLocalToGlobal();
                

                //TrkToDet
                Hep3Matrix trk2det = CoordinateTransformations.getMatrixInverse();
                
                
                System.out.println("Original prediction: " + 
                                   stripData.get(i).getTrackPos().x() + "," +
                                   stripData.get(i).getTrackPos().y() + "," +
                                   stripData.get(i).getTrackPos().z());
                System.out.println("Original 3D path: " + stripData.get(i).getPath3D());
                System.out.println("Original path: " + stripData.get(i).getPath());
                
                //These are from the original points.
                double s_original = HelixUtils.PathToXPlane(htf, stripData.get(i).getTrackPosTrk().x(),0.,0).get(0);
                Hep3Vector dir = HelixUtils.Direction(htf,s_original);

                System.out.println("Position in Tracking ref: " + stripData.get(i).getTrackPosTrk().x() +","+stripData.get(i).getTrackPosTrk().y()+","+stripData.get(i).getTrackPosTrk().z());

                Hep3Vector posGlobalRef = VecOp.mult(trk2det,stripData.get(i).getTrackPosTrk());

                System.out.println("Position in Global ref: " + posGlobalRef.x() +","+posGlobalRef.y()+","+posGlobalRef.z());
                
                Hep3Vector posLocalRef = trans.inverse().transformed(posGlobalRef);
                System.out.println("Position in Local ref 1: " + posLocalRef.x() +","+posLocalRef.y()+","+posLocalRef.z());
                
                
                HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit) h.getRawHits().get(0)).getDetectorElement();
                Hep3Vector posLocalRef_2 = sensor.getGeometry().getGlobalToLocal().transformed(posGlobalRef);
                System.out.println("Position in Local ref 2: " + posLocalRef_2.x() +","+posLocalRef_2.y()+","+posLocalRef_2.z());

                System.out.println("Recomputed path to position: " + s_original);

                //These are from the corrected helix - without changing the location
                double s_corr = HelixUtils.PathToXPlane(corr_htf, stripData.get(i).getTrackPosTrk().x(),0.,0).get(0);
                Hep3Vector dir_corr = HelixUtils.Direction(corr_htf,s_corr);
                System.out.println("Recomputed path to position with corrected track: " + s_corr);
                
                //Now update the location and check the paths
                GblData gbldata = fit.get_traj().theData.get(i);
                double[] results = new double[5];
                List<Integer> indLocal = new ArrayList<Integer>();
                List<Double> derLocal = new ArrayList<Double>();
                gbldata.getResidual(results,indLocal,derLocal);
                double xprime = (strip.getTrackPos().x() + results[4]);
                double yprime = strip.getTrackPos().y();
                double zprime = strip.getTrackPos().z();
                System.out.println("New pred = (" + xprime + "," +yprime + "," + zprime+")");
                
                //New prediction in local coordinates
                Hep3Vector new_pred_l = new BasicHep3Vector(xprime, yprime, zprime);
                
                //New prediction in global coordinates
                Hep3Vector new_pred_g = trans.transformed(new_pred_l);
                
                //New prediction in tracking coordinates 
                Hep3Vector new_pred_trk = CoordinateTransformations.transformVectorToTracking(new_pred_g);

                
                //Check the position on the sensor made from the extrapolated corrected track
                //In global to the sensor
                Hep3Vector global_extr_pos       = TrackStateUtils.getLocationAtSensor(corr_htf,sensor,bfield);

                System.out.println("New Position in global ref: " + new_pred_g.x() +","+new_pred_g.y()+","+new_pred_g.z());
                System.out.println("New Position in global ref (track extr): " + global_extr_pos.x() +","+global_extr_pos.y()+","+global_extr_pos.z());

                System.out.println("New Position in trk ref: " + new_pred_trk.x() +","+new_pred_trk.y()+","+new_pred_trk.z());
                
                //These are from the corrected helix - without changing the location
                double s_corr_2 = HelixUtils.PathToXPlane(corr_htf, new_pred_trk.x(),0.,0).get(0);
                Hep3Vector dir_corr_2 = HelixUtils.Direction(corr_htf,s_corr_2);
                System.out.println("Recomputed path to corrected position with corrected track: " + s_corr_2);
                
                //These recomputed positions come from the track state on surfaces
                Hep3Vector ts_extr_pos = TrackStateUtils.getLocationAtSensor(newTrack.getFirst(),sensor,bfield);
                if (ts_extr_pos != null) {
                    System.out.println("New Position in global ref (ts): " + ts_extr_pos.x() +","+ts_extr_pos.y()+","+ts_extr_pos.z());
                }
            }
        }
        */
        return fit;
    }
        
    /**
     * Create a list of {@link GBLStripClusterData} objects that can be used as
     * input to the GBL fitter.
     *
     * @param htf
     * @param stripHits
     * @param _scattering
     * @param _B
     * @param _debug
     * @return the list of GBL strip cluster data
     */
    public List<GBLStripClusterData> makeStripData(HelicalTrackFit htf, List<TrackerHit> stripHits, MultipleScattering _scattering, double _B, int _debug,boolean includeMS) {
        List<GBLStripClusterData> stripClusterDataList = new ArrayList<GBLStripClusterData>();

        // Find scatter points along the path
        //In principle I could use this to add the hits - TODO ?
        MultipleScattering.ScatterPoints scatters = _scattering.FindHPSScatterPoints(htf);
        
        //Loop over the Scatters
        
        //Two things can happen here at each iteration:
        //1) Nscatters >= Nhits on tracks =>
        //   Build GBLStripClusterData for both measurements and scatters
        //2) Nscatters < Nhits on track =>
        //   The hit has been associated to the track but couldn't find the scatter on the volume => create a scatter at that sensor
        //   Does this makes sense?

        //Case (1)
        //Check if the scatter has a measurement => then use the usual way to build the HelicalTrackStripGbl
        
        int nhits = 0;
        int nscatters = 0;
        boolean addScatters = true;
            
        if (scatters.getPoints().size() >= stripHits.size() && includeMS) {
            
            for (MultipleScattering.ScatterPoint scatter : scatters.getPoints()) {
            
                boolean MeasOnScatter = false;
            
                //This is bit inefficient as it always loop on all the stripHits
                //TODO: optimize and only check hit once if is in the scatters.
                
                HpsSiSensor scatter_sensor = (HpsSiSensor) scatter.getDet();
            
                for (TrackerHit stripHit : stripHits) {
                
                    if (MeasOnScatter)
                        continue;
                
                    IDetectorElement det_element = ((RawTrackerHit) stripHit.getRawHits().get(0)).getDetectorElement();
                
                    if (det_element.equals(scatter.getDet())) {
                        MeasOnScatter = true;
                        HpsSiSensor hit_sensor = (HpsSiSensor) det_element;
                        nhits+=1;
                        HelicalTrackStripGbl gbl_strip;
                        if (stripHit instanceof SiTrackerHitStrip1D) {
                            gbl_strip = new HelicalTrackStripGbl(makeDigiStrip((SiTrackerHitStrip1D) stripHit),true);
                        } else {
                            SiTrackerHitStrip1D newHit = new SiTrackerHitStrip1D(stripHit);
                            gbl_strip = new HelicalTrackStripGbl(makeDigiStrip(newHit), true);
                        }
                        //hit_sensor or scatter_sensor, they are the same here
                        GBLStripClusterData stripData = makeStripData(hit_sensor, gbl_strip, htf, scatter);
                        if (stripData != null) {
                            stripClusterDataList.add(stripData); 
                        }
                        else {
                            System.out.printf("WARNING::MakeGblTracks::Couldn't make stripClusterData for hps sensor. Skipping scatter");
                        }
                    }
                }
            
                if (!MeasOnScatter && addScatters) {
                    //No measurement
                    nscatters+=1;
                    //If the scatter has no measurement => then only build a scatter point 
                    //Here only scatter sensor is available
                    GBLStripClusterData stripData = makeScatterOnlyData(scatter_sensor,htf,scatter);
                    
                    if (stripData != null) {
                        stripClusterDataList.add(stripData);
                    }
                }
            } // loop on scatters
        } // more scatters than hits
        
        else { //more hits than scatters 
            
   
            for (TrackerHit stripHit : stripHits) {
                HelicalTrackStripGbl strip;
                if (stripHit instanceof SiTrackerHitStrip1D) {
                    strip = new HelicalTrackStripGbl(makeDigiStrip((SiTrackerHitStrip1D) stripHit), true);
                } else {
                    SiTrackerHitStrip1D newHit = new SiTrackerHitStrip1D(stripHit);
                    strip = new HelicalTrackStripGbl(makeDigiStrip(newHit),true);
                }
                HpsSiSensor sensor = (HpsSiSensor) ((RawTrackerHit) stripHit.getRawHits().get(0)).getDetectorElement();
                MultipleScattering.ScatterPoint temp = scatters.getScatterPoint(((RawTrackerHit) strip.getStrip().rawhits().get(0)).getDetectorElement());
                
                //This is done to correct the fact that the helical fit might not hit a volume
                //But hit is associated to the track => must scatter
                if (temp == null){
                    temp = getScatterPointGbl(sensor, strip, htf, _scattering, _B);
                    if(temp == null){
                        return null;
                    }
                }
                GBLStripClusterData stripData = makeStripData(sensor, strip, htf, temp);
                if (stripData != null)
                    stripClusterDataList.add(stripData);
            }
        }
        
        return stripClusterDataList;
    }

    public  MultipleScattering.ScatterPoint getScatterPointGbl(HpsSiSensor sensor, HelicalTrackStripGbl strip, HelicalTrackFit htf, MultipleScattering _scattering, double _B) {

        MultipleScattering.ScatterPoint temp = null;

        Hep3Vector pos = TrackUtils.getHelixPlaneIntercept(htf, strip, Math.abs(_B));
        if (pos == null) {
            System.out.println("Can't find track intercept; aborting Track refit");
            return null;
        }
        ScatterAngle scatAngle = new ScatterAngle((HelixUtils.PathToXPlane(htf, pos.x(), 0, 0).get(0)), GblUtils.estimateScatter(sensor, htf, _scattering, _B));
        temp = new MultipleScattering.ScatterPoint(((RawTrackerHit) strip.getStrip().rawhits().get(0)).getDetectorElement(), scatAngle);
        temp.setPosition(pos);
        temp.setDirection(HelixUtils.Direction(htf, scatAngle.PathLen()));

        return temp;
    }
    
    public  GBLStripClusterData makeScatterOnlyData(HpsSiSensor sensor, HelicalTrackFit htf, MultipleScattering.ScatterPoint temp) {
        
        if (temp ==null)
            return null;
        
        // find Millepede layer definition from DetectorElement
        int millepedeId = sensor.getMillepedeId();
        // find volume of the sensor (top or bottom)
        int volume = sensor.isTopLayer() ? 0 : 1;
                                                
        // GBLDATA
        GBLStripClusterData stripData = new GBLStripClusterData(millepedeId);
        stripData.setVolume(volume);
        
        //This GBLStripData doesn't hold a measurement
        stripData.setScatterOnly(1);

        double s3D = temp.getScatterAngle().PathLen() / Math.cos(Math.atan(htf.slope()));
        stripData.setPath(temp.getScatterAngle().PathLen());
        stripData.setPath3D(s3D);
        
        //Do not set U,V,W
        
        // Print track direction at intercept
        double phi = htf.phi0() - temp.getScatterAngle().PathLen() / htf.R();
        double lambda = Math.atan(htf.slope());
        
        stripData.setTrackDir(temp.getDirection());
        stripData.setTrackPhi(phi);
        stripData.setTrackLambda(lambda);
        
        //Do not set the measurement. 
        //Set a negative large error
        stripData.setMeasErr(-9999);

        stripData.setScatterAngle(temp.getScatterAngle().Angle());
        
        return stripData;
    }

    public  GBLStripClusterData makeStripData(HpsSiSensor sensor, HelicalTrackStripGbl strip, HelicalTrackFit htf, MultipleScattering.ScatterPoint temp) {
        if (temp == null)
            return null;

        // find Millepede layer definition from DetectorElement
        int millepedeId = sensor.getMillepedeId();
        // find volume of the sensor (top or bottom)
        int volume = sensor.isTopLayer() ? 0 : 1;

        // Center of the sensor
        Hep3Vector origin = strip.origin();

        // GBLDATA
        GBLStripClusterData stripData = new GBLStripClusterData(millepedeId);

        // Add the volume
        stripData.setVolume(volume);

        //This GBLStripData holds a measurement
        stripData.setScatterOnly(0);
        
        //Set Origin of the sensor (in tracking frame) of this gbl strip cluster data
        stripData.setOrigin(origin);

        // Add to output list

        // path length to intercept
        double s3D = temp.getScatterAngle().PathLen() / Math.cos(Math.atan(htf.slope()));

        // GBLDATA
        stripData.setPath(temp.getScatterAngle().PathLen());
        stripData.setPath3D(s3D);

        // GBLDATA
        stripData.setU(strip.u());
        stripData.setV(strip.v());
        stripData.setW(strip.w());

        // Print track direction at intercept
        double phi = htf.phi0() - temp.getScatterAngle().PathLen() / htf.R();
        double lambda = Math.atan(htf.slope());

        // GBLDATA
        stripData.setTrackDir(temp.getDirection());
        stripData.setTrackPhi(phi);
        stripData.setTrackLambda(lambda);

        // Print residual in measurement system
        // start by find the distance vector between the center and the track position
        Hep3Vector vdiffTrk = VecOp.sub(temp.getPosition(), origin);

        // then find the rotation from tracking to measurement frame
        Hep3Matrix trkToStripRot = getTrackToStripRotation(sensor);

        // then rotate that vector into the measurement frame to get the predicted measurement position
        Hep3Vector trkpos_meas = VecOp.mult(trkToStripRot, vdiffTrk);

        // GBLDATA
        stripData.setMeas(strip.umeas());
        stripData.setTrackPos(trkpos_meas);
        stripData.setTrackPosTrk(temp.getPosition());
        stripData.setMeasErr(strip.du());

        // GBLDATA
        stripData.setScatterAngle(temp.getScatterAngle().Angle());

        return stripData;
    }

    private  Hep3Matrix getTrackToStripRotation(SiSensor sensor) {
        // This function transforms the hit to the sensor coordinates

        // Transform from JLab frame to sensor frame (done through the RawTrackerHit)
        SiSensorElectrodes electrodes = sensor.getReadoutElectrodes(ChargeCarrier.HOLE);
        ITransform3D detToStrip = electrodes.getGlobalToLocal();
        // Get rotation matrix
        Hep3Matrix detToStripMatrix = detToStrip.getRotation().getRotationMatrix();

        return VecOp.mult(detToStripMatrix, CoordinateTransformations.getMatrixInverse());
    }

    private  HelicalTrackStrip makeDigiStrip(SiTrackerHitStrip1D h) {
        SiTrackerHitStrip1D local = h.getTransformedHit(TrackerHitType.CoordinateSystem.SENSOR);
        SiTrackerHitStrip1D global = h.getTransformedHit(TrackerHitType.CoordinateSystem.GLOBAL);

        ITransform3D trans = local.getLocalToGlobal();
        Hep3Vector org = trans.transformed(new BasicHep3Vector(0., 0., 0.));
        Hep3Vector u = global.getMeasuredCoordinate();
        Hep3Vector v = global.getUnmeasuredCoordinate();

        // rotate to tracking frame
        Hep3Vector neworigin = CoordinateTransformations.transformVectorToTracking(org);
        Hep3Vector newu = CoordinateTransformations.transformVectorToTracking(u);
        Hep3Vector newv = CoordinateTransformations.transformVectorToTracking(v);

        double umeas = local.getPosition()[0];
        double vmin = VecOp.dot(local.getUnmeasuredCoordinate(), local.getHitSegment().getStartPoint());
        double vmax = VecOp.dot(local.getUnmeasuredCoordinate(), local.getHitSegment().getEndPoint());
        double du = Math.sqrt(local.getCovarianceAsMatrix().diagonal(0));
        double dEdx = h.getdEdx();
        double time = h.getTime();
        List<RawTrackerHit> rawhits = h.getRawHits();
        HelicalTrackStrip strip = new HelicalTrackStrip(neworigin, newu, newv, umeas, du, vmin, vmax, dEdx, time, rawhits, null, -1, null);

        return strip;
    }
}
