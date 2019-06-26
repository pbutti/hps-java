package org.hps.analysis.MC;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.conditions.hodoscope.HodoscopeChannel;
import org.hps.conditions.hodoscope.HodoscopeChannel.HodoscopeChannelCollection;
import org.hps.detector.hodoscope.HodoscopeDetectorElement;
import org.hps.readout.triggerstudies.PlotToTextModule;
import org.hps.record.triggerbank.TriggerModule;
import org.hps.util.Pair;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.MCParticle;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.FieldMap;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.lcsim.util.swim.VectorArithmetic;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;

public class TriggerTuningDriver extends Driver {
    private int totalEvents = 0;
    private int goodEvents = 0;
    private int singlesEvents = 0;
    private int pairEvents = 0;
    private int tripletEvents = 0;
    private int coptEvents = 0;
    private int hodoscopeEvents = 0;
    
    private double coptXLowerBound = 90.0;
    private double chiSquaredUpperBound = Double.MAX_VALUE;
    private double hodoscopeHitEnergyLowerBound = 1.3 / 1000.0;
    
    private String gblTrackCollectionName = "GBLTracks";
    private String gtpClusterCollectionName = "EcalClustersGTP";
    private String hodoscopeHitCollectionName = "HodoscopeCorrectedHits";
    
    private double[] boundsPositronP = new double[2];
    private double[] boundsElectronP = new double[2];
    private double[][] boundsTopPositron = new double[2][0];
    private double[][] boundsBotPositron = new double[2][0];
    private double[][] boundsTopElectron = new double[2][0];
    private double[][] boundsBotElectron = new double[2][0];
    
    private String outputDirectory = ".";
    private FieldMap fieldMap = null;
    private HodoscopeDetectorElement hodoscopeDetectorElement;
    private Map<Long, HodoscopeChannel> hodoscopeChannelMap = new HashMap<Long, HodoscopeChannel>();
    
    private static final String CHI_SQUARED = "Cluster\\Track-Matching/Chi Squared Distribution";
    private static final String DEBUG_TRACK_POSITION = "Debug/All Track Position Distribution";
    private static final String DEBUG_CLUSTER_POSITION = "Debug/All Cluster Position Distribution";
    
    private static final String HODOSCOPE_SCINTILLATOR_ENERGY = "Hodoscope/Scintillator Energy Distribution";
    
    private static final String SINGLES_POSITION = "Singles Trigger Cuts/Cluster Position Distribution";
    private static final String SINGLES_SEED_ENERGY = "Singles Trigger Cuts/Cluster Seed Energy Distribution";
    private static final String SINGLES_MULTIPLICITY = "Singles Trigger Cuts/Cluster Multiplicity Distribution";
    private static final String SINGLES_TOTAL_ENERGY = "Singles Trigger Cuts/Cluster Total Energy Distribution";
    
    private static final String PAIR_POSITION = "Pair Trigger Cuts/Pair Cluster Position Distribution";
    private static final String PAIR_MULTIPLICITY = "Pair Trigger Cuts/Pair Cluster Multiplicity Distribution";
    private static final String PAIR_TOTAL_ENERGY = "Pair Trigger Cuts/Pair Cluster Total Energy Distribution";
    private static final String PAIR_ENERGY_SUM = "Pair Trigger Cuts/Pair Energy Sum Distribution";
    private static final String PAIR_ENERGY_DIFF = "Pair Trigger Cuts/Pair Energy Difference Distribution";
    private static final String PAIR_ENERGY_SLOPE = "Pair Trigger Cuts/Pair Energy Slope Distribution";
    private static final String PAIR_COPLANARITY = "Pair Trigger Cuts/Pair Coplanarity Distribution";
    
    private static final String COPT_POSITION = "COP Trigger Cuts/Cluster Position Distribution";
    private static final String COPT_SEED_ENERGY = "COP Trigger Cuts/Cluster Seed Energy Distribution";
    private static final String COPT_MULTIPLICITY = "COP Trigger Cuts/Cluster Multiplicity Distribution";
    private static final String COPT_TOTAL_ENERGY = "COP Trigger Cuts/Cluster Total Energy Distribution";
    private static final String COPT_ENERGY_POSITION = "COP Trigger Cuts/Cluster Total Energy vs. Cluster Position Distribution";
    private static final String DEBUG_MOMENTUM_POSITION = "Debug/Track Momentum vs. Cluster Position Distribution";
    private static final String DEBUG_PHI_POSITION = "Debug/Track Phi vs. Cluster Position Distribution";
    
    private static final String INV_MASS_TRUTH = "Invariant Mass/Invariant Mass Distribution (Truth)";
    private static final String INV_MASS_REALLY_NO_CUTS = "Invariant Mass/Invariant Mass Distribution (Really No Cuts)";
    private static final String INV_MASS_NO_CUTS = "Invariant Mass/Invariant Mass Distribution (No Cuts)";
    private static final String INV_MASS_95 = "Invariant Mass/Invariant Mass Distribution (95% COPT Threshold)";
    private static final String INV_MASS_97 = "Invariant Mass/Invariant Mass Distribution (97% COPT Threshold)";
    private static final String INV_MASS_ME = "Invariant Mass/Invariant Mass Distribution (99% COPT Threshold)";
    
    private static final String HODOSCOPE_TRUTH_ENERGY = "Hodoscope/Truth Hit Energy Distribution";
    private static final String HODOSCOPE_TRUTH_COMP_ENERGY = "Hodoscope/Truth Hit Energy Distribution (Compiled)";
    
    @Override
    public void detectorChanged(Detector detector) {
        // Get the hodoscope channels.
        final DatabaseConditionsManager conditions = DatabaseConditionsManager.getInstance();
        final HodoscopeChannelCollection channels = conditions.getCachedConditions(HodoscopeChannelCollection.class, "hodo_channels").getCachedData();
        
        // Map them to their IDs.
        for(HodoscopeChannel channel : channels) {
            hodoscopeChannelMap.put(Long.valueOf(channel.getChannelId().longValue()), channel);
        }
        
        // Update the hodoscope detector object.
        hodoscopeDetectorElement = (HodoscopeDetectorElement) detector.getSubdetector("Hodoscope").getDetectorElement();
        
        // Update the fieldmap.
        fieldMap = detector.getFieldMap();
    }
    
    @Override
    public void endOfData() {
        // Print the number of each type of event that were seen.
        System.out.println("Total Events Processed: " + totalEvents);
        System.out.println("\tAnalyzable Events :: " + goodEvents);
        System.out.println("\tSingles Events    :: " + singlesEvents);
        System.out.println("\tPair Events       :: " + pairEvents);
        System.out.println("\tTriplet Events    :: " + tripletEvents);
        System.out.println("\tCOPT Events       :: " + coptEvents);
        System.out.println("\tHodoscope Events  :: " + hodoscopeEvents);
        
        // Output the Mathematica dat files.
        try {
            // Preliminary plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(CHI_SQUARED), new File(outputDirectory + File.separator + "tuning_chiSquared.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(HODOSCOPE_SCINTILLATOR_ENERGY), new File(outputDirectory + File.separator + "tuning_hodoHitEnergy.dat"));
            
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(DEBUG_MOMENTUM_POSITION),
                    new File(outputDirectory + File.separator + "debug_copt_momentumVsPosition.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(DEBUG_PHI_POSITION),
                    new File(outputDirectory + File.separator + "debug_copt_phiVsPosition.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(DEBUG_TRACK_POSITION),
                    new File(outputDirectory + File.separator + "debug_track_position.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(DEBUG_CLUSTER_POSITION),
                    new File(outputDirectory + File.separator + "debug_cluster_position.dat"));
            
            // Cluster/track matching plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(true, true, false)),
                        new File(outputDirectory + File.separator + "tuning_clusterTrackMatching_TPR.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(true, false, false)),
                    new File(outputDirectory + File.separator + "tuning_clusterTrackMatching_TNR.dat"));
        
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(false, true, false)),
                    new File(outputDirectory + File.separator + "tuning_clusterTrackMatching_BPR.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(false, false, false)),
                    new File(outputDirectory + File.separator + "tuning_clusterTrackMatching_BNR.dat"));
            
            // Cluster/track matching debug recon plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(true, true, true)),
                        new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_TPR.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(true, false, true)),
                    new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_TNR.dat"));
        
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(false, true, true)),
                    new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_BPR.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(false, false, true)),
                    new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_BNR.dat"));
            
            
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingXPlotName(true, true, true)),
                        new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_TPX.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingXPlotName(true, false, true)),
                    new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_TNX.dat"));
        
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingXPlotName(false, true, true)),
                    new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_BPX.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingXPlotName(false, false, true)),
                    new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_BNX.dat"));
            
            
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingYPlotName(true, true, true)),
                        new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_TPY.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingYPlotName(true, false, true)),
                    new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_TNY.dat"));
        
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingYPlotName(false, true, true)),
                    new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_BPY.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingYPlotName(false, false, true)),
                    new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_BNY.dat"));
            
            // Singles trigger plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(SINGLES_POSITION),
                    new File(outputDirectory + File.separator + "tuning_singles_position.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(SINGLES_MULTIPLICITY),
                    new File(outputDirectory + File.separator + "tuning_singles_multiplicity.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(SINGLES_SEED_ENERGY),
                    new File(outputDirectory + File.separator + "tuning_singles_seedEnergy.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(SINGLES_TOTAL_ENERGY),
                    new File(outputDirectory + File.separator + "tuning_singles_totalEnergy.dat"));
            
            // Pair trigger plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(PAIR_POSITION),
                    new File(outputDirectory + File.separator + "tuning_pair_position.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(PAIR_MULTIPLICITY),
                    new File(outputDirectory + File.separator + "tuning_pair_multiplicity.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(PAIR_TOTAL_ENERGY),
                    new File(outputDirectory + File.separator + "tuning_pair_totalEnergy.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(PAIR_ENERGY_SUM),
                    new File(outputDirectory + File.separator + "tuning_pair_energySum.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(PAIR_ENERGY_DIFF),
                    new File(outputDirectory + File.separator + "tuning_pair_energyDifference.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(PAIR_ENERGY_SLOPE),
                    new File(outputDirectory + File.separator + "tuning_pair_energySlope.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(PAIR_COPLANARITY),
                    new File(outputDirectory + File.separator + "tuning_pair_coplanarity.dat"));
            
            // Hodoscope layer-to-layer plots.
            for(int ix = 0; ix < 9; ix++) {
                PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(getLayerToLayerPlotName(ix, true)),
                        new File(outputDirectory + File.separator + "tuning_hodoLayerToLayer_T" + (ix + 1) + ".dat"));
                PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(getLayerToLayerPlotName(ix, false)),
                        new File(outputDirectory + File.separator + "tuning_hodoLayerToLayer_B" + (ix + 1) + ".dat"));
                
                PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(getLayerToCalorimeterPlotName(ix, true, 0)),
                        new File(outputDirectory + File.separator + "tuning_hodoLayerToCalorimeter_L1_T" + (ix + 1) + ".dat"));
                PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(getLayerToCalorimeterPlotName(ix, false, 0)),
                        new File(outputDirectory + File.separator + "tuning_hodoLayerToCalorimeter_L1_B" + (ix + 1) + ".dat"));
                PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(getLayerToCalorimeterPlotName(ix, true, 1)),
                        new File(outputDirectory + File.separator + "tuning_hodoLayerToCalorimeter_L2_T" + (ix + 1) + ".dat"));
                PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(getLayerToCalorimeterPlotName(ix, false, 1)),
                        new File(outputDirectory + File.separator + "tuning_hodoLayerToCalorimeter_L2_B" + (ix + 1) + ".dat"));
            }
            
            // COP-trigger plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(COPT_POSITION),
                    new File(outputDirectory + File.separator + "tuning_copt_position.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(COPT_ENERGY_POSITION),
                    new File(outputDirectory + File.separator + "tuning_copt_energyVsPosition.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(COPT_MULTIPLICITY),
                    new File(outputDirectory + File.separator + "tuning_copt_multiplicity.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(COPT_SEED_ENERGY),
                    new File(outputDirectory + File.separator + "tuning_copt_seedEnergy.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(COPT_TOTAL_ENERGY),
                    new File(outputDirectory + File.separator + "tuning_copt_totalEnergy.dat"));
            
            // Invariant mass plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(INV_MASS_TRUTH),
                    new File(outputDirectory + File.separator + "tuning_invaritanMass_truth.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(INV_MASS_REALLY_NO_CUTS),
                    new File(outputDirectory + File.separator + "tuning_invaritanMass_reallyNoCuts.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(INV_MASS_NO_CUTS),
                    new File(outputDirectory + File.separator + "tuning_invaritanMass_noCuts.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(INV_MASS_95),
                    new File(outputDirectory + File.separator + "tuning_invaritanMass_95.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(INV_MASS_97),
                    new File(outputDirectory + File.separator + "tuning_invaritanMass_97.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(INV_MASS_ME),
                    new File(outputDirectory + File.separator + "tuning_invaritanMass_99.dat"));
            
            // Write the hodoscope truth hit debug plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(HODOSCOPE_TRUTH_ENERGY), new File(outputDirectory + File.separator + "debug_hodoTruthEnergy.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(HODOSCOPE_TRUTH_COMP_ENERGY), new File(outputDirectory + File.separator + "debug_hodoTruthEnergy_Compiled.dat"));
        } catch(IOException e) {
            throw new RuntimeException(e);
        }
    }
    
    @Override
    public void process(EventHeader event) {
        // Track the total number of processed events.
        totalEvents++;
        
        // Populate the chi squared plots. This distribution is used
        // to determine the chi squared cut for "good events," so it
        // is populated before the "good event" cut.
        fillChiSquaredDistro(event, gblTrackCollectionName);
        
        // Populate the debugging cluster position plot.
        fillClusterPositionDistro(event, gtpClusterCollectionName);
        fillTrackPositionDistro(event, gblTrackCollectionName);
        
        // Perform the cluster/track matching analysis. This does not
        // require that an event be analyzable to be useful.
        performTrackAnalysis(event);
        
        // Get cluster/track matched pairs.
        List<Pair<Cluster, Track>> pairList = TriggerTuningUtilityModule.getClusterTrackMatchedPairs(
                TriggerTuningUtilityModule.getCollection(event, gtpClusterCollectionName, Cluster.class),
                TriggerTuningUtilityModule.getCollection(event, gblTrackCollectionName, Track.class),
                fieldMap, boundsTopPositron, boundsBotPositron, boundsTopElectron, boundsBotElectron,
                boundsPositronP, boundsElectronP);
        
        // Perform the absolutely no cuts invariant mass analysis.
        performInvariantMassAnalysisTruth(TriggerTuningUtilityModule.getCollection(event, "MCParticle", MCParticle.class));
        performInvariantMassAnalysisNoClusters(TriggerTuningUtilityModule.getCollection(event, gblTrackCollectionName, Track.class));
        //performInvariantMassAnalysisMatched(pairList);
        
        // Check if this is a good event. If it isn't do nothing.
        if(!TriggerTuningUtilityModule.isGoodEvent(event, gblTrackCollectionName, chiSquaredUpperBound)) {
            return;
        }
        goodEvents++;
        
        // Track whether this meets the additional criteria for each
        // of the triggers. The singles trigger criteria is always
        // met by default.
        boolean meetsSinglesCriteria = true;
        boolean meetsPairCriteria = TriggerTuningUtilityModule.isPairEvent(event, gtpClusterCollectionName);
        boolean meetsTripletCriteria = TriggerTuningUtilityModule.isTripletEvent(event, gtpClusterCollectionName);
        boolean meetsCOPTCriteria = TriggerTuningUtilityModule.isCOPTEvent(event, gtpClusterCollectionName, coptXLowerBound);
        boolean meetsHodoscopeCriteria = TriggerTuningUtilityModule.isHodoscopeEvent(event, hodoscopeHitCollectionName, gtpClusterCollectionName, coptXLowerBound);
        
        // Perform singles trigger tuning.
        if(meetsSinglesCriteria) {
            singlesEvents++;
            performSinglesTriggerAnalysis(pairList);
        }
        
        // Perform pair trigger tuning.
        if(meetsPairCriteria) {
            pairEvents++;
            performPairTriggerAnalysis(pairList);
        }
        
        // Perform triplet trigger tuning.
        if(meetsTripletCriteria) {
            tripletEvents++;
        }
        
        // Perform COPT trigger tuning.
        if(meetsCOPTCriteria) {
            coptEvents++;
            performCOPTriggerAnalysis(pairList);
        }
        
        // Perform hodoscope trigger tuning.
        if(meetsHodoscopeCriteria) {
            hodoscopeEvents++;
            performHodoscopeAnalysis(event);
        }
        
        performInvariantMassAnalysis(pairList);
    }
    
    @Override
    public void startOfData() {
        // Preliminary plots.
        AIDA.defaultInstance().histogram1D(CHI_SQUARED, 50, 0.0, 50.0);
        AIDA.defaultInstance().histogram1D(HODOSCOPE_SCINTILLATOR_ENERGY, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram2D(DEBUG_CLUSTER_POSITION, 350, -350, 350, 100, -100, 100);
        AIDA.defaultInstance().histogram2D(DEBUG_PHI_POSITION, 23, 0.5, 23.5, 100, -100, 100);
        AIDA.defaultInstance().histogram2D(DEBUG_MOMENTUM_POSITION, 23, 0.5, 23.5, 250, 0.000, 5.0);
        AIDA.defaultInstance().histogram2D(DEBUG_TRACK_POSITION, 350, -350, 350, 100, -100, 100);
        
        // Cluster/Track matching plots.
        for(int i = 0; i < 2; i++) {
            boolean isTop = (i == 0);
            for(int j = 0; j < 2; j++) {
                boolean isPositive = (j == 0);
                for(int k = 0; k < 2; k++) {
                    boolean isRecon = (k == 1);
                    AIDA.defaultInstance().histogram2D(getClusterTrackMatchingXPlotName(isTop, isPositive, isRecon), 250, 0.000, 5.000, 70, 0, 70);
                    AIDA.defaultInstance().histogram2D(getClusterTrackMatchingYPlotName(isTop, isPositive, isRecon), 250, 0.000, 5.000, 70, 0, 70);
                    AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(isTop, isPositive, isRecon), 250, 0.000, 5.000, 70, 0, 70);
                }
            }
        }
        
        // Hodoscope Layer-to-Layer correlation plots.
        for(int clusterIndex = 0; clusterIndex < 9; clusterIndex++) {
            AIDA.defaultInstance().histogram1D(getLayerToLayerPlotName(clusterIndex, true), 9, 0.5, 9.5);
            AIDA.defaultInstance().histogram1D(getLayerToLayerPlotName(clusterIndex, false), 9, 0.5, 9.5);
            
            AIDA.defaultInstance().histogram1D(getLayerToCalorimeterPlotName(clusterIndex, true, 0), 23, 0.5, 23.5);
            AIDA.defaultInstance().histogram1D(getLayerToCalorimeterPlotName(clusterIndex, false, 0), 23, 0.5, 23.5);
            AIDA.defaultInstance().histogram1D(getLayerToCalorimeterPlotName(clusterIndex, true, 1), 23, 0.5, 23.5);
            AIDA.defaultInstance().histogram1D(getLayerToCalorimeterPlotName(clusterIndex, false, 1), 23, 0.5, 23.5);
        }
        
        // Singles trigger plots.
        AIDA.defaultInstance().histogram1D(SINGLES_MULTIPLICITY, 9, 0.5, 9.5);
        AIDA.defaultInstance().histogram1D(SINGLES_SEED_ENERGY, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram1D(SINGLES_TOTAL_ENERGY, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram2D(SINGLES_POSITION, 47, -23.5, 23.5, 11, -5.5, 5.5);
        
        // Pair trigger plots.
        AIDA.defaultInstance().histogram1D(PAIR_MULTIPLICITY, 9, 0.5, 9.5);
        AIDA.defaultInstance().histogram1D(PAIR_TOTAL_ENERGY, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram2D(PAIR_POSITION, 47, -23.5, 23.5, 11, -5.5, 5.5);
        AIDA.defaultInstance().histogram1D(PAIR_ENERGY_SUM, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram1D(PAIR_ENERGY_DIFF, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram2D(PAIR_ENERGY_SLOPE, 200, -400, 400, 250, 0.000, 5.000);
        AIDA.defaultInstance().histogram1D(PAIR_COPLANARITY, 90, 0.000, 180);
        
        // COP-trigger plots.
        AIDA.defaultInstance().histogram1D(COPT_MULTIPLICITY, 9, 0.5, 9.5);
        AIDA.defaultInstance().histogram1D(COPT_SEED_ENERGY, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram1D(COPT_TOTAL_ENERGY, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram2D(COPT_POSITION, 47, -23.5, 23.5, 11, -5.5, 5.5);
        AIDA.defaultInstance().histogram2D(COPT_ENERGY_POSITION, 23, 0.5, 23.5, 250, 0.000, 5.0);
        
        // Hodoscope Hit Truth Comparison Debug Plots
        AIDA.defaultInstance().histogram1D(HODOSCOPE_TRUTH_ENERGY, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram1D(HODOSCOPE_TRUTH_COMP_ENERGY, 500, 0.000, 5.000);
        
        // Invariant mass plots.
        AIDA.defaultInstance().histogram1D(INV_MASS_TRUTH, 750, 0.0000, 0.3000);
        AIDA.defaultInstance().histogram1D(INV_MASS_REALLY_NO_CUTS, 750, 0.0000, 0.3000);
        AIDA.defaultInstance().histogram1D(INV_MASS_NO_CUTS, 750, 0.0000, 0.3000);
        AIDA.defaultInstance().histogram1D(INV_MASS_95, 750, 0.0000, 0.3000);
        AIDA.defaultInstance().histogram1D(INV_MASS_97, 750, 0.0000, 0.3000);
        AIDA.defaultInstance().histogram1D(INV_MASS_ME, 750, 0.0000, 0.3000);
    }
    
    public void setHodoscopeXLowerBound(double bound) {
        coptXLowerBound = bound;
    }
    
    public void setChiSquaredUpperBound(double bound) {
        chiSquaredUpperBound = bound;
    }
    
    public void setGtpClusterCollectionName(String collectionName) {
        gtpClusterCollectionName = collectionName;
    }
    
    public void setGblTrackCollectionName(String collectionName) {
        gblTrackCollectionName = collectionName;
    }
    
    public void setHodoscopeHitCollectionName(String collectionName) {
        hodoscopeHitCollectionName = collectionName;
    }
    
    public void setHodoscopeHitEnergyLowerBound(double bound) {
        hodoscopeHitEnergyLowerBound = bound / 1000.0;
    }
    
    public void setMatchingBoundsPositronMomentum(String coeffs) {
        boundsPositronP = parseDoubleArray(coeffs);
    }
    
    public void setMatchingBoundsElectronMomentum(String coeffs) {
        boundsElectronP = parseDoubleArray(coeffs);
    }
    
    public void setMatchingBoundsTopPositronLowerBound(String coeffs) {
        boundsTopPositron[0] = parseDoubleArray(coeffs);
    }
    
    public void setMatchingBoundsTopPositronUpperBound(String coeffs) {
        boundsTopPositron[1] = parseDoubleArray(coeffs);
    }
    
    public void setMatchingBoundsBottomPositronLowerBound(String coeffs) {
        boundsBotPositron[0] = parseDoubleArray(coeffs);
    }
    
    public void setMatchingBoundsBottomPositronUpperBound(String coeffs) {
        boundsBotPositron[1] = parseDoubleArray(coeffs);
    }
    
    public void setMatchingBoundsTopElectronLowerBound(String coeffs) {
        boundsTopElectron[0] = parseDoubleArray(coeffs);
    }
    
    public void setMatchingBoundsTopElectronUpperBound(String coeffs) {
        boundsTopElectron[1] = parseDoubleArray(coeffs);
    }
    
    public void setMatchingBoundsBottomElectronLowerBound(String coeffs) {
        boundsBotElectron[0] = parseDoubleArray(coeffs);
    }
    
    public void setMatchingBoundsBottomElectronUpperBound(String coeffs) {
        boundsBotElectron[1] = parseDoubleArray(coeffs);
    }
    
    private static final double[] parseDoubleArray(String text) {
        StringBuffer readBuffer = new StringBuffer();
        List<Double> valList = new ArrayList<Double>();
        for(int i = 0; i < text.length(); i++) {
            if(text.charAt(i) == ',') {
                valList.add(Double.parseDouble(readBuffer.toString()));
                readBuffer.delete(0, readBuffer.length());
            } else {
                if(!Character.isWhitespace(text.charAt(i))) { readBuffer.append(text.charAt(i)); }
            }
        }
        valList.add(Double.parseDouble(readBuffer.toString()));
        
        double[] array = new double[valList.size()];
        for(int i = 0; i < valList.size(); i++) {
            array[i] = valList.get(i).doubleValue();
        }
        
        return array;
    }
    
    public void setOutputDirectory(String dir) {
        outputDirectory = dir;
    }
    
    private static final void fillChiSquaredDistro(EventHeader event, String gblTrackCollectionName) {
        List<Track> gblTracks = TriggerTuningUtilityModule.getCollection(event, gblTrackCollectionName, Track.class);
        for(Track track : gblTracks) {
            AIDA.defaultInstance().histogram1D(CHI_SQUARED).fill(track.getChi2());
        }
    }
    
    private static final void fillClusterPositionDistro(EventHeader event, String gtpClusterCollectionName) {
        List<Cluster> gtpClusters = TriggerTuningUtilityModule.getCollection(event, gtpClusterCollectionName, Cluster.class);
        for(Cluster cluster : gtpClusters) {
            AIDA.defaultInstance().histogram2D(DEBUG_CLUSTER_POSITION).fill(cluster.getPosition()[0], cluster.getPosition()[1]);
        }
    }
    
    private static final void fillHodoscopeHitEnergyDistro(double[][][] hodoHits) {
        for(int ix = 0; ix < 5; ix++) {
            for(int iy = 0; iy < 2; iy++) {
                for(int iz = 0; iz < 2; iz++) {
                    if(hodoHits[ix][iy][iz] != 0) {
                        AIDA.defaultInstance().histogram1D(HODOSCOPE_SCINTILLATOR_ENERGY).fill(hodoHits[ix][iy][iz] * 1000);
                    }
                }
            }
        }
    }
    
    private static final void fillHodoscopeLayerToLayerDistro(double[][][] hodoHits, double hodoscopeHitEnergyLowerBound) {
        // Get the hodoscope cluster energies.
        double[][][] clusterEnergies = TriggerTuningUtilityModule.getHodoscopeClusterEnergies(hodoHits);
        
        // Iterate over each layer 1 hit.
        for(int l1y = 0; l1y < 2; l1y++) {
            boolean isTop = (l1y == 0);
            
            for(int l1x = 0; l1x < 9; l1x++) {
                // If the hit is below the energy cut threshold, do
                // not plot it.
                if(clusterEnergies[l1x][l1y][0] < hodoscopeHitEnergyLowerBound) { continue; }
                
                // Iterate over each layer 2 hit and plot the
                // correlation. Only consider hits on the top if the
                // L1 hit is on the top, and vice versa if it is on
                // the bottom.
                for(int l2x = 0; l2x < 9; l2x++) {
                    // Only plot for this hit if it exceeds the
                    // energy lower bound threshold.
                    if(clusterEnergies[l2x][l1y][1] >= hodoscopeHitEnergyLowerBound) {
                        AIDA.defaultInstance().histogram1D(getLayerToLayerPlotName(l1x, isTop)).fill(l2x + 1);
                    }
                }
            }
        }
    }
    
    private static final void fillHodoscopeLayerToCalorimeterDistro(double[][][] hodoHits, List<Cluster> gtpClusters, double hodoscopeHitEnergyLowerBound) {
        // Get the hodoscope cluster energies.
        double[][][] clusterEnergies = TriggerTuningUtilityModule.getHodoscopeClusterEnergies(hodoHits);
        
        // Get the highest energy calorimeter clusters for the top
        // and bottom on the positron side of the calorimeter.
        Cluster bestTopCluster = null;
        Cluster bestBottomCluster = null;
        for(Cluster cluster : gtpClusters) {
            // Skip clusters that are not on the positron side of the
            // calorimeter.
            if(!TriggerTuningUtilityModule.isPositronSideCluster(cluster)) { continue; }
            
            // Track the cluster which has the highest energy for the
            // top and bottom  of the calorimeter.
            if(TriggerTuningUtilityModule.isTopCluster(cluster)) {
                bestTopCluster = TriggerTuningUtilityModule.getHighestEnergyCluster(bestTopCluster, cluster);
            } else {
                bestBottomCluster = TriggerTuningUtilityModule.getHighestEnergyCluster(bestBottomCluster, cluster);
            }
        }
        
        // Iterate over all hodoscope hits/clusters.
        for(int layer = 0; layer < 2; layer++) {
            for(int iy = 0; iy < 2; iy++) {
                // If there does not exist a best cluster for the
                // current half of the calorimeter, skip this step.
                boolean isTop = (iy == 0);
                if(isTop && bestTopCluster == null) { continue; }
                else if(!isTop && bestBottomCluster == null) { continue; }
                
                // Correlate each hodoscope cluster with the
                // appropriate calorimeter cluster. Top clusters
                // match with top clusters and vice versa. Skip
                // clusters that are below the threshold.
                for(int ix = 0; ix < 9; ix++) {
                    // If the cluster is below the energy cut
                    // threshold, do not plot it.
                    if(clusterEnergies[ix][iy][layer] < hodoscopeHitEnergyLowerBound) { continue; }
                    
                    // Plot the correlation with the best cluster.
                    if(isTop) { AIDA.defaultInstance().histogram1D(getLayerToCalorimeterPlotName(ix, isTop, layer)).fill(TriggerModule.getClusterXIndex(bestTopCluster)); }
                    else { AIDA.defaultInstance().histogram1D(getLayerToCalorimeterPlotName(ix, isTop, layer)).fill(TriggerModule.getClusterXIndex(bestBottomCluster)); }
                    
                    /*
                    for(Cluster cluster : gtpClusters) {
                        if(isTop) { AIDA.defaultInstance().histogram1D(getLayerToCalorimeterPlotName(ix, isTop, layer)).fill(TriggerModule.getClusterXIndex(cluster)); }
                        else { AIDA.defaultInstance().histogram1D(getLayerToCalorimeterPlotName(ix, isTop, layer)).fill(TriggerModule.getClusterXIndex(cluster)); }
                    }
                    */
                }
            }
        }
    }
    
    private static final void fillHodoscopeTruthEnergyDistros(List<SimTrackerHit> truthHits, HodoscopeDetectorElement hodoscopeDetectorElement) {
        // Plot the hits individually.
        for(org.lcsim.event.SimTrackerHit hit : truthHits) {
            AIDA.defaultInstance().histogram1D(HODOSCOPE_TRUTH_ENERGY).fill(hit.getdEdx() * 1000);
        }
        
        // Plot the compiled hits.
        double[][][] compiledTruthEnergies = TriggerTuningUtilityModule.getCompiledHodoscopeEnergies(truthHits, hodoscopeDetectorElement);
        for(int ix = 0; ix < 5; ix++) {
            for(int iy = 0; iy < 2; iy++) {
                for(int iz = 0; iz < 2; iz++) {
                    if(compiledTruthEnergies[ix][iy][iz] != 0) {
                        AIDA.defaultInstance().histogram1D(HODOSCOPE_TRUTH_COMP_ENERGY).fill(compiledTruthEnergies[ix][iy][iz] * 1000);
                    }
                }
            }
        }
    }
    
    private final void fillTrackPositionDistro(EventHeader event, String gblTrackCollectionName) {
        List<Track> gblTracks = TriggerTuningUtilityModule.getCollection(event, gblTrackCollectionName, Track.class);
        for(Track track : gblTracks) {
            double[] r = TriggerTuningUtilityModule.getTrackPositionAtCalorimeterFace(track, fieldMap);
            AIDA.defaultInstance().histogram2D(DEBUG_TRACK_POSITION).fill(r[0], r[1]);
        }
    }
    
    private static final String getClusterTrackMatchingBasePlotName(boolean isTop, boolean isPositive, boolean isRecon) {
        if(isRecon) {
            return "Debug/Cluster-Track Matching/" + (isPositive ? "Positive/" : "Negative/") + (isTop ? "Top/" : "Bottom/") + "Momentum vs. #Delta ";
        } else {
            return "Cluster-Track Matching/" + (isPositive ? "Positive/" : "Negative/") + (isTop ? "Top/" : "Bottom/") + "Momentum vs. #Delta ";
        }
    }
    
    private static final String getClusterTrackMatchingXPlotName(boolean isTop, boolean isPositive, boolean isRecon) {
        return getClusterTrackMatchingBasePlotName(isTop, isPositive, isRecon) + "x";
    }
    
    private static final String getClusterTrackMatchingYPlotName(boolean isTop, boolean isPositive, boolean isRecon) {
        return getClusterTrackMatchingBasePlotName(isTop, isPositive, isRecon) + "y";
    }
    
    private static final String getClusterTrackMatchingRPlotName(boolean isTop, boolean isPositive, boolean isRecon) {
        return getClusterTrackMatchingBasePlotName(isTop, isPositive, isRecon) + "r";
    }
    
    private static final String getLayerToLayerPlotName(int clusterIndex, boolean isTop) {
        final String[] clusterNames = { "1", "1&2", "2", "2&3", "3", "3&4", "4", "4&5", "5" };
        return "Hodoscope/Layer Correlation/" + (isTop ? "Top/" : "Bottom/") + "Scintillator " + clusterNames[clusterIndex] + " Layer Correlation Distribution";
    }
    
    private static final String getLayerToCalorimeterPlotName(int clusterIndex, boolean isTop, int layer) {
        final String[] clusterNames = { "1", "1&2", "2", "2&3", "3", "3&4", "4", "4&5", "5" };
        return "Hodoscope/Calorimeter Correlation/" + (isTop ? "Top/" : "Bottom/") + "Layer " + (layer + 1) + "/" + "Scintillator " + clusterNames[clusterIndex]
                + " Calorimeter Correlation Distribution";
    }
    
    private void performHodoscopeAnalysis(EventHeader event) {
        // Get the hodoscope hits.
        List<CalorimeterHit> hodoHits = TriggerTuningUtilityModule.getCollection(event, hodoscopeHitCollectionName, CalorimeterHit.class);
        List<Cluster> gtpClusters = TriggerTuningUtilityModule.getCollection(event, gtpClusterCollectionName, Cluster.class);
        
        // Some hodoscope scintillators have multiple channels, which
        // allows for multiple hits per channel. Combine these so
        // that there is only one energy per scintillator.
        double[][][] compiledHodoHits = TriggerTuningUtilityModule.getCompiledHodoscopeEnergies(hodoHits, hodoscopeChannelMap);
        
        // Plot the distribution of hodoscope hit energies.
        fillHodoscopeHitEnergyDistro(compiledHodoHits);
        
        // Fill the debug truth comparison distributions.
        fillHodoscopeTruthEnergyDistros(TriggerTuningUtilityModule.getCollection(event, "HodoscopeHits", SimTrackerHit.class), hodoscopeDetectorElement);
        
        // Fill the layer-to-layer correlation plots.
        fillHodoscopeLayerToLayerDistro(compiledHodoHits, hodoscopeHitEnergyLowerBound);
        
        // Fill the layer-to-calorimeter correlation plots.
        fillHodoscopeLayerToCalorimeterDistro(compiledHodoHits, gtpClusters, hodoscopeHitEnergyLowerBound);
    }
    
    private void performTrackAnalysis(EventHeader event) {
        // Get the necessary collections.
        List<Track> gblTracks = TriggerTuningUtilityModule.getCollection(event, gblTrackCollectionName, Track.class);
        List<Cluster> gtpClusters = TriggerTuningUtilityModule.getCollection(event, gtpClusterCollectionName, Cluster.class);
        List<Cluster> reconClusters = TriggerTuningUtilityModule.getCollection(event, "EcalClustersCorr", Cluster.class);
        
        // Iterate over the tracks and plot their momenta based on
        // whether a given track is positive or negative and a top or
        // bottom track.
        for(Track gblTrack : gblTracks) {
            // TODO: These methods really need to be updated to the correct version.
            // Get the track position at the calorimeter face.
            double[] trackR = TriggerTuningUtilityModule.getTrackPositionAtCalorimeterFace(gblTrack, fieldMap);
            double trackP = TriggerTuningUtilityModule.getMomentumMagnitude(gblTrack, fieldMap);
            
            // Plot the difference in x and y between the track and
            // all possible clusters.
            for(Cluster gtpCluster : gtpClusters) {
                double[] clusterR = gtpCluster.getPosition();
                double deltaX = clusterR[0] - trackR[0];
                double deltaY = clusterR[1] - trackR[1];
                double deltaR = Math.sqrt(Math.pow(deltaX, 2) + Math.pow(deltaY, 2));
                
                boolean isTop = trackR[1] > 0;
                boolean isPositive = TriggerTuningUtilityModule.isPositive(gblTrack);
                
                AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(isTop, isPositive, false)).fill(trackP, deltaR);
                AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(isTop, isPositive, false)).fill(trackP, Math.abs(deltaX));
                AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(isTop, isPositive, false)).fill(trackP, Math.abs(deltaY));
            }
            
            // As a debugging step, repeat this for recon clusters.
            for(Cluster reconCluster : reconClusters) {
                double[] clusterR = reconCluster.getPosition();
                double deltaX = clusterR[0] - trackR[0];
                double deltaY = clusterR[1] - trackR[1];
                double deltaR = Math.sqrt(Math.pow(deltaX, 2) + Math.pow(deltaY, 2));
                
                boolean isTop = trackR[1] > 0;
                boolean isPositive = TriggerTuningUtilityModule.isPositive(gblTrack);
                
                AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(isTop, isPositive, true)).fill(trackP, deltaR);
                AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(isTop, isPositive, true)).fill(trackP, Math.abs(deltaX));
                AIDA.defaultInstance().histogram2D(getClusterTrackMatchingRPlotName(isTop, isPositive, true)).fill(trackP, Math.abs(deltaY));
            }
        }
    }
    
    private void performSinglesTriggerAnalysis(List<Pair<Cluster, Track>> pairList) {
        // Populate the singles cut distribution for all clusters.
        for(Pair<Cluster, Track> pair : pairList) {
            AIDA.defaultInstance().histogram1D(SINGLES_MULTIPLICITY).fill(TriggerModule.getClusterHitCount(pair.getFirstElement()));
            AIDA.defaultInstance().histogram1D(SINGLES_SEED_ENERGY).fill(TriggerModule.getValueClusterSeedEnergy(pair.getFirstElement()));
            AIDA.defaultInstance().histogram1D(SINGLES_TOTAL_ENERGY).fill(TriggerModule.getValueClusterTotalEnergy(pair.getFirstElement()));
            AIDA.defaultInstance().histogram2D(SINGLES_POSITION).fill(TriggerModule.getClusterXIndex(pair.getFirstElement()),
                    TriggerModule.getClusterYIndex(pair.getFirstElement()));
        }
    }
    
    private void performPairTriggerAnalysis(List<Pair<Cluster, Track>> pairList) {
        // Only plot clusters once for singles cuts.
        Set<Cluster> plottedClusterSet = new HashSet<Cluster>();
        
        // Consider all possible pairs of clusters, and plot for
        // those that meet the pair conditions.
        for(int i = 0; i < pairList.size(); i++) {
            // Get the pair.
            Pair<Cluster, Track> pair1 = pairList.get(i);
            
            // Determine the charge of the track and the position of
            // the cluster.
            boolean isTop = TriggerModule.getClusterYIndex(pair1.getFirstElement()) > 0;
            boolean isPositive = TriggerTuningUtilityModule.isPositive(pair1.getSecondElement());
            
            // Loop over all remaining pairs.
            for(int j = i + 1; j < pairList.size(); j++) {
                // Get the pair.
                Pair<Cluster, Track> pair2 = pairList.get(j);
                
                // Only consider pairs of positive/negative tracks
                // and top/bottom clusters.
                if((TriggerTuningUtilityModule.isPositive(pair2.getSecondElement()) && isPositive)
                        || (TriggerTuningUtilityModule.isNegative(pair2.getSecondElement()) && !isPositive)) {
                    continue;
                }
                if((TriggerModule.getClusterYIndex(pair2.getFirstElement()) > 0 && isTop)
                        || (TriggerModule.getClusterYIndex(pair2.getFirstElement()) < 0 && !isTop)) {
                    continue;
                }
                
                // Get the cluster pair.
                Cluster[] pair = new Cluster[] { pair1.getFirstElement(), pair2.getFirstElement() };
                
                // Fill the pair trigger singles plots.
                for(Cluster cluster : pair) {
                    if(!plottedClusterSet.contains(cluster)) {
                        plottedClusterSet.add(cluster);
                        AIDA.defaultInstance().histogram1D(PAIR_MULTIPLICITY).fill(TriggerModule.getClusterHitCount(cluster));
                        AIDA.defaultInstance().histogram1D(PAIR_TOTAL_ENERGY).fill(TriggerModule.getValueClusterTotalEnergy(cluster));
                        AIDA.defaultInstance().histogram2D(PAIR_POSITION).fill(TriggerModule.getClusterXIndex(cluster), TriggerModule.getClusterYIndex(cluster));
                    }
                }
                
                // Fill the pair trigger pair plots.
                AIDA.defaultInstance().histogram1D(PAIR_ENERGY_SUM).fill(TriggerModule.getValueEnergySum(pair));
                AIDA.defaultInstance().histogram1D(PAIR_ENERGY_DIFF).fill(TriggerModule.getValueEnergyDifference(pair));
                AIDA.defaultInstance().histogram1D(PAIR_COPLANARITY).fill(TriggerModule.getValueCoplanarity(pair));
                
                // Get the lowest energy cluster.
                Cluster lowestEnergyCluster = pair[0];
                if(TriggerModule.getValueClusterTotalEnergy(pair[1]) < TriggerModule.getValueClusterTotalEnergy(pair[0])) {
                    lowestEnergyCluster = pair[1];
                }
                
                // Get deltaR for the lowest energy cluster.
                double r = Math.sqrt(Math.pow(TriggerModule.getClusterX(lowestEnergyCluster), 2) + Math.pow(TriggerModule.getClusterY(lowestEnergyCluster), 2));
                
                // Fill the energy slope plot.
                AIDA.defaultInstance().histogram2D(PAIR_ENERGY_SLOPE).fill(r, TriggerModule.getValueClusterTotalEnergy(lowestEnergyCluster));
            }
        }
        
        for(Pair<Cluster, Track> pair : pairList) {
            AIDA.defaultInstance().histogram1D(SINGLES_MULTIPLICITY).fill(TriggerModule.getClusterHitCount(pair.getFirstElement()));
            AIDA.defaultInstance().histogram1D(SINGLES_SEED_ENERGY).fill(TriggerModule.getValueClusterSeedEnergy(pair.getFirstElement()));
            AIDA.defaultInstance().histogram1D(SINGLES_TOTAL_ENERGY).fill(TriggerModule.getValueClusterTotalEnergy(pair.getFirstElement()));
            AIDA.defaultInstance().histogram2D(SINGLES_POSITION).fill(TriggerModule.getClusterXIndex(pair.getFirstElement()),
                    TriggerModule.getClusterYIndex(pair.getFirstElement()));
        }
    }
    
    private void performCOPTriggerAnalysis(List<Pair<Cluster, Track>> pairList) {
        // Populate the cut distributions for all positron clusters.
        for(Pair<Cluster, Track> pair : pairList) {
            if(TriggerTuningUtilityModule.isPositive(pair.getSecondElement())) {
                AIDA.defaultInstance().histogram1D(COPT_MULTIPLICITY).fill(TriggerModule.getClusterHitCount(pair.getFirstElement()));
                AIDA.defaultInstance().histogram1D(COPT_SEED_ENERGY).fill(TriggerModule.getValueClusterSeedEnergy(pair.getFirstElement()));
                AIDA.defaultInstance().histogram1D(COPT_TOTAL_ENERGY).fill(TriggerModule.getValueClusterTotalEnergy(pair.getFirstElement()));
                AIDA.defaultInstance().histogram2D(COPT_POSITION).fill(TriggerModule.getClusterXIndex(pair.getFirstElement()),
                        TriggerModule.getClusterYIndex(pair.getFirstElement()));
                AIDA.defaultInstance().histogram2D(COPT_ENERGY_POSITION).fill(TriggerModule.getClusterXIndex(pair.getFirstElement()),
                        TriggerModule.getValueClusterTotalEnergy(pair.getFirstElement()));
                
                AIDA.defaultInstance().histogram2D(DEBUG_MOMENTUM_POSITION).fill(TriggerModule.getClusterXIndex(pair.getFirstElement()),
                        TriggerTuningUtilityModule.getMomentumMagnitude(pair.getSecondElement(), fieldMap));
                AIDA.defaultInstance().histogram2D(DEBUG_PHI_POSITION).fill(TriggerModule.getClusterXIndex(pair.getFirstElement()),
                        1000 * TriggerTuningUtilityModule.getTrackPhi(pair.getSecondElement()));
            }
        }
    }
    
    private void performInvariantMassAnalysis(List<Pair<Cluster, Track>> pairList) {
        // If there are fewer than two tracks, then no analysis may
        // be performed.
        if(pairList.size() < 2) { return; }
        
        // Plot top/bottom, positive/negative track pairs.
        for(int i = 0; i < pairList.size(); i++) {
            for(int j = i + 1; j < pairList.size(); j++) {
                if(!isValidPair(pairList.get(i).getSecondElement(), pairList.get(j).getSecondElement())) { continue; };
                List<Hep3Vector> momenta = new ArrayList<Hep3Vector>(2);
                momenta.add(new BasicHep3Vector(TriggerTuningUtilityModule.getMomentum(pairList.get(i).getSecondElement(), fieldMap)));
                momenta.add(new BasicHep3Vector(TriggerTuningUtilityModule.getMomentum(pairList.get(j).getSecondElement(), fieldMap)));
                double invariantMass = getInvariantMass(momenta);
                AIDA.defaultInstance().histogram1D(INV_MASS_NO_CUTS).fill(invariantMass);
                
                // Define the COPT energy vs. position cut
                // coefficients.
                final double[] firCoeff95 = { 2.186811145510821,  -0.18388028895768568, 0.006550567595459063,  -0.00007997936016511498 };
                final double[] firCoeff97 = { 1.9004592363261057, -0.1715872606352472,  0.005825593395252756,  -0.00005360623781676235 };
                final double[] firCoeff99 = { 1.2412860734037259, -0.15417121549474627, 0.0074262012497307685, -0.00011773085302497325 };
                
                // Get the clusters.
                Cluster[] cluster = { pairList.get(i).getFirstElement(), pairList.get(j).getFirstElement() };
                int[] ix = { TriggerModule.getClusterXIndex(cluster[0]), TriggerModule.getClusterXIndex(cluster[1]) };
                
                // Get the positron cluster. This will be the one
                // with ix >= 2.
                int positiveIndex = -1;
                if(ix[0] >= 2 && ix[0] > ix[1]) { positiveIndex = 0; }
                else if(ix[1] >= 2 && ix[1] > ix[0]) { positiveIndex = 1; }
                
                // If there is a positive cluster, plot it as appropriate.
                if(positiveIndex != -1) {
                    // Perform the COPT cuts and plot the track pair
                    // if it passes.
                    double threshold95 = TriggerTuningUtilityModule.polynomial(firCoeff95, ix[positiveIndex]);
                    double threshold97 = TriggerTuningUtilityModule.polynomial(firCoeff97, ix[positiveIndex]);
                    double threshold99 = TriggerTuningUtilityModule.polynomial(firCoeff99, ix[positiveIndex]);
                    
                    // Fill the relevant plots.
                    if(cluster[positiveIndex].getEnergy() >= threshold95) {
                        AIDA.defaultInstance().histogram1D(INV_MASS_95).fill(invariantMass);
                    }
                    if(cluster[positiveIndex].getEnergy() >= threshold97) {
                        AIDA.defaultInstance().histogram1D(INV_MASS_97).fill(invariantMass);
                    }
                    if(cluster[positiveIndex].getEnergy() >= threshold99) {
                        AIDA.defaultInstance().histogram1D(INV_MASS_ME).fill(invariantMass);
                    }
                }
            }
        }
    }
    
    private void performInvariantMassAnalysisMatched(List<Pair<Cluster, Track>> pairList) {
        // If there are fewer than two tracks, then no analysis may
        // be performed.
        if(pairList.size() < 2) { return; }
        
        // Plot top/bottom, positive/negative track pairs.
        for(int i = 0; i < pairList.size(); i++) {
            for(int j = i + 1; j < pairList.size(); j++) {
                if(!isValidPair(pairList.get(i).getSecondElement(), pairList.get(j).getSecondElement())) { continue; };
                List<Hep3Vector> momenta = new ArrayList<Hep3Vector>(2);
                momenta.add(new BasicHep3Vector(TriggerTuningUtilityModule.getMomentum(pairList.get(i).getSecondElement(), fieldMap)));
                momenta.add(new BasicHep3Vector(TriggerTuningUtilityModule.getMomentum(pairList.get(j).getSecondElement(), fieldMap)));
                double invariantMass = getInvariantMass(momenta);
                AIDA.defaultInstance().histogram1D(INV_MASS_NO_CUTS).fill(invariantMass);
            }
        }
    }
    
    private static final double getInvariantMass(List<Hep3Vector> momenta) {
        // Calculate gamma * m of each origin particle and sum them.
        double energySum = 0.0;
        final double m = 0.000511;
        for(Hep3Vector momentum : momenta) {
            energySum += Math.sqrt(Math.pow(m, 2) + momentum.magnitudeSquared());
        }
        
        // Get the vector sum of the particle momenta.
        Hep3Vector momentumSum = new BasicHep3Vector(0, 0, 0);
        for(Hep3Vector momentum : momenta) {
            momentumSum = VectorArithmetic.add(momentumSum, momentum);
        }
        
        // Calculate and plot the invariant mass.
        return Math.sqrt(Math.pow(energySum, 2) - momentumSum.magnitudeSquared());
    }
    
    private static final boolean isValidPair(Track track0, Track track1) {
        boolean[] isTop = { TriggerTuningUtilityModule.isTopTrack(track0), TriggerTuningUtilityModule.isTopTrack(track1) };
        boolean[] isPositive = { TriggerTuningUtilityModule.isPositive(track0), TriggerTuningUtilityModule.isPositive(track1) };
        
        boolean isTopBot = (isTop[0] && !isTop[1]) || (!isTop[0] && isTop[1]);
        boolean isPosNeg = (isPositive[0] && !isPositive[1]) || (!isPositive[0] && isPositive[1]);
        return isTopBot && isPosNeg;
    }
    
    private static final boolean isValidPair(MCParticle part0, MCParticle part1) {
        return (part0.getCharge() > 0 && part1.getCharge() < 0) || (part0.getCharge() < 0 && part1.getCharge() > 0);
    }
    
    private void performInvariantMassAnalysisNoClusters(List<Track> tracks) {
        // If there are fewer than two tracks, then no analysis may
        // be performed.
        if(tracks.size() < 2) { return; }
        
        // Plot top/bottom, positive/negative track pairs.
        for(int i = 0; i < tracks.size(); i++) {
            for(int j = i + 1; j < tracks.size(); j++) {
                if(!isValidPair(tracks.get(i), tracks.get(j))) { continue; };
                List<Hep3Vector> momenta = new ArrayList<Hep3Vector>(2);
                momenta.add(new BasicHep3Vector(TriggerTuningUtilityModule.getMomentum(tracks.get(i), fieldMap)));
                momenta.add(new BasicHep3Vector(TriggerTuningUtilityModule.getMomentum(tracks.get(j), fieldMap)));
                double invariantMass = getInvariantMass(momenta);
                AIDA.defaultInstance().histogram1D(INV_MASS_REALLY_NO_CUTS).fill(invariantMass);
            }
        }
    }
    
    private boolean isOriginParticle(MCParticle particle) {
        boolean isPosEle = (particle.getPDGID() == 11) || (particle.getPDGID() == -11);
        boolean hasAPrimeParent = (!particle.getParents().isEmpty() && particle.getParents().get(0).getPDGID() == 622);
        return isPosEle && hasAPrimeParent;
    }
    
    private void performInvariantMassAnalysisTruth(List<MCParticle> particles) {
        // Get all of the origin particles.
        List<MCParticle> originParticles = new ArrayList<MCParticle>();
        for(MCParticle particle : particles) {
            if(isOriginParticle(particle)) {
                originParticles.add(particle);
            }
        }
        
        // Plot top/bottom, positive/negative track pairs.
        for(int i = 0; i < originParticles.size(); i++) {
            for(int j = i + 1; j < originParticles.size(); j++) {
                if(!isValidPair(originParticles.get(i), originParticles.get(j))) { continue; };
                List<Hep3Vector> momenta = new ArrayList<Hep3Vector>(2);
                momenta.add(originParticles.get(i).getMomentum());
                momenta.add(originParticles.get(j).getMomentum());
                double invariantMass = getInvariantMass(momenta);
                AIDA.defaultInstance().histogram1D(INV_MASS_TRUTH).fill(invariantMass);
            }
        }
    }
}