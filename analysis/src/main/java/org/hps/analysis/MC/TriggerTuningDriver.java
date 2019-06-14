package org.hps.analysis.MC;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.conditions.hodoscope.HodoscopeChannel;
import org.hps.conditions.hodoscope.HodoscopeChannel.HodoscopeChannelCollection;
import org.hps.detector.hodoscope.HodoscopeDetectorElement;
import org.hps.readout.triggerstudies.Coordinate;
import org.hps.readout.triggerstudies.PlotToTextModule;
import org.hps.record.triggerbank.TriggerModule;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.FieldMap;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

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
    
    private String outputDirectory = ".";
    private FieldMap fieldMap = null;
    private HodoscopeDetectorElement hodoscopeDetectorElement;
    private Map<Long, HodoscopeChannel> hodoscopeChannelMap = new HashMap<Long, HodoscopeChannel>();
    
    private static final String CHI_SQUARED = "Cluster\\Track-Matching/Chi Squared Distribution";
    private static final String HODOSCOPE_SCINTILLATOR_ENERGY = "Hodoscope/Scintillator Energy Distribution";
    
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
            
            // Cluster/track matching plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(true, true)),
                        new File(outputDirectory + File.separator + "tuning_clusterTrackMatching_TPR.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(true, false)),
                    new File(outputDirectory + File.separator + "tuning_clusterTrackMatching_TNR.dat"));
        
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(false, true)),
                    new File(outputDirectory + File.separator + "tuning_clusterTrackMatching_BPR.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(false, false)),
                    new File(outputDirectory + File.separator + "tuning_clusterTrackMatching_BNR.dat"));
            
            // Define the cluster/track matching fit threshold.
            double clusterTrackMatchingFitThreshold = 0.99;
            
            // Get the upper and lower bounds of the distribution for
            // each data type.
            // Top Positrons
            Coordinate[] positronTopXLowerBound =
                    TriggerTuningUtilityModule.getEnergyThreshold(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(true, true)),
                    0, 0, Integer.MAX_VALUE, clusterTrackMatchingFitThreshold, true);
            Coordinate[] positronTopXUpperBound =
                    TriggerTuningUtilityModule.getEnergyThreshold(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(true, true)),
                    0, 0, Integer.MAX_VALUE, clusterTrackMatchingFitThreshold, false);
            
            // Bottom Positrons
            Coordinate[] positronBotXLowerBound =
                    TriggerTuningUtilityModule.getEnergyThreshold(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(false, true)),
                    0, 0, Integer.MAX_VALUE, clusterTrackMatchingFitThreshold, true);
            Coordinate[] positronBotXUpperBound =
                    TriggerTuningUtilityModule.getEnergyThreshold(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(false, true)),
                    0, 0, Integer.MAX_VALUE, clusterTrackMatchingFitThreshold, false);
            
            // Top Electrons
            Coordinate[] electronTopXLowerBound =
                    TriggerTuningUtilityModule.getEnergyThreshold(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(true, false)),
                    0, Integer.MAX_VALUE, 0, clusterTrackMatchingFitThreshold, true);
            Coordinate[] electronTopXUpperBound =
                    TriggerTuningUtilityModule.getEnergyThreshold(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(true, false)),
                    0, Integer.MAX_VALUE, 0, clusterTrackMatchingFitThreshold, false);
            
            // Bottom Electrons
            Coordinate[] electronBotXLowerBound =
                    TriggerTuningUtilityModule.getEnergyThreshold(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(false, false)),
                    0, Integer.MAX_VALUE, 0, clusterTrackMatchingFitThreshold, true);
            Coordinate[] electronBotXUpperBound =
                    TriggerTuningUtilityModule.getEnergyThreshold(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(false, false)),
                    0, Integer.MAX_VALUE, 0, clusterTrackMatchingFitThreshold, false);
            
            // Write the cluster/track matching boundary files.
            PlotToTextModule.writeCoordinateSet(positronTopXLowerBound, new File(outputDirectory + File.separator + "fit_clusterTrackMatching_TPRL.dat"));
            PlotToTextModule.writeCoordinateSet(positronTopXUpperBound, new File(outputDirectory + File.separator + "fit_clusterTrackMatching_TPRU.dat"));
            
            PlotToTextModule.writeCoordinateSet(positronBotXLowerBound, new File(outputDirectory + File.separator + "fit_clusterTrackMatching_BPRL.dat"));
            PlotToTextModule.writeCoordinateSet(positronBotXUpperBound, new File(outputDirectory + File.separator + "fit_clusterTrackMatching_BPRU.dat"));
            
            PlotToTextModule.writeCoordinateSet(electronTopXLowerBound, new File(outputDirectory + File.separator + "fit_clusterTrackMatching_TNRL.dat"));
            PlotToTextModule.writeCoordinateSet(electronTopXUpperBound, new File(outputDirectory + File.separator + "fit_clusterTrackMatching_TNRU.dat"));
            
            PlotToTextModule.writeCoordinateSet(electronBotXLowerBound, new File(outputDirectory + File.separator + "fit_clusterTrackMatching_BNRL.dat"));
            PlotToTextModule.writeCoordinateSet(electronBotXUpperBound, new File(outputDirectory + File.separator + "fit_clusterTrackMatching_BNRU.dat"));
            
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
        
        // Perform the cluster/track matching analysis. This does not
        // require that an event be analyzable to be useful.
        performTrackAnalysis(event);
        
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
            // Increment the event count.
            singlesEvents++;
        }
        
        // Perform pair trigger tuning.
        if(meetsPairCriteria) {
            // Increment the event count.
            pairEvents++;
        }
        
        // Perform triplet trigger tuning.
        if(meetsTripletCriteria) {
            // Increment the event count.
            tripletEvents++;
        }
        
        // Perform COPT trigger tuning.
        if(meetsCOPTCriteria) {
            // Increment the event count.
            coptEvents++;
        }
        
        // Perform hodoscope trigger tuning.
        if(meetsHodoscopeCriteria) {
            // Increment the event count.
            hodoscopeEvents++;
            
            // Perform hodoscope analysis.
            performHodoscopeAnalysis(event);
        }
    }
    
    @Override
    public void startOfData() {
        // Preliminary plots.
        AIDA.defaultInstance().histogram1D(CHI_SQUARED, 50, 0.0, 50.0);
        AIDA.defaultInstance().histogram1D(HODOSCOPE_SCINTILLATOR_ENERGY, 500, 0.000, 5.000);
        
        // Cluster/Track matching plots.
        for(int i = 0; i < 2; i++) {
            boolean isTop = (i == 0);
            for(int j = 0; j < 2; j++) {
                boolean isPositive = (j == 0);
                AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(isTop, isPositive), 250, 0.000, 5.000, 70, 0, 70);
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
        
        // Hodoscope Hit Truth Comparison Debug Plots
        AIDA.defaultInstance().histogram1D(HODOSCOPE_TRUTH_ENERGY, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram1D(HODOSCOPE_TRUTH_COMP_ENERGY, 500, 0.000, 5.000);
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
    
    public void setOutputDirectory(String dir) {
        outputDirectory = dir;
    }
    
    private static final void fillChiSquaredDistro(EventHeader event, String gblTrackCollectionName) {
        List<Track> gblTracks = TriggerTuningUtilityModule.getCollection(event, gblTrackCollectionName, Track.class);
        for(Track track : gblTracks) {
            AIDA.defaultInstance().histogram1D(CHI_SQUARED).fill(track.getChi2());
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
    
    private static final String getClusterTrackMatchingPlotName(boolean isTop, boolean isPositive) {
        return "Cluster-Track Matching/" + (isPositive ? "Positive/" : "Negative/") + (isTop ? "Top/" : "Bottom/") + "Momentum vs. #Delta r";
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
        
        // Iterate over the tracks and plot their momenta based on
        // whether a given track is positive or negative and a top or
        // bottom track.
        for(Track gblTrack : gblTracks) {
            // TODO: These methods really need to be updated to the correct version.
            // Get the track position at the calorimeter face.
            double[] trackR = TriggerTuningUtilityModule.getTrackPositionAtCalorimeterFace(gblTrack);
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
                
                AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(isTop, isPositive)).fill(trackP, deltaR);
            }
        }
    }
}