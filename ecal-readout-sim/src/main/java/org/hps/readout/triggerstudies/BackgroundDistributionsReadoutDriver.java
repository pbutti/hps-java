package org.hps.readout.triggerstudies;

import java.awt.Point;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

import org.hps.readout.ReadoutDataManager;
import org.hps.readout.ReadoutDriver;
import org.hps.record.triggerbank.TriggerModule;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.geometry.IDDecoder;
import org.lcsim.util.aida.AIDA;

public class BackgroundDistributionsReadoutDriver extends ReadoutDriver {
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
    
    private int instance = 0;
    private double localTime = 0.0;
    private String outputDirectory = ".";
    private String gtpClusterCollectionName = "EcalClustersGTP";
    private Set<Cluster> pairPlottedClusters = new HashSet<Cluster>();
    private LinkedList<Collection<Cluster>> clusterBuffer = new LinkedList<Collection<Cluster>>();
    
    @Override
    public void endOfData() {
        // Output the Mathematica dat files.
        try {
            // Singles trigger plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(SINGLES_POSITION),
                    new File(outputDirectory + File.separator + "tuning_singles_position_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(SINGLES_MULTIPLICITY),
                    new File(outputDirectory + File.separator + "tuning_singles_multiplicity_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(SINGLES_SEED_ENERGY),
                    new File(outputDirectory + File.separator + "tuning_singles_seedEnergy_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(SINGLES_TOTAL_ENERGY),
                    new File(outputDirectory + File.separator + "tuning_singles_totalEnergy_" + Integer.toString(instance) + ".dat"));
            
            // Pair trigger plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(PAIR_POSITION),
                    new File(outputDirectory + File.separator + "tuning_pair_position_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(PAIR_MULTIPLICITY),
                    new File(outputDirectory + File.separator + "tuning_pair_multiplicity_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(PAIR_TOTAL_ENERGY),
                    new File(outputDirectory + File.separator + "tuning_pair_totalEnergy_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(PAIR_ENERGY_SUM),
                    new File(outputDirectory + File.separator + "tuning_pair_energySum_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(PAIR_ENERGY_DIFF),
                    new File(outputDirectory + File.separator + "tuning_pair_energyDifference_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(PAIR_ENERGY_SLOPE),
                    new File(outputDirectory + File.separator + "tuning_pair_energySlope_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(PAIR_COPLANARITY),
                    new File(outputDirectory + File.separator + "tuning_pair_coplanarity_" + Integer.toString(instance) + ".dat"));
            
            // COP-trigger plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(COPT_POSITION),
                    new File(outputDirectory + File.separator + "tuning_copt_position_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(COPT_ENERGY_POSITION),
                    new File(outputDirectory + File.separator + "tuning_copt_energyVsPosition_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(COPT_MULTIPLICITY),
                    new File(outputDirectory + File.separator + "tuning_copt_multiplicity_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(COPT_SEED_ENERGY),
                    new File(outputDirectory + File.separator + "tuning_copt_seedEnergy_" + Integer.toString(instance) + ".dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(COPT_TOTAL_ENERGY),
                    new File(outputDirectory + File.separator + "tuning_copt_totalEnergy_" + Integer.toString(instance) + ".dat"));
        } catch(IOException e) {
            throw new RuntimeException(e);
        }
    }
    
    @Override
    public void process(EventHeader event) {
        // If no new clusters are ready, do nothing.
        if(!ReadoutDataManager.checkCollectionStatus(gtpClusterCollectionName, localTime + 4.0)) {
            return;
        }
        
        // Otherwise, get the clusters from the current clock-cycle.
        Collection<Cluster> gtpClusters = ReadoutDataManager.getData(localTime, localTime + 4.0, gtpClusterCollectionName, Cluster.class);
        
        // Add the cluster collection to the buffer. If the buffer is
        // now larger than 8 (+/- 16 ns), remove the oldest entry.
        // Also remove all of the removed buffer entry clusters from
        // the pair singles cuts plotted clusters set.
        clusterBuffer.addFirst(gtpClusters);
        if(clusterBuffer.size() > 8) {
            Collection<Cluster> oldBufferEntry = clusterBuffer.removeLast();
            pairPlottedClusters.removeAll(oldBufferEntry);
        }
        
        // Populate the trigger distributions.
        performSinglesTriggerAnalysis();
        performPairTriggerAnalysis();
        performCOPTriggerAnalysis();
        
        // Increment the local time.
        localTime += 4.0;
    }
    
    public void setInstance(int instance) {
        this.instance = instance;
    }
    
    public void setOutputDirectory(String dir) {
        outputDirectory = dir;
    }
    
    @Override
    public void startOfData() {
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
    }
    
    @Override
    protected double getTimeDisplacement() {
        return 0;
    }
    
    @Override
    protected double getTimeNeededForLocalOutput() {
        return 0;
    }
    
    private void performSinglesTriggerAnalysis() {
        // The singles trigger should always use the newest cluster
        // buffer entry.
        Collection<Cluster> clusters = clusterBuffer.getFirst();
        
        // Populate the singles cut distribution for all clusters.
        IDDecoder decoder = ReadoutDataManager.getIDDecoder("EcalHits");
        for(Cluster cluster : clusters) {
            decoder.setID(cluster.getCalorimeterHits().get(0).getCellID());
            AIDA.defaultInstance().histogram1D(SINGLES_MULTIPLICITY).fill(TriggerModule.getClusterHitCount(cluster));
            AIDA.defaultInstance().histogram1D(SINGLES_SEED_ENERGY).fill(TriggerModule.getValueClusterSeedEnergy(cluster));
            AIDA.defaultInstance().histogram1D(SINGLES_TOTAL_ENERGY).fill(TriggerModule.getValueClusterTotalEnergy(cluster));
            AIDA.defaultInstance().histogram2D(SINGLES_POSITION).fill(decoder.getValue("ix"), decoder.getValue("iy"));
            System.out.printf("<%3d, %2d>%n", decoder.getValue("ix"), decoder.getValue("iy"));
        }
    }
    
    private void performPairTriggerAnalysis() {
        // All cluster pairs are formed using one cluster from the
        // most recent cluster buffer entry to avoid duplicates.
        Collection<Cluster> pairingClusters = clusterBuffer.getFirst();
        
        // Form all possible pairs consisting of one cluster from the
        // most recent entry of the cluster buffer and one cluster
        // from any of the remaining buffers.
        IDDecoder decoder = ReadoutDataManager.getIDDecoder("EcalHits");
        for(Cluster pairingCluster : pairingClusters) {
            // Track whether the cluster is a top or bottom cluster.
            decoder.setID(pairingCluster.getCalorimeterHits().get(0).getCellID());
            Point pairingClusterIR = new Point(decoder.getValue("ix"), decoder.getValue("iy"));
            boolean pairingIsTop = (pairingClusterIR.y > 0);
            
            // Iterate over each entry in the cluster buffer. It is
            // allowed that clusters pair with others in the most
            // recent buffer entry, so long as it is not the same
            // cluster.
            for(Collection<Cluster> bufferEntry : clusterBuffer) {
                for(Cluster bufferCluster : bufferEntry) {
                    // Do not form "pairs" from the one cluster.
                    if(pairingCluster == bufferCluster) { continue; }
                    
                    // Only consider pairs where one cluster is a top
                    // cluster and one cluster is a bottom cluster.
                    decoder.setID(bufferCluster.getCalorimeterHits().get(0).getCellID());
                    Point bufferClusterIR = new Point(decoder.getValue("ix"), decoder.getValue("iy"));
                    boolean bufferIsTop = (bufferClusterIR.y > 0);
                    if((pairingIsTop && bufferIsTop) || (!pairingIsTop && !bufferIsTop)) { continue; }
                    
                    // Populate the singles plots with each cluster,
                    // but only if the cluster has not already been
                    // plotted.
                    if(!pairPlottedClusters.contains(pairingCluster)) {
                        pairPlottedClusters.add(pairingCluster);
                        AIDA.defaultInstance().histogram1D(PAIR_MULTIPLICITY).fill(TriggerModule.getClusterHitCount(pairingCluster));
                        AIDA.defaultInstance().histogram1D(PAIR_TOTAL_ENERGY).fill(TriggerModule.getValueClusterTotalEnergy(pairingCluster));
                        AIDA.defaultInstance().histogram2D(PAIR_POSITION).fill(pairingClusterIR.x, pairingClusterIR.y);
                    }
                    if(!pairPlottedClusters.contains(bufferCluster)) {
                        pairPlottedClusters.add(bufferCluster);
                        AIDA.defaultInstance().histogram1D(PAIR_MULTIPLICITY).fill(TriggerModule.getClusterHitCount(bufferCluster));
                        AIDA.defaultInstance().histogram1D(PAIR_TOTAL_ENERGY).fill(TriggerModule.getValueClusterTotalEnergy(bufferCluster));
                        AIDA.defaultInstance().histogram2D(PAIR_POSITION).fill(bufferClusterIR.x, bufferClusterIR.y);
                    }
                    
                    // Populate the pair cuts plots.
                    Cluster[] pair = new Cluster[] { pairingCluster, bufferCluster };
                    AIDA.defaultInstance().histogram1D(PAIR_ENERGY_SUM).fill(TriggerModule.getValueEnergySum(pair));
                    AIDA.defaultInstance().histogram1D(PAIR_ENERGY_DIFF).fill(TriggerModule.getValueEnergyDifference(pair));
                    AIDA.defaultInstance().histogram1D(PAIR_COPLANARITY).fill(TriggerModule.getValueCoplanarity(pairingClusterIR, bufferClusterIR));
                    
                    // Get the lowest energy cluster.
                    Cluster lowestEnergyCluster = pair[0];
                    Point lowestIR = pairingClusterIR;
                    if(TriggerModule.getValueClusterTotalEnergy(pair[1]) < TriggerModule.getValueClusterTotalEnergy(pair[0])) {
                        lowestEnergyCluster = pair[1];
                        lowestIR = bufferClusterIR;
                    }
                    
                    // Get deltaR for the lowest energy cluster.
                    double r = Math.sqrt(Math.pow(TriggerModule.getClusterX(lowestIR), 2) + Math.pow(TriggerModule.getClusterY(lowestIR), 2));
                    
                    // Fill the energy slope plot.
                    AIDA.defaultInstance().histogram2D(PAIR_ENERGY_SLOPE).fill(r, TriggerModule.getValueClusterTotalEnergy(lowestEnergyCluster));
                }
            }
        }
    }
    
    private void performCOPTriggerAnalysis() {
        // The COP-trigger should always use the most recent cluster
        // buffer entry.
        Collection<Cluster> clusters = clusterBuffer.getFirst();
        
        // Populate the cut distributions for all positron clusters.
        IDDecoder decoder = ReadoutDataManager.getIDDecoder("EcalHits");
        for(Cluster cluster : clusters) {
            decoder.setID(cluster.getCalorimeterHits().get(0).getCellID());
            Point ir = new Point(decoder.getValue("ix"), decoder.getValue("iy"));
            if(ir.x > 4) {
                AIDA.defaultInstance().histogram1D(COPT_MULTIPLICITY).fill(TriggerModule.getClusterHitCount(cluster));
                AIDA.defaultInstance().histogram1D(COPT_SEED_ENERGY).fill(TriggerModule.getValueClusterSeedEnergy(cluster));
                AIDA.defaultInstance().histogram1D(COPT_TOTAL_ENERGY).fill(TriggerModule.getValueClusterTotalEnergy(cluster));
                AIDA.defaultInstance().histogram2D(COPT_POSITION).fill(ir.x, ir.y);
                AIDA.defaultInstance().histogram2D(COPT_ENERGY_POSITION).fill(ir.x, TriggerModule.getValueClusterTotalEnergy(cluster));
            }
        }
    }
}