package org.hps.readout.triggerstudies;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import org.hps.readout.ReadoutDataManager;
import org.hps.readout.ReadoutDriver;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.util.aida.AIDA;

import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;

/**
 * Driver <code>ClusterEnergyDistributionDriver</code> is designed to
 * plot the distribution of cluster energies during readout.
 * 
 * @author Kyle McCarty <mccarty@jlab.org>
 */
public class ClusterEnergyDistributionDriver extends ReadoutDriver {
    private double localTime = 0.0;
    private String outputFilepath = "";
    private String clusterCollectionName = "EcalClustersGTP";
    private final static String ENERGY_PLOT_NAME = "Trigger Studies/Pre-Readout/All Cluster Total Energy Distribution";
    private final static String SEED_PLOT_NAME = "Trigger Studies/Pre-Readout/All Cluster Seed Energy Distribution";
    private final static String SEED_VS_N_PLOT_NAME = "Trigger Studies/Pre-Readout/All Cluster Multiplicity vs. Seed Energy Distribution";
    private final static String SEED_VS_ENERGY_PLOT_NAME = "Trigger Studies/Pre-Readout/All Cluster Total Energy vs. Energy Distribution";
    
    private final static String CLUSTER_TOTAL_PLOT_SUFFIX = "-clusterEnergyDistro.dat";
    private final static String CLUSTER_SEED_PLOT_SUFFIX = "-clusterSeedDistro.dat";
    private final static String CLUSTER_SEED_VS_N_PLOT_SUFFIX = "-seedVsNDistro.dat";
    private final static String CLUSTER_SEED_VS_TOTAL_PLOT_SUFFIX = "-seedVsEnergyDistro.dat";
    
    @Override
    public void endOfData() {
        // Write the plot as a Mathematica-compatible file.
        try {
            final String[] plotNames = { ENERGY_PLOT_NAME, SEED_PLOT_NAME, SEED_VS_N_PLOT_NAME, SEED_VS_ENERGY_PLOT_NAME };
            final String[] plotSuffixes = { CLUSTER_TOTAL_PLOT_SUFFIX, CLUSTER_SEED_PLOT_SUFFIX, CLUSTER_SEED_VS_N_PLOT_SUFFIX, CLUSTER_SEED_VS_TOTAL_PLOT_SUFFIX };
            final Class<?>[] type = { IHistogram1D.class, IHistogram1D.class, IHistogram2D.class, IHistogram2D.class};
            
            FileWriter writer = null;
            String outputText = null;
            for(int i = 0; i < plotNames.length; i++) {
                writer = new FileWriter(outputFilepath + plotSuffixes[i]);
                if(type[i] == IHistogram1D.class) {
                    outputText = PlotToTextModule.aidaToMathematica(AIDA.defaultInstance().histogram1D(plotNames[i]));
                } else {
                    outputText = PlotToTextModule.aidaToMathematica(AIDA.defaultInstance().histogram2D(plotNames[i]));
                }
                writer.write(outputText);
                writer.close();
            }
        } catch(IOException e) {
            throw new RuntimeException(e);
        }
    }
    
    @Override
    public void process(EventHeader event) {
        // Do nothing until the clusters are ready.
        if(!ReadoutDataManager.checkCollectionStatus(clusterCollectionName, localTime + 4.0)) {
            return;
        }
        
        // Get the clusters from the data manager.
        Collection<Cluster> clusters = ReadoutDataManager.getData(localTime, localTime + 4.0, clusterCollectionName, Cluster.class);
        
        // Plot the cluster energies.
        for(Cluster cluster : clusters) {
            double seedEnergy = cluster.getCalorimeterHits().get(0).getCorrectedEnergy();
            AIDA.defaultInstance().histogram1D(SEED_PLOT_NAME).fill(seedEnergy);
            AIDA.defaultInstance().histogram1D(ENERGY_PLOT_NAME).fill(cluster.getEnergy());
            AIDA.defaultInstance().histogram2D(SEED_VS_ENERGY_PLOT_NAME).fill(seedEnergy, cluster.getEnergy());
            AIDA.defaultInstance().histogram2D(SEED_VS_N_PLOT_NAME).fill(seedEnergy, cluster.getCalorimeterHits().size());
        }
        
        // Increment the local time.
        localTime += 4.0;
    }
    
    /**
     * Sets the name of the cluster collection that is used to
     * populate the output plot.
     * @param collectionName - The cluster collection name.
     */
    public void setClusterCollectionName(String collectionName) {
        clusterCollectionName = collectionName;
    }
    
    /**
     * Defines the root name of the output files. All files will use
     * this as a prefix.
     * @param filepath - The root name of the output files.
     */
    public void setOutputFile(String filepath) {
        outputFilepath = filepath;
    }
    
    @Override
    public void startOfData() {
        // Add the driver dependencies.
        addDependency(clusterCollectionName);
        
        // Instantiate the cluster plot.
        AIDA.defaultInstance().histogram1D(SEED_PLOT_NAME, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram1D(ENERGY_PLOT_NAME, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram2D(SEED_VS_N_PLOT_NAME, 500, 0.000, 5.000, 10, -0.5, 9.5);
        AIDA.defaultInstance().histogram2D(SEED_VS_ENERGY_PLOT_NAME, 500, 0.000, 5.000, 500, 0.000, 5.000);
    }
    
    @Override
    protected double getTimeDisplacement() {
        return 0;
    }
    
    @Override
    protected double getTimeNeededForLocalOutput() {
        return 0;
    }
}