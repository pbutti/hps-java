package org.hps.readout.triggerstudies;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import org.hps.readout.ReadoutDataManager;
import org.hps.readout.ReadoutDriver;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.util.aida.AIDA;

import hep.aida.IHistogram1D;

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
    private final static String PLOT_NAME = "Trigger Studies/Pre-Readout/All Cluster Energy Distribution";
    
    @Override
    public void endOfData() {
        // Write the plot as a Mathematica-compatible file.
        try {
            FileWriter writer = new FileWriter(outputFilepath);
            String outputText = aidaToMathematica(AIDA.defaultInstance().histogram1D(PLOT_NAME));
            writer.write(outputText);
            writer.close();
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
            AIDA.defaultInstance().histogram1D(PLOT_NAME).fill(cluster.getEnergy());
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
     * Defines the output file where the plot data will be saved.
     * @param filepath - The filepath to where the plot should be
     * saved.
     */
    public void setOutputFilepath(String filepath) {
        outputFilepath = filepath;
    }
    
    @Override
    public void startOfData() {
        // Add the driver dependencies.
        addDependency(clusterCollectionName);
        
        // Instantiate the cluster plot.
        AIDA.defaultInstance().histogram1D(PLOT_NAME, 500, 0.000, 5.000);
        
        // Validate that an appropriate output filepath was given.
        File outputFile = new File(outputFilepath);
        if(outputFile.isDirectory()) {
            throw new IllegalArgumentException("Output file can not be a directory.");
        } else if(outputFile.getParentFile() == null || outputFile.getParentFile().exists()) {
            if(outputFile.exists()) { outputFile.delete(); }
        } else {
            throw new IllegalArgumentException("Output file directory not found!");
        }
    }

    @Override
    protected double getTimeDisplacement() {
        return 0;
    }
    
    @Override
    protected double getTimeNeededForLocalOutput() {
        return 0;
    }
    
    /**
     * Converts an AIDA histogram to a text file formatted for import
     * into Wolfram Mathematica as a list of points. Note that this
     * method does not handle writing the output file.
     * @param plot - The plot to convert.
     * @return Returns a {@link java.util.String} object that can be
     * read into Mathematica as a list of points.
     */
    private static final String aidaToMathematica(IHistogram1D plot) {
        StringBuffer outputBuffer = new StringBuffer();
        
        int bins = plot.axis().bins();
        
        for(int x = 0; x < bins; x++) {
            outputBuffer.append(String.format("%f,%f%n",plot.axis().binCenter(x), plot.binHeight(x)));
        }
        
        return outputBuffer.toString();
    }
}