package org.hps.readout.triggerstudies;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

public class ReadoutTriggerDistributionsDriver extends Driver {
    private static final String clusterEnergyPlot = "Verification/GTP Cluster Energy Distribution";
    private static final String clusterSeedPlot = "Verification/GTP Cluster Seed Energy Distribution";
    private static final String clusterPositionPlot = "Verification/GTP Cluster Position Distribution";
    private static final String clusterMultiplicityPlot = "Verification/GTP Cluster Multiplicity Distribution";
    
    private static final String hodoHitTruthEnergyPlot = "Verification/Hodoscope Hit Truth Energy Distribution";
    
    private static final String hodoHitEnergyPlot = "Verification/Hodoscope Hit Energy Distribution";
    private static final String hodoHitL1PositionPlot = "Verification/Hodoscope Hit Layer 1 Position Distribution";
    private static final String hodoHitL2PositionPlot = "Verification/Hodoscope Hit Layer 2 Position Distribution";
    
    private String outputDirectory = ".";
    
    @Override
    public void endOfData() {
        try {
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(clusterEnergyPlot), new File(outputDirectory + File.separator + "veri_clusterTotalEnergy.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(clusterSeedPlot), new File(outputDirectory + File.separator + "veri_clusterSeedEnergy.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(clusterPositionPlot), new File(outputDirectory + File.separator + "veri_clusterPosition.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(clusterMultiplicityPlot), new File(outputDirectory + File.separator + "veri_clusterMultiplicity.dat"));
            
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(hodoHitTruthEnergyPlot), new File(outputDirectory + File.separator + "veri_hodoHitTruthEnergy.dat"));
            
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram1D(hodoHitEnergyPlot), new File(outputDirectory + File.separator + "veri_hodoHitEnergy.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(hodoHitL1PositionPlot), new File(outputDirectory + File.separator + "veri_hodoHitLayer1Position.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(hodoHitL2PositionPlot), new File(outputDirectory + File.separator + "veri_hodoHitLayer2Position.dat"));
        } catch(IOException e) {
            throw new RuntimeException(e);
        }
        
        System.out.println("Saw " + AIDA.defaultInstance().histogram1D(clusterEnergyPlot).entries() + " clusters.");
        System.out.println("Saw " + AIDA.defaultInstance().histogram1D(hodoHitTruthEnergyPlot).entries() + " hodoscope truth hits.");
        System.out.println("Saw " + AIDA.defaultInstance().histogram1D(hodoHitEnergyPlot).entries() + " hodoscope readout hits.");
    }
    
    @Override
    public void process(EventHeader event) {
        /*
         * GTP Clusters
         */
        
        // Get the cluster collection.
        List<Cluster> gtpClusters = null;
        if(event.hasCollection(Cluster.class, "EcalClustersGTP")) {
            gtpClusters = event.get(Cluster.class, "EcalClustersGTP");
        } else {
            gtpClusters = new ArrayList<Cluster>(0);
        }
        
        // Plot the clusters.
        for(Cluster cluster : gtpClusters) {
            AIDA.defaultInstance().histogram1D(clusterEnergyPlot).fill(cluster.getEnergy());
            AIDA.defaultInstance().histogram1D(clusterSeedPlot).fill(cluster.getCalorimeterHits().get(0).getCorrectedEnergy());
            AIDA.defaultInstance().histogram1D(clusterMultiplicityPlot).fill(cluster.getCalorimeterHits().size());
            AIDA.defaultInstance().histogram2D(clusterPositionPlot).fill(cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("ix"),
                    cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy"));
        }
        
        
        /*
         * Hodoscope Truth Hits
         */
        
        // Get the hodoscope truth hit collection.
        List<SimTrackerHit> hodoscopeTruthHits = null;
        if(event.hasCollection(SimTrackerHit.class, "HodoscopeHits")) {
            hodoscopeTruthHits = event.get(SimTrackerHit.class, "HodoscopeHits");
        } else {
            hodoscopeTruthHits = new ArrayList<SimTrackerHit>(0);
        }
        
        // Plot the hodoscope hits.
        for(SimTrackerHit hit : hodoscopeTruthHits) {
            AIDA.defaultInstance().histogram1D(hodoHitTruthEnergyPlot).fill(hit.getdEdx());
        }
        
        
        /*
         * Hodoscope Hits
         */
        
        // Get the hodoscope hit collection.
        List<CalorimeterHit> hodoscopeHits = null;
        if(event.hasCollection(CalorimeterHit.class, "HodoscopeCorrectedHits")) {
            hodoscopeHits = event.get(CalorimeterHit.class, "HodoscopeCorrectedHits");
        } else {
            hodoscopeHits = new ArrayList<CalorimeterHit>(0);
        }
        
        // Plot the hodoscope hits.
        for(CalorimeterHit hit : hodoscopeHits) {
            AIDA.defaultInstance().histogram1D(hodoHitEnergyPlot).fill(hit.getCorrectedEnergy());
            if(hit.getLayerNumber() == 0) {
                AIDA.defaultInstance().histogram2D(hodoHitL1PositionPlot).fill(hit.getIdentifierFieldValue("ix"), hit.getIdentifierFieldValue("iy"));
            } else if(hit.getLayerNumber() == 1) {
                AIDA.defaultInstance().histogram2D(hodoHitL2PositionPlot).fill(hit.getIdentifierFieldValue("ix"), hit.getIdentifierFieldValue("iy"));
            } else {
                throw new RuntimeException("Unrecognized layer number \"" + hit.getLayerNumber() + "\"!");
            }
        }
    }
    
    @Override
    public void startOfData() {
        AIDA.defaultInstance().histogram1D(clusterSeedPlot, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram1D(clusterEnergyPlot, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram1D(clusterMultiplicityPlot, 10, -0.5, 9.5);
        AIDA.defaultInstance().histogram2D(clusterPositionPlot, 27, -23.5, 23.5, 11, -5.5, 5.5);
        
        AIDA.defaultInstance().histogram1D(hodoHitTruthEnergyPlot, 500, 0.000, 5.000);
        
        AIDA.defaultInstance().histogram1D(hodoHitEnergyPlot, 500, 0.000, 5.000);
        AIDA.defaultInstance().histogram2D(hodoHitL1PositionPlot, 6, -0.5, 5.5, 3, -1.5, 1.5);
        AIDA.defaultInstance().histogram2D(hodoHitL2PositionPlot, 6, -0.5, 5.5, 3, -1.5, 1.5);
    }
    
    public void setOutputDirectory(String dir) {
        outputDirectory = dir;
    }
}