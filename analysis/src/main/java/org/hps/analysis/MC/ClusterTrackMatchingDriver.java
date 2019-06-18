package org.hps.analysis.MC;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.hps.readout.triggerstudies.PlotToTextModule;
import org.hps.recon.tracking.CoordinateTransformations;
import org.hps.recon.tracking.TrackUtils;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.FieldMap;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;

public class ClusterTrackMatchingDriver extends Driver {
    private FieldMap fieldMap = null;
    private String outputDirectory = ".";
    
    @Override
    public void detectorChanged(Detector detector) {
        fieldMap = detector.getFieldMap();
    }
    
    @Override
    public void endOfData() {
        try {
            // Cluster/track matching plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(true, true, false)),
                        new File(outputDirectory + File.separator + "tuning_clusterTrackMatching_TPR.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(true, false, false)),
                    new File(outputDirectory + File.separator + "tuning_clusterTrackMatching_TNR.dat"));
        
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(false, true, false)),
                    new File(outputDirectory + File.separator + "tuning_clusterTrackMatching_BPR.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(false, false, false)),
                    new File(outputDirectory + File.separator + "tuning_clusterTrackMatching_BNR.dat"));
            
            // Cluster/track matching debug recon plots.
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(true, true, true)),
                        new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_TPR.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(true, false, true)),
                    new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_TNR.dat"));
        
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(false, true, true)),
                    new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_BPR.dat"));
            PlotToTextModule.writePlot(AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(false, false, true)),
                    new File(outputDirectory + File.separator + "debug_reconClusterTrackMatching_BNR.dat"));
        } catch(IOException e) {
            throw new RuntimeException(e);
        }
    }
    
    @Override
    public void process(EventHeader event) {
        // Get the necessary collections.
        List<Track> gblTracks = TriggerTuningUtilityModule.getCollection(event, "GBLTracks", Track.class);
        List<Cluster> gtpClusters = TriggerTuningUtilityModule.getCollection(event, "EcalClustersGTP", Cluster.class);
        List<Cluster> reconClusters = TriggerTuningUtilityModule.getCollection(event, "EcalClustersCorr", Cluster.class);
        
        // Iterate over the tracks and plot their momenta based on
        // whether a given track is positive or negative and a top or
        // bottom track.
        for(Track gblTrack : gblTracks) {
            TrackState stateEcal = TrackUtils.getTrackExtrapAtEcalRK(gblTrack, fieldMap);
            Hep3Vector trackRVector = new BasicHep3Vector(stateEcal.getReferencePoint());
            trackRVector = CoordinateTransformations.transformVectorToDetector(trackRVector);
            double trackR[] = trackRVector.v();
            
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
                
                AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(isTop, isPositive, false)).fill(trackP, deltaR);
            }
            
            // As a debugging step, repeat this for recon clusters.
            for(Cluster reconCluster : reconClusters) {
                double[] clusterR = reconCluster.getPosition();
                double deltaX = clusterR[0] - trackR[0];
                double deltaY = clusterR[1] - trackR[1];
                double deltaR = Math.sqrt(Math.pow(deltaX, 2) + Math.pow(deltaY, 2));
                
                boolean isTop = trackR[1] > 0;
                boolean isPositive = TriggerTuningUtilityModule.isPositive(gblTrack);
                
                AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(isTop, isPositive, true)).fill(trackP, deltaR);
            }
        }
    }
    
    public void setOutputDirectory(String dir) {
        outputDirectory = dir;
    }
    
    @Override
    public void startOfData() {
        for(int i = 0; i < 2; i++) {
            boolean isTop = (i == 0);
            for(int j = 0; j < 2; j++) {
                boolean isPositive = (j == 0);
                for(int k = 0; k < 2; k++) {
                    boolean isRecon = (k == 1);
                    AIDA.defaultInstance().histogram2D(getClusterTrackMatchingPlotName(isTop, isPositive, isRecon), 250, 0.000, 5.000, 70, 0, 70);
                }
            }
        }
    }
    
    private static final String getClusterTrackMatchingPlotName(boolean isTop, boolean isPositive, boolean isRecon) {
        if(isRecon) {
            return "Debug/Cluster-Track Matching/" + (isPositive ? "Positive/" : "Negative/") + (isTop ? "Top/" : "Bottom/") + "Momentum vs. #Delta r";
        } else {
            return "Cluster-Track Matching/" + (isPositive ? "Positive/" : "Negative/") + (isTop ? "Top/" : "Bottom/") + "Momentum vs. #Delta r";
        }
    }
}
