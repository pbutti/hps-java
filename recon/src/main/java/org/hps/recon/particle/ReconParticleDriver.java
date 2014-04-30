package org.hps.recon.particle;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.BasicHepLorentzVector;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.HepLorentzVector;
import hep.physics.vec.VecOp;

import java.util.ArrayList;
import java.util.List;

import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Track;
import org.lcsim.event.base.BaseReconstructedParticle;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.hps.recon.ecal.HPSEcalCluster;
import org.hps.recon.tracking.TrackUtils;

/**
 * Driver that matches SVT Tracks and Ecal Clusters and creates
 * ReconstructedParticles.
 *
 * @author Mathew Graham <mgraham@slac.stanford.edu>
 * @author Omar Moreno <omoreno1@ucsc.edu>
 * @version $Id$ 
 */
public abstract class ReconParticleDriver extends Driver {

	
	// Flags
	boolean debug = false; 
	
	// Reconstructed particle collections
	List<ReconstructedParticle> finalStateParticles;
    List<ReconstructedParticle> candidates;
    List<ReconstructedParticle> candidatesBeamCon;
    List<ReconstructedParticle> candidatesTargetCon;
    List<ReconstructedParticle> electrons;
    List<ReconstructedParticle> positrons;
    
    // Collections
    String ecalClustersCollectionName 			= "EcalClusters";
    String tracksCollectionName       			= "MatchedTracks";
    String finalStateParticlesColName 		    = "FinalStateParticles";
    String candidatesCollectionName 			= null; 
    String candidatesBeamConCollectionName 		= null;
    String candidatesTargetConCollectionName	= null;
    String vertexCandidatesCollectionName       = null;
    String vertexBeamConsCandidatesName         = null;
    
    // The beamsize array is in the tracking frame
    double[] beamsize = {0.001, 0.2, 0.02};
    double maxTrackClusterDistance = 10000; // [mm] 
	double bField; 
    
     
    //  flipSign is a kludge...
    //  HelicalTrackFitter doesn't deal with B-fields in -ive Z correctly
	//  so we set the B-field in +iveZ and flip signs of fitted tracks
    //  
    //  Note:  This should be -1 for test run configurations and +1 for 
    //         prop-2014 configurations 
	int flipSign = 1;

    public ReconParticleDriver() {
    }

    public void setDebug(boolean debug){
    	this.debug = debug; 
    }
    
    public void setMaxTrackClusterDistance(double maxTrackClusterDistance){
    	this.maxTrackClusterDistance = maxTrackClusterDistance; 
    }

    public void setEcalClusterCollectionName(String ecalClustersCollectionName) {
        this.ecalClustersCollectionName = ecalClustersCollectionName;
    }

    public void setTrackCollectoinName(String tracksCollectionName) {
        this.tracksCollectionName = tracksCollectionName;
    }

    public void setBeamSigmaX(double sigma_x) {
        beamsize[1] = sigma_x; 
    }

    public void setBeamSigmaY(double sigma_y) {
        beamsize[2] = sigma_y;  
    }
    
    
	@Override
	protected void detectorChanged(Detector detector){
		
		Hep3Vector ip = new BasicHep3Vector(0., 0., 1.);
		bField = detector.getFieldMap().getField(ip).y(); 
		if(bField < 0) flipSign = -1; 
	
	}

    public void process(EventHeader event) {
       
        // All events should have a collection of Ecal clusters.  If the event 
    	// doesn't have one, skip the event.
        if (!event.hasCollection(HPSEcalCluster.class, ecalClustersCollectionName)) return;
        
        // Get the collection of Ecal clusters from the event. A triggered 
        // event should have Ecal clusters.  If it doesn't, skip the event.
        List<HPSEcalCluster> clusters = event.get(HPSEcalCluster.class, ecalClustersCollectionName);
        if(clusters.isEmpty()) return;  
        this.printDebug("Number of Ecal clusters: " + clusters.size());
        
        // Get the collection of tracks from the event
        List<Track> tracks = event.get(Track.class, tracksCollectionName);
        this.printDebug("Number of Tracks: " + tracks.size()); 

        finalStateParticles		= new ArrayList<ReconstructedParticle>();
        electrons 				= new ArrayList<ReconstructedParticle>();
        positrons 				= new ArrayList<ReconstructedParticle>();
        candidates 				= new ArrayList<ReconstructedParticle>();
        candidatesBeamCon 		= new ArrayList<ReconstructedParticle>(); 
        candidatesTargetCon 	= new ArrayList<ReconstructedParticle>(); 

        // 
        finalStateParticles = this.makeReconstructedParticles(clusters, tracks);

        // Put all the reconstructed particles in the event
        event.put(finalStateParticlesColName, finalStateParticles, ReconstructedParticle.class, 0);

        
        // Loop through the list of final state particles and separate the
        // charged particles to either electrons or positrons.  These lists
        // will be used for vertexing purposes.
        for(ReconstructedParticle finalStateParticle : finalStateParticles){
        	if(finalStateParticle.getCharge() > 0) positrons.add(finalStateParticle);
        	else if(finalStateParticle.getCharge() < 0) electrons.add(finalStateParticle);
        }
        
        // Vertex electron and positron candidates 
        vertexParticles(electrons, positrons);
       
        // If the list exist, put the vertexed candidates into the event
        if(candidatesCollectionName != null)
        	event.put(candidatesCollectionName, candidates, ReconstructedParticle.class, 0);
        if(candidatesBeamConCollectionName != null)
        	event.put(candidatesBeamConCollectionName, candidatesBeamCon, ReconstructedParticle.class, 0);
        if(candidatesTargetConCollectionName != null)
        	event.put(candidatesTargetConCollectionName, candidatesTargetCon, ReconstructedParticle.class, 0);
    }

    /**
     * 
     */
    abstract void vertexParticles(List<ReconstructedParticle> electrons, List<ReconstructedParticle> positrons);
    
    /**
     * 
     */
    private List<ReconstructedParticle> makeReconstructedParticles(List<HPSEcalCluster> clusters, List<Track> tracks){
    	
    	// Instantiate the list of reconstructed particles
    	List<ReconstructedParticle> particles = new ArrayList<ReconstructedParticle>();
    	
       	// Instantiate the list of matched tracks and clusters. This list 
       	// will be used to check if a track has already been previously matched
       	List<Track> unmatchedTracks = new ArrayList<Track>(tracks); 
    	
  	    Track matchedTrack = null;
  	    // Loop over all of the Ecal clusters and pair them with tracks
    	for(HPSEcalCluster cluster : clusters){
    	
    		// Get the position of the Ecal cluster
    		Hep3Vector clusterPosition = new BasicHep3Vector(cluster.getPosition());
    		this.printDebug("Ecal cluster position: " + clusterPosition.toString());
    		
            double rMax = Double.MAX_VALUE;
    		// Loop over all of the tracks in the event
    		for(Track track : tracks){
    		
    			// Check if the Ecal cluster and track are within the same 
    			// detector volume i.e. both top or bottom
    			if(clusterPosition.y()*track.getTrackStates().get(0).getZ0() < 0) continue; 
    			
    			// Extrapolate the track to the Ecal cluster position
    			Hep3Vector trackPosAtEcal = TrackUtils.extrapolateTrack(track, clusterPosition.z());
    			
       			// Check if any of the extrapolated values are invalid.
       			// TODO: There are some track whose extrapolated coordinates
       			//        are NaN. The problem seems to be that the y-coordinate
       			//        of the extrapolated helix is found to be non-real. This
       			//        needs to be fixed.
        		if(Double.isNaN(trackPosAtEcal.x()) || Double.isNaN(trackPosAtEcal.y())) continue; 
        		this.printDebug("Track position at shower max: " + trackPosAtEcal.toString());
        		
       			double r = VecOp.sub(trackPosAtEcal, clusterPosition).magnitude();
       			this.printDebug("Distance between Ecal cluster and track position: " + r + " mm");
        		
                if (r < rMax && r <= maxTrackClusterDistance) {
                	rMax = r;
                	matchedTrack = track; 
                }
    		}
    		
    		// Create a reconstructed particle and add it to the 
    		// collection of particles
    		ReconstructedParticle particle = new BaseReconstructedParticle(); 
    		HepLorentzVector fourVector = new BasicHepLorentzVector(0, 0, 0, 0); 
    		particle.addCluster(cluster); 
    		((BasicHepLorentzVector) fourVector).setT(cluster.getEnergy());
    		if(matchedTrack != null){
    			particle.addTrack(matchedTrack);
    			Hep3Vector momentum = new BasicHep3Vector(matchedTrack.getTrackStates().get(0).getMomentum());
    			((BasicHepLorentzVector) fourVector).setV3(fourVector.t(), momentum);
    			((BaseReconstructedParticle) particle).setCharge(matchedTrack.getCharge()*flipSign);
    		
    			if(unmatchedTracks.contains(matchedTrack)) unmatchedTracks.remove(matchedTrack);
    		}
    		((BaseReconstructedParticle) particle).set4Vector(fourVector);
    	
    		particles.add(particle);
    	}
    	
    	if(!unmatchedTracks.isEmpty()){
    		for(Track  track : unmatchedTracks){
    			// Create a reconstructed particle and add it to the 
    			// collection of particles
    			ReconstructedParticle particle = new BaseReconstructedParticle(); 
    			HepLorentzVector fourVector = new BasicHepLorentzVector(0, 0, 0, 0); 
    			
    			particle.addTrack(matchedTrack);
    			Hep3Vector momentum = new BasicHep3Vector(track.getTrackStates().get(0).getMomentum());
    			((BasicHepLorentzVector) fourVector).setV3(fourVector.t(), momentum);
    			((BaseReconstructedParticle) particle).setCharge(track.getCharge()*flipSign);
    			((BaseReconstructedParticle) particle).set4Vector(fourVector);
    			
    			particles.add(particle);
    		}
    	}
        
    	return particles;
    }
       
    /**
     * 
     * @param debugMessage
     */
    private void printDebug(String debugMessage){
    	if(debug)
    		System.out.println(this.getClass().getSimpleName() + ": " + debugMessage); 
    }
}
