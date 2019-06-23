package org.hps.analysis.MC;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.hps.conditions.hodoscope.HodoscopeChannel;
import org.hps.detector.hodoscope.HodoscopeDetectorElement;
import org.hps.readout.triggerstudies.Coordinate;
import org.hps.recon.tracking.TrackUtils;
import org.hps.record.triggerbank.TriggerModule;
import org.hps.util.Pair;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.fit.helicaltrack.HelicalTrackFit;
import org.lcsim.geometry.FieldMap;

import hep.aida.IHistogram2D;
import hep.physics.vec.BasicHep3Vector;

public class TriggerTuningUtilityModule {
    
    public static final List<Pair<Cluster, Track>> getClusterTrackMatchedPairs(List<Cluster> clusters, List<Track> tracks, FieldMap fieldMap,
            double[][] boundsTopPos, double[][] boundsBotPos, double[][] boundsTopEle, double[][] boundsBotEle, double[] boundsPositronP, double[] boundsElectronP) {
        // Clusters may only be matched to one track. Track which
        // clusters have already been matched.
        Set<Cluster> matchedClusters = new HashSet<Cluster>();
        
        // Store matched pairs.
        List<Pair<Cluster, Track>> pairList = new ArrayList<Pair<Cluster, Track>>();
        
        // Iterate over the tracks.
        trackLoop:
        for(Track track : tracks) {
            // Matching criteria depends on the charge of the track
            // and whether it is a top or bottom track.
            boolean isPositive = isPositive(track);
            boolean isTop = isTopTrack(track);
            
            // Get the track momentum and check that it is within the
            // allowed momentum range for its charge. Tracks outside
            // this range are very unlikely to be A' tracks.
            double trackP = getMomentumMagnitude(track, fieldMap);
            if(isPositive && (trackP < boundsPositronP[0] || trackP > boundsPositronP[1])) { continue trackLoop; }
            else if(!isPositive && (trackP < boundsElectronP[0] || trackP > boundsElectronP[1])) { continue trackLoop; }
            
            // Get the track position at the calorimeter face.
            double[] trackR = getTrackPositionAtCalorimeterFace(track, fieldMap);
            
            // Get the correct set of matching parameters for this
            // particular track.
            double[][] bounds = null;
            if(isPositive && isTop)        { bounds = boundsTopPos; }
            else if(isPositive && !isTop)  { bounds = boundsBotPos; }
            else if(!isPositive && isTop)  { bounds = boundsTopEle; }
            else if(!isPositive && !isTop) { bounds = boundsBotEle; }
            
            // Check each cluster to see if it matches the track. It
            // is assumed, given pure signal data, that there will
            // only be one plausible match, so the first is accepted.
            clusterLoop:
            for(Cluster cluster : clusters) {
                if(matches(cluster, trackP, trackR, bounds[0], bounds[1])) {
                    matchedClusters.add(cluster);
                    pairList.add(new Pair<Cluster, Track>(cluster, track));
                    break clusterLoop;
                }
            }
        }
        
        // Return the matched pairs.
        return pairList;
    }
    
    /**
     * Gets a the collection containing objects of the specified type
     * from the event.
     * @param event - The event.
     * @param collectionName - The name of the collection.
     * @param type - The object type of the collection.
     * @return Returns the collection if it exists, and an empty list
     * if it does not.
     */
    public static final <T> List<T> getCollection(EventHeader event, String collectionName, Class<T> type) {
        if(event.hasCollection(type, collectionName)) {
            return event.get(type, collectionName);
        } else {
            return new ArrayList<T>(0);
        }
    }
    
    /**
     * Adds the energies for all hodoscope hits that occur on the
     * same scintillator together. Energies are stored in an array
     * that follows the format <code>double[Scintillator Number -
     * 1][Top = 0; Bottom = 1][L1 = 0; L2 = 1]</code>.
     * @param hodoscopeHits - The list of hits representing energy
     * depositions on an FADC channel.
     * @param hodoscopeChannelMap - A map which maps a channel ID to
     * its {@link org.hps.conditions.hodoscope.HodoscopeChannel
     * HodoscopeChannel} object.
     * @return Returns the compiled hodoscope hit energies as an
     * array.
     */
    public static final double[][][] getCompiledHodoscopeEnergies(List<? extends CalorimeterHit> hodoscopeHits, Map<Long, HodoscopeChannel> hodoscopeChannelMap) {
        // Instantiate the output array. The first index corresponds
        // to the scintillator number, the second to top/bottom, and
        // the third to the layer.
        // double[Scintillator; 0 - 4][Top = 0; Bottom = 1][L1 = 0; L2 = 1]
        double[][][] energies = new double[5][2][2];
        
        // Hodoscope hits that occur on the same scintillator add
        // their energies together.
        for(CalorimeterHit hit : hodoscopeHits) {
            // Get the array index versions of the scintillator
            // position.
            HodoscopeChannel channel = hodoscopeChannelMap.get(Long.valueOf(hit.getCellID()));
            int arrayX = channel.getIX();
            int arrayY = (channel.getIY() == HodoscopeChannel.TOP ? 0 : 1);
            int arrayZ = (channel.getLayer() == HodoscopeChannel.LAYER_1 ? 0 : 1);
            
            // Add the hit energy to appropriate array index.
            energies[arrayX][arrayY][arrayZ] += hit.getCorrectedEnergy();
        }
        
        // Return the results.
        return energies;
    }
    
    /**
     * Adds the energies for all hodoscope hits that occur on the
     * same scintillator together. Energies are stored in an array
     * that follows the format <code>double[Scintillator Number -
     * 1][Top = 0; Bottom = 1][L1 = 0; L2 = 1]</code>.
     * @param hodoscopeHits - The list of hits representing energy
     * depositions on a SLIC scintillator.
     * @param hodoscopeDetectorElement - The detector element for
     * the hodoscope.
     * @return Returns the compiled hodoscope hit energies as an
     * array.
     */
    public static final double[][][] getCompiledHodoscopeEnergies(List<SimTrackerHit> hodoscopeHits, HodoscopeDetectorElement hodoscopeDetectorElement) {
        // Instantiate the output array. The first index corresponds
        // to the scintillator number, the second to top/bottom, and
        // the third to the layer.
        // double[Scintillator; 0 - 4][Top = 0; Bottom = 1][L1 = 0; L2 = 1]
        double[][][] energies = new double[5][2][2];
        
        // Hodoscope hits that occur on the same scintillator add
        // their energies together.
        for(SimTrackerHit hit : hodoscopeHits) {
            // Get the array index versions of the scintillator
            // position.
            int[] arrayIndices = hodoscopeDetectorElement.getHodoscopeIndices(hit);
            int arrayX = arrayIndices[0];
            int arrayY = (arrayIndices[1] == HodoscopeChannel.TOP ? 0 : 1);
            int arrayZ = arrayIndices[2];
            
            // Add the hit energy to appropriate array index.
            energies[arrayX][arrayY][arrayZ] += hit.getdEdx();
        }
        
        // Return the results.
        return energies;
    }
    
    public static final Coordinate[] getEnergyThreshold(IHistogram2D plot, double xMin, double yMin, double yMax, double threshold, boolean fromTop) {
        // Get the number of bins in each direction.
        int xBins = plot.xAxis().bins();
        int yBins = plot.yAxis().bins();
        
        // Get the x-axis values.
        double[] xVals = new double[xBins];
        for(int i = 0; i < xBins; i++) {
            xVals[i] = plot.xAxis().binCenter(i);
        }
        
        double[] yVals = new double[yBins];
        for(int i = 0; i < yBins; i++) {
            yVals[i] = plot.yAxis().binCenter(i);
        }
        
        // Determine the total number of counts for all y-bins for a
        // given x-bin.
        double[] totalCount = new double[xBins];
        for(int x = 0; x < xBins; x++) {
            yLoop:
            for(int y = 0; y < yBins; y++) {
                // Ignore y-values that are below the minimum
                // allowed y-value or above the maximum y-value.
                if(yVals[y] < yMin || yVals[y] > yMax) { continue yLoop; }
                
                // Increment the total count for this x.
                totalCount[x] += plot.binHeight(x, y);
            }
        }
        
        // Determine which y-bin is below [threshold] percent of all
        // entries for the current x-value. Ignore x-values below
        // [xMin].
        List<Coordinate> thresholdList = new ArrayList<Coordinate>();
        
        // If the threshold is cumulative from the top...
        if(fromTop) {
            xLoop:
            for(int x = 0; x < xBins; x++) {
                // Ignore x-values that are below the minimum specified
                // x-value.
                if(xVals[x] < xMin) { continue xLoop; }
                
                // Otherwise, start from the top y-bin and check whether
                // the cumulative total entries from that bin to the
                // current bin is equal to or greater than the threshold.
                int yTotal = 0;
                yLoop:
                for(int y = yBins - 1; y >= 0; y--) {
                    // Ignore y-values that are below the minimum
                    // allowed y-value or above the maximum y-value.
                    if(yVals[y] < yMin || yVals[y] > yMax) { continue yLoop; }
                    
                    // Increment the cumulative y-entries.
                    yTotal += plot.binHeight(x, y);
                    
                    // Calculate the cumulative percentage.
                    double percentage = yTotal / totalCount[x];
                    
                    // If it exceeds the threshold, store it.
                    if(percentage >= threshold) {
                        thresholdList.add(new Coordinate(xVals[x], yVals[y]));
                        break yLoop;
                    }
                }
            }
        } else {
            xLoop:
            for(int x = 0; x < xBins; x++) {
                // Ignore x-values that are below the minimum specified
                // x-value.
                if(xVals[x] < xMin) { continue xLoop; }
                
                // Otherwise, start from the top y-bin and check whether
                // the cumulative total entries from that bin to the
                // current bin is equal to or greater than the threshold.
                int yTotal = 0;
                yLoop:
                for(int y = 0; y < yBins; y++) {
                    // Ignore y-values that are below the minimum
                    // allowed y-value or above the maximum y-value.
                    if(yVals[y] < yMin || yVals[y] > yMax) { continue yLoop; }
                    
                    // Increment the cumulative y-entries.
                    yTotal += plot.binHeight(x, y);
                    
                    // Calculate the cumulative percentage.
                    double percentage = yTotal / totalCount[x];
                    
                    // If it exceeds the threshold, store it.
                    if(percentage >= threshold) {
                        thresholdList.add(new Coordinate(xVals[x], yVals[y]));
                        break yLoop;
                    }
                }
            }
        }
        
        // Return the list of threshold coordinates.
        return thresholdList.toArray(new Coordinate[thresholdList.size()]);
    }
    
    /**
     * Gets the highest energy cluster of a set of two. If one
     * cluster is <code>null</code>, it always considered the
     * lowest-energy cluster. If both are <code>null</code>, a value
     * of <code>null</code> is returned.
     * @param cluster1 - The first cluster.
     * @param cluster2 - The second cluster.
     * @return Returns <code>null</code> if both clusters are not
     * defined, the defined clusters if only one cluster is not
     * defined, or the cluster with the highest energy if both
     * clusters are defined.
     */
    public static final Cluster getHighestEnergyCluster(Cluster cluster1, Cluster cluster2) {
        if(cluster1 == null) { return cluster2; }
        else if(cluster2 == null) { return cluster1; }
        else if(cluster1.getEnergy() > cluster2.getEnergy()) { return cluster1; }
        else { return cluster2; }
    }
    
    /**
     * Creates an array of all possible hodoscope cluster energies.
     * The array follows the form
     * <code>hodoscopeClusterEnergies[ix][iy][iz]</code> where
     * <code>ix</code> maps to the following:
     * <ul>
     * <li><code>ix = 0</code> --> S1</li>
     * <li><code>ix = 1</code> --> S1 + S2</li>
     * <li><code>ix = 2</code> --> S2</li>
     * <li><code>ix = 3</code> --> S2 + S3</li>
     * <li><code>ix = 4</code> --> S3</li>
     * <li><code>ix = 5</code> --> S3 + S4</li>
     * <li><code>ix = 6</code> --> S4</li>
     * <li><code>ix = 7</code> --> S4 + S5</li>
     * <li><code>ix = 8</code> --> S5</li>
     * </ul>
     * @param hodoscopeHitEnergies - The compiled hodoscope energies
     * for each individual scintillator. This may be obtained from
     * event objects through the method {@link
     * org.hps.analysis.MC.TriggerTuningUtilityModule#getCompiledHodoscopeEnergies(List, Map)
     * getCompiledHodoscopeEnergies(List<CalorimeterHit>,
     * Map<Long, HodoscopeChannel>)}.
     * @return Returns an array containing all possible hodoscope
     * cluster energies.
     */
    public static final double[][][] getHodoscopeClusterEnergies(double[][][] hodoscopeHitEnergies) {
        // Create all possible hodoscope clusters. These consist of
        // the energies of each individual scintillator, as well as
        // the energy sum of all adjacent scintillators, if both have
        // a defined energy individually.
        double[][][] hodoscopeClusterEnergies = new double[9][2][2];
        for(int iy = 0; iy < 2; iy++) {
            for(int iz = 0; iz < 2; iz++) {
                hodoscopeClusterEnergies[0][iy][iz] = hodoscopeHitEnergies[0][iy][iz];
                if(hodoscopeHitEnergies[0][iy][iz] != 0 && hodoscopeHitEnergies[1][iy][iz] != 0) {
                    hodoscopeClusterEnergies[1][iy][iz] = hodoscopeHitEnergies[0][iy][iz] + hodoscopeHitEnergies[1][iy][iz];
                }
                hodoscopeClusterEnergies[2][iy][iz] = hodoscopeHitEnergies[1][iy][iz];
                if(hodoscopeHitEnergies[1][iy][iz] != 0 && hodoscopeHitEnergies[2][iy][iz] != 0) {
                    hodoscopeClusterEnergies[3][iy][iz] = hodoscopeHitEnergies[1][iy][iz] + hodoscopeHitEnergies[2][iy][iz];
                }
                hodoscopeClusterEnergies[4][iy][iz] = hodoscopeHitEnergies[2][iy][iz];
                if(hodoscopeHitEnergies[2][iy][iz] != 0 && hodoscopeHitEnergies[3][iy][iz] != 0) {
                    hodoscopeClusterEnergies[5][iy][iz] = hodoscopeHitEnergies[2][iy][iz] + hodoscopeHitEnergies[3][iy][iz];
                }
                hodoscopeClusterEnergies[6][iy][iz] = hodoscopeHitEnergies[3][iy][iz];
                if(hodoscopeHitEnergies[3][iy][iz] != 0 && hodoscopeHitEnergies[4][iy][iz] != 0) {
                    hodoscopeClusterEnergies[7][iy][iz] = hodoscopeHitEnergies[3][iy][iz] + hodoscopeHitEnergies[4][iy][iz];
                }
                hodoscopeClusterEnergies[8][iy][iz] = hodoscopeHitEnergies[4][iy][iz];
            }
        }
        
        // Return the results.
        return hodoscopeClusterEnergies;
    }
    
    /**
     * Gets the invariant mass of the originating particle of a
     * positron/electron two particle decay.
     * @param positron - The positron track.
     * @param electron - The electron track.
     * @param fieldMap - The magnetic fieldmap in which the particles
     * are moving.
     * @return Returns the invariant mass.
     */
    public static final double getInvarientMass(Track positron, Track electron, FieldMap fieldMap) {
        // Define the mass of the particles.
        final double m = 0.51099895000 / 1000;
        
        // Get the momenta of the particles.
        double[][] p = new double[][] { getMomentum(positron, fieldMap), getMomentum(electron, fieldMap) };
        double[] magP = new double[] { getMomentumMagnitude(positron, fieldMap), getMomentumMagnitude(electron, fieldMap) };
        
        // Get the Lorentz factor for each particle.
        double[] gamma = new double[] { getLorentzFactor(magP[0], m), getLorentzFactor(magP[1], m) };
        
        // Calculate the originating particle invariant mass.
        return Math.sqrt(Math.pow(gamma[0] * m + gamma[1] * m, 2) - Math.pow(getMagnitude(getVectorSum(p[0], p[1])), 2));
    }
    
    /**
     * Calculates the Lorentz factor for a particle with mass m and
     * with momentum p.
     * @param p - The particle momentum.
     * @param m - The particle rest mass.
     * @return Returns the Lortentz factor.
     */
    public static final double getLorentzFactor(double p, double m) {
        // Calculate the lorentz factor.
        return Math.sqrt(1 + Math.pow(p / m, 2));
    }
    
    /**
     * Gets the magnitude of a vector.
     * @param v - The vector.
     * @return Returns the magnitude of the vector.
     */
    public static final double getMagnitude(double[] v) {
        double squareSum = 0.0;
        for(double vi : v) {
            squareSum += Math.pow(vi, 2);
        }
        return Math.sqrt(squareSum);
    }
    
    /**
     * Gets the vector momentum of a track.
     * @param track - The track.
     * @return Returns the track momentum as a size three array of
     * type <code>double</code>.
     */
    @Deprecated
    public static final double[] getMomentum(Track track, FieldMap fieldMap) {
        double phi = TrackUtils.getTrackStateAtLocation(track, TrackState.AtIP).getPhi();
        double tanLambda = TrackUtils.getTrackStateAtLocation(track, TrackState.AtIP).getTanLambda();
        //double phi = track.getTrackStates().get(0).getPhi();
        //double tanLambda = track.getTrackStates().get(0).getTanLambda();
        
        double magP = getMomentumMagnitude(track, fieldMap);
        
        double px = magP * Math.cos(phi);
        double py = magP * Math.sin(phi);
        double pz = magP * tanLambda;
        
        return new double[] { px, py, pz };
    }
    
    /**
     * Gets the magnitude of the momentum of a track.
     * @param track - The track.
     * @return Returns the track momentum as a <code>double</code>.
     */
    public static final double getMomentumMagnitude(Track track, FieldMap fieldMap) {
        // Matt Solt code.
        HelicalTrackFit helicalTrackFit = TrackUtils.getHTF(track);
        double bFieldY = fieldMap.getField(new BasicHep3Vector(0, 0, 500)).y();
        return Math.abs(helicalTrackFit.p(bFieldY));
    }
    
    /**
     * Gets the track phi at the vertex.
     * @param track - The track.
     * @return Returns the track phi angle at the vertex.
     */
    public static final double getTrackPhi(Track track) {
        return TrackUtils.getTrackStateAtLocation(track, TrackState.AtIP).getPhi();
    }
    
    /**
     * Extrapolates the track position at the calorimeter face.
     * <br/><br/>
     * <b>NOTE:</b> This uses a deprecated method and needs to be
     * replaced when a correct version is known.
     * @param track - The track.
     * @return Returns the extrapolated track position at the face of
     * the calorimeter as a size 3 <code>double</code> array.
     */
    public static final double[] getTrackPositionAtCalorimeterFace(Track track, FieldMap fieldMap) {
        double[] tempP = TrackUtils.getTrackExtrapAtEcalRK(track, fieldMap).getReferencePoint();
        //double[] tempP = TrackUtils.getTrackStateAtECal(track).getReferencePoint();
        return new double[] { tempP[1], tempP[2], tempP[0] };
    }
    
    /**
     * Gets the component sum of two vectors.
     * @param v0 - A vector.
     * @param v1 - Another vector.
     * @return Returns the vector sum of the vectors.
     * @throws IllegalArgumentException Occurs if the vectors are of
     * different lengths.
     */
    public static final double[] getVectorSum(double[] v0, double[] v1) throws IllegalArgumentException {
        // Vectors of different lengths may not be summed.
        if(v0.length != v1.length) {
            throw new IllegalArgumentException("Vectors can only be added if they are the same length.");
        }
        
        // Sum the vectors component by component.
        double[] vSum = new double[v0.length];
        for(int i = 0; i < v0.length; i++) {
            vSum[i] = v0[i] + v1[i];
        }
        
        // Return the resultant sum.
        return vSum;
    }
    
    /**
     * Specifies whether or not a cluster is a bottom cluster.
     * @param cluster - The cluster.
     * @return Returns <code>true</code> for bottom clusters and
     * <code>false</code> for top clusters.
     */
    public static final boolean isBottomCluster(Cluster cluster) {
        return !isTopCluster(cluster);
    }
    
    /**
     * Specifies whether or not a track is a bottom track. This is
     * based on the tan(Lambda) of the track.
     * @param track - The track.
     * @return Returns <code>true</code> for bottom tracks and
     * <code>false</code> for top tracks.
     */
    public static final boolean isBottomTrack(Track track) {
        return !isTopTrack(track);
    }
    
    /**
     * Specifies whether an event meets the minimum requirements for
     * a calorimeter-only positron event. This requires the event to
     * have a cluster at a position of <code>x >= xLowerBound</code>.
     * @param event - The event to test.
     * @param gtpClusterCollectionName - The name of the collection
     * containing the clusters.
     * @return Returns <code>true</code> if there exists a cluster
     * pair that meets the COPT conditions and otherwise returns
     * <code>false</code>.
     */
    public static final boolean isCOPTEvent(EventHeader event, String clusterCollectionName, double xLowerBound) {
        // Get the cluster collection.
        List<Cluster> clusters = getCollection(event, clusterCollectionName, Cluster.class);
        
        // Track whether there exists a cluster on the positron side
        // of the calorimeter. If there is, this event passes.
        for(Cluster cluster : clusters) {
            if(TriggerModule.getClusterX(cluster) >= xLowerBound) {
                return true;
            }
        }
        
        // If no valid cluster was found, the event fails.
        return false;
    }
    
    /**
     * Specifies whether a cluster is located on the electron side of
     * the calorimeter.
     * @param cluster - The cluster.
     * @return Returns <code>true</code> if the cluster is on the
     * electron side of the calorimeter and <code>false</code> if it
     * is not.
     */
    public static final boolean isElectronSideCluster(Cluster cluster) {
        return !isPositronSideCluster(cluster);
    }
    
    /**
     * Specifies whether or not an event meets the conditions for use
     * in data analysis. This requires that the event has a positive
     * and negative track, one of which must be a top track and one
     * which must be a bottom track. The tracks are also required to
     * pass a chi squared upper bound cut.
     * @param event - The event to test.
     * @param chiSquaredUpperBound - The chi squared upper bound cut.
     * @return Returns <code>true</code> if there exists a track pair
     * that meets the conditions for analyzability and otherwise
     * returns <code>false</code>.
     */
    public static final boolean isGoodEvent(EventHeader event, String gblTrackCollectionName, double chiSquaredUpperBound) {
        // Get all GBL tracks.
        List<Track> gblTracks = getCollection(event, gblTrackCollectionName, Track.class);
        
        // Process the tracks. This should check for all permutations
        // of charge and position.
        boolean sawTopPositiveTrack = false;
        boolean sawBotPositiveTrack = false;
        boolean sawTopNegativeTrack = false;
        boolean sawBotNegativeTrack = false;
        for(Track track : gblTracks) {
            // Permutations of tracks are only considered if they
            // pass the chi squared cut.
            if(track.getChi2() <= chiSquaredUpperBound) {
                if(isPositive(track)) {
                    if(isTopTrack(track)) { sawTopPositiveTrack = true; }
                    else { sawBotPositiveTrack = true; }
                } else {
                    if(isTopTrack(track)) { sawTopNegativeTrack = true; }
                    else { sawBotNegativeTrack = true; }
                }
            }
        }
        
        // A good event has a positive and negative track, where one
        // is a top track and the other a bottom track.
        if((sawTopPositiveTrack && sawBotNegativeTrack) || (sawTopNegativeTrack && sawBotPositiveTrack)) {
            return true;
        } else { return false; }
    }
    
    /**
     * Specifies whether an event meets the minimum requirements for
     * a hodoscope event. This simply requires that a hodoscope hit
     * actually exist and that there exists a cluster in the allowed
     * positron region.
     * @param event - The event to test.
     * @param hodoscopeHitCollectionName - The name of the collection
     * containing the hodoscope hits.
     * @param xLowerBound - The lower bound that defines the allowed
     * positron region.
     * @return Returns <code>true</code> if there is a hodoscope hit
     * and returns <code>false</code> otherwise.
     */
    public static final boolean isHodoscopeEvent(EventHeader event, String hodoscopeHitCollectionName, String clusterCollectionName, double xLowerBound) {
        // Check that there exists a cluster in the positron region.
        if(!isCOPTEvent(event, clusterCollectionName, xLowerBound)) { return false; }
        
        // Get the cluster collection.
        List<CalorimeterHit> hodoHits = getCollection(event, hodoscopeHitCollectionName, CalorimeterHit.class);
        
        // A hodoscope is any event that has a hodoscope hit.
        return !hodoHits.isEmpty();
    }
    
    /**
     * Indicates whether or not a track is negatively charged.
     * @param track - The track.
     * @return Returns <code>true</code> true if the track is
     * negatively charged and <code>false</code> otherwise.
     */
    public static final boolean isNegative(Track track) {
        return (-Math.signum(track.getTrackStates().get(0).getOmega())) < 0;
    }
    
    /**
     * Specifies whether an event meets the minimum requirements for
     * a pair event. These requirements are that an event has a pair
     * of clusters where one cluster is on the top of the calorimeter
     * and one cluster is on the bottom.
     * @param event - The event to test.
     * @param gtpClusterCollectionName - The name of the collection
     * containing the clusters.
     * @return Returns <code>true</code> if there exists a cluster
     * pair that meets the pair conditions and otherwise returns
     * <code>false</code>.
     */
    public static final boolean isPairEvent(EventHeader event, String clusterCollectionName) {
        // Get the cluster collection.
        List<Cluster> clusters = getCollection(event, clusterCollectionName, Cluster.class);
        
        // Track whether there exists any top and bottom cluster.
        boolean sawTopCluster = false;
        boolean sawBotCluster = false;
        for(Cluster cluster : clusters) {
            if(TriggerModule.getClusterYIndex(cluster) > 0) {
                sawTopCluster = true;
            } else { sawBotCluster = true; }
        }
        
        // If there exists a top and bottom cluster, this is a pair
        // event. Otherwise, it is not.
        return (sawTopCluster && sawBotCluster);
    }
    
    /**
     * Indicates whether or not a track is positively charged.
     * @param track - The track.
     * @return Returns <code>true</code> true if the track is
     * positively charged and <code>false</code> otherwise.
     */
    public static final boolean isPositive(Track track) {
        return (-Math.signum(track.getTrackStates().get(0).getOmega())) > 0;
    }
    
    /**
     * Specifies whether a cluster is located on the positron side of
     * the calorimeter.
     * @param cluster - The cluster.
     * @return Returns <code>true</code> if the cluster is on the
     * positron side of the calorimeter and <code>false</code> if it
     * is not.
     */
    public static final boolean isPositronSideCluster(Cluster cluster) {
        return TriggerModule.getClusterXIndex(cluster) > 0;
    }
    
    /**
     * Specifies whether or not a cluster is a top cluster.
     * @param cluster - The cluster.
     * @return Returns <code>true</code> for top clusters and
     * <code>false</code> for bottom clusters.
     */
    public static final boolean isTopCluster(Cluster cluster) {
        return TriggerModule.getClusterYIndex(cluster) > 0;
    }
    
    /**
     * Specifies whether or not a track is a top track. This is based
     * on the tan(Lambda) of the track.
     * @param track - The track.
     * @return Returns <code>true</code> for top tracks and
     * <code>false</code> for bottom tracks.
     */
    public static final boolean isTopTrack(Track track) {
        return track.getTrackStates().get(0).getTanLambda() > 0;
    }
    
    /**
     * Specifies whether an event meets the minimum requirements for
     * a triplet event. These requirements are that an event has a
     * pair of clusters where one cluster is on the top of the
     * calorimeter and one cluster is on the bottom. There must also
     * be one additional cluster anywhere in the calorimeter.
     * @param event - The event to test.
     * @param gtpClusterCollectionName - The name of the collection
     * containing the clusters.
     * @return Returns <code>true</code> if there exists a cluster
     * pair that meets the triplet conditions and otherwise returns
     * <code>false</code>.
     */
    public static final boolean isTripletEvent(EventHeader event, String clusterCollectionName) {
        // Get the cluster collection.
        List<Cluster> clusters = getCollection(event, clusterCollectionName, Cluster.class);
        
        // The event must contain at least three clusters.
        if(clusters.size() < 3) { return false; }
        
        // Track whether there exists any top and bottom cluster.
        boolean sawTopCluster = false;
        boolean sawBotCluster = false;
        for(Cluster cluster : clusters) {
            if(TriggerModule.getClusterYIndex(cluster) > 0) {
                sawTopCluster = true;
            } else { sawBotCluster = true; }
        }
        
        // If there exists a top and bottom cluster, this is a pair
        // event. Otherwise, it is not.
        return (sawTopCluster && sawBotCluster);
    }
    
    /**
     * Calculates the value of <code>x</code> for a polynomial with
     * coefficients defined by <code>coeff</code>. The coefficients
     * must be in order of increasing power of <code>x</code>.
     * <br/><br/>
     * For example, consider <code>coeff = { 1, 2, 3 }</code>. The
     * polynomial will take the form:
     * <br/>
     * <code>1 + 2x + 3x<sup>2</sup></code>
     * @param coeff - The coefficients, in increasing order of powers
     * of <code>x</code>.
     * @param x - The value at which the polynomial is to be
     * evaluated.
     * @return Returns the value of the polynomial at <code>x</code>.
     */
    public static final double polynomial(double[] coeff, double x) {
        double total = 0.0;
        for(int i = 0; i < coeff.length; i++) {
            total += coeff[i] * Math.pow(x, i);
        }
        return total;
    }
    
    /**
     * Checks if a given deltaR is within the allowed matching bounds
     * defined by two polynomials with coefficients as defined in
     * <code>lowCoeff</code> and <code>uppCoeff</code> (lower and
     * upper bounds respectively) for a track with momentum
     * <code>momentum</code>.
     * @param lowCoeff - The coefficients for the lower bound fit
     * function.
     * @param uppCoeff - The coefficients for the upper bound fit
     * function.
     * @param momentum - The track momentum. Defines the point in the
     * fit function at which the bounds are to be checked.
     * @param deltaR - The difference between the cluster and track
     * positions.
     * @return Returns <code>true</code> if deltaR is a match and
     * <code>false</code> otherwise.
     * @see org.hps.analysis.MC.TriggerTuningUtilityModule#polynomial(double[], double)
     */
    private static final boolean inRange(double[] lowCoeff, double[] uppCoeff, double momentum, double deltaR) {
        // Get the fit value.
        double lowFitVal = polynomial(lowCoeff, momentum);
        double uppFitVal = polynomial(uppCoeff, momentum);
        
        // The cluster matches if deltaR is between these two
        // values.
        return (deltaR >= lowFitVal && deltaR <=  uppFitVal);
    }
    
    /**
     * Checks if a cluster matches with a track with the specified
     * parameters.
     * @param cluster - The cluster.
     * @param isPositive - Whether or not the track is positive.
     * @param isTop - Whether or not the track is a top track.
     * @param momentum - The track momentum.
     * @param trackR - The track position at the calorimeter face.
     * @return Returns <code>true</code> if the track matches with
     * the cluster and <code>false</code> otherwise.
     */
    private static final boolean matches(Cluster cluster, double momentum, double[] trackR, double[] lowCoeff, double[] uppCoeff) {
        // Get the cluster position.
        double[] clusterR = cluster.getPosition();
        
        // Calculate the difference in position.
        double deltaX = clusterR[0] - trackR[0];
        double deltaY = clusterR[1] - trackR[1];
        double deltaR = Math.sqrt(Math.pow(deltaX, 2) + Math.pow(deltaY, 2));
        
        // The cluster matches if deltaR is between the
        // values defined by each fit.
        return inRange(lowCoeff, uppCoeff, momentum, deltaR);
    }
}