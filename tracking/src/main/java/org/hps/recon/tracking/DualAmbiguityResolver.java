package org.hps.recon.tracking;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;

public class DualAmbiguityResolver extends SimpleAmbiguityResolver {
    protected List<Track> dualDups;

    public DualAmbiguityResolver(List<List<Track>> inputTracks, AmbiMode mode, int share, double score) {
        super(inputTracks, mode, share, score);
        dualDups = new ArrayList<Track>();
    }

    public List<Track> getDualDuplicateTracks() {
        return dualDups;
    }

    public boolean isShared(Track trk) {
        if (!sharedTracksMap.containsKey(trk))
            return false;
        List<Track> shared = sharedTracksMap.get(trk);
        if (shared.isEmpty())
            return false;
        return true;
    }

    @Override
    protected boolean areShared(Track trk1, Track trk2) {
        if (numSharedHits(trk1, trk2) > getShareThreshold())
            return true;
        return false;
    }

    public static int numSharedHits(Track trk1, Track trk2) {
        int counter = 0;
        List<TrackerHit> hits1 = trk1.getTrackerHits();
        List<TrackerHit> hits2 = trk2.getTrackerHits();

        for (TrackerHit hit1 : hits1) {
            for (TrackerHit hit2 : hits2) {
                if (isEqual(hit1, hit2)) {
                    counter++;
                    break;
                }
            }
        }
        return counter;
    }

    private boolean isEqual(Set<TrackerHit> h1List, Set<TrackerHit> h2List) {
        if (h1List.size() != h2List.size())
            return false;
        return isPartial(h1List, h2List);
    }

    private boolean isPartial(Set<TrackerHit> h1List, Set<TrackerHit> h2List) {
        for (TrackerHit h1 : h1List) {
            boolean found = false;
            for (TrackerHit h2 : h2List) {
                if (isEqual(h1, h2)) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                return false;
            }
        }
        return true;
    }

    public static boolean isEqual(TrackerHit h1, TrackerHit h2) {
        double[] pos1 = h1.getPosition();
        double[] pos2 = h2.getPosition();
        double TOL = 0.001;

        if (Math.abs(pos1[0] - pos2[0]) > TOL)
            return false;
        if (Math.abs(pos1[1] - pos2[1]) > TOL)
            return false;
        if (Math.abs(pos1[2] - pos2[2]) > TOL)
            return false;

        return true;
    }

    public void addToTrackList(List<Track> tracklist) {
        utils.makeTrackScoreMap(tracklist);
        dualDups.clear();
        for (Track trk : tracklist) {
            //List<TrackerHit> mapEntry = trk.getTrackerHits();
            Set<TrackerHit> mapEntrySet = new HashSet<TrackerHit>(trk.getTrackerHits());

            // update shared map
            List<Track> newShared = new ArrayList<Track>();
            for (Track otherTrack : this.tracks) {
                if (areShared(trk, otherTrack)) {
                    newShared.add(otherTrack);
                    if (sharedTracksMap.containsKey(otherTrack)) {
                        sharedTracksMap.get(otherTrack).add(trk);
                    }
                }

                // update partials list
                // is this a partial of an existing track?
                // or is an existing track a partial of this one?
                if (otherTrack.getTrackerHits().size() > mapEntrySet.size()) {
                    if (isPartial(mapEntrySet, new HashSet<TrackerHit>(otherTrack.getTrackerHits()))) {
                        partials.add(trk);
                    }
                } else if (otherTrack.getTrackerHits().size() < mapEntrySet.size()) {
                    if (isPartial(new HashSet<TrackerHit>(otherTrack.getTrackerHits()), mapEntrySet)) {
                        partials.add(otherTrack);
                    }
                }
            }
            sharedTracksMap.put(trk, newShared);

            // update map for duplicates
            if (hitsToTracksMap.containsValue(mapEntrySet)) {
                List<Track> dups = hitsToTracksMap.get(mapEntrySet);
                dups.add(trk);
                for (Track dup : dups) {
                    if (!dualDups.contains(dup))
                        dualDups.add(dup);
                }
            } else {
                Set<TrackerHit> dualEntry = null;
                Set<Set<TrackerHit>> existingHitLists = hitsToTracksMap.keySet();
                for (Set<TrackerHit> existingHitList : existingHitLists) {
                    if (isEqual(existingHitList, mapEntrySet)) {
                        dualEntry = existingHitList;
                        break;
                    }
                }
                if (dualEntry != null) {
                    List<Track> dups = hitsToTracksMap.get(dualEntry);
                    dups.add(trk);
                    for (Track dup : dups) {
                        if (!dualDups.contains(dup))
                            dualDups.add(dup);
                    }
                } else {
                    List<Track> newList = new ArrayList<Track>();
                    newList.add(trk);
                    hitsToTracksMap.put(mapEntrySet, newList);
                }
            }
        }

        // update operable tracks
        tracks.addAll(tracklist);
    }

}
