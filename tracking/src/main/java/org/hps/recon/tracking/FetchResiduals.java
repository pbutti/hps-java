package org.hps.recon.tracking;

import java.util.List;

import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.event.LCRelation;
import org.lcsim.event.Track;
import org.lcsim.util.Driver;

public class FetchResiduals extends Driver {
    private String trackCollectionName = "MatchedTracks";
    private String relationsCollectionName = "TrackResidualsRelations";

    @Override
    public void process(EventHeader event) {
        // fetching track and LCRelation collections
        if (!event.hasCollection(Track.class, trackCollectionName)) {
            System.out.println(trackCollectionName + " does not exist; skipping event");
            return;
        }
        List<Track> tracks = event.get(Track.class, trackCollectionName);

        if (!event.hasCollection(LCRelation.class, relationsCollectionName)) {
            System.out.println(relationsCollectionName + " does not exist; skipping event");
            return;
        }
        List<LCRelation> trackresRelation = event.get(LCRelation.class, relationsCollectionName);

        // printing LCRelations
        System.out.println("\n Printing tracks in LCRelations");
        for (LCRelation relation : trackresRelation) {
            if (relation == null)
                continue;
            Track trkTest = (Track) relation.getTo();
            System.out.printf("    Track %d \n", System.identityHashCode(trkTest));
        }
        System.out.println("");

        // printing tracks
        System.out.println("\n Printing tracks in track collection");
        for (Track trk : tracks) {
            System.out.printf("track %d \n", System.identityHashCode(trk));
            GenericObject trackRes = null;
            for (LCRelation relation : trackresRelation) {
                if (relation == null)
                    continue;
                Track trkTest = (Track) relation.getTo();
                if (trkTest == trk) {
                    System.out.printf("Relation match found for track %d \n", System.identityHashCode(trk));
                    trackRes = (GenericObject) relation.getFrom();
                }
            }
            if (trackRes == null) {
                System.out.printf("null TrackResidualsData for track %d \n", System.identityHashCode(trk));
            }
        }

        System.out.println("\n End of event");

    }
}
