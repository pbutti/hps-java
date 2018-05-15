package org.hps.recon.tracking;

// /nfs/slac/g/hps3/data/engrun2015/pass8/svtTests/

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

import org.lcsim.event.EventHeader;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.util.loop.LCIODriver;
import org.lcsim.util.loop.LCSimLoop;
import org.lcsim.util.test.TestUtil.TestOutputFile;

public class AddTrackCollections extends TestCase {
    // /nfs/slac/work/mdiamond
    protected String extraFileName = "/nfs/slac/g/hps3/data/engrun2015/pass8/svtTests/hps_005772.0.twoClusterSkim.evio_timingcut24.slcio";
    protected String mainFileName = "/nfs/slac/g/hps3/data/engrun2015/pass8/svtTests/hps_005772.0.twoClusterSkim.evio_timingcut15.slcio";
    protected String outputFileName = "tst_combined.slcio";
    protected String TrackCollName = "GBLTracks";
    protected String HitCollName = "RotatedHelicalTrackHits";
    protected String suffix = "-24s";
    protected long nEvents = 5000;

    public void testClear() throws Exception {
        String[] extraFileNameList = { extraFileName };

        File mainFile = new File(mainFileName);
        File outputFile = new TestOutputFile(outputFileName);
        //        outputFile.getParentFile().mkdirs();

        LCSimLoop loop = new LCSimLoop();
        loop.setLCIORecordSource(mainFile);

        AddTrackCollectionsDriver atc = new AddTrackCollectionsDriver();
        atc.setOverlayFiles(extraFileNameList);
        atc.setTrackCollName(TrackCollName);
        atc.setHitCollName(HitCollName);
        atc.setSuffix(suffix);
        loop.add(atc);

        loop.add(new LCIODriver(outputFile));

        try {
            loop.loop(nEvents);
        } catch (Exception e) {
            System.out.println(e.toString());
        }
        loop.dispose();
    }

    protected class AddTrackCollectionsDriver extends org.lcsim.util.OverlayDriver {

        protected String trackCollName = "GBLTracks";
        protected String hitCollName = "RotatedHelicalTrackHits";
        protected String suffix = "-20s";

        @Override
        protected void detectorChanged(Detector detector) {
            // nothing to do
        }

        public void setTrackCollName(String input) {
            trackCollName = input;
        }

        public void setHitCollName(String input) {
            hitCollName = input;
        }

        public void setSuffix(String input) {
            suffix = input;
        }

        @Override
        protected void process(EventHeader event) {
            EventHeader extraEvent = null;
            extraEvent = getNextEvent(overlayEvents);
            List<Track> extraTracks = null;
            List<TrackerHit> extraHits = null;
            // List<Track> extraTracksOutput = null;
            if (extraEvent != null) {
                //System.out.println("Got extra event");
                extraTracks = extraEvent.get(Track.class, trackCollName);
                extraHits = extraEvent.get(TrackerHit.class, hitCollName);
                // extraTracksOutput = new ArrayList<Track>();
            }
            int flags = 0;
            if (extraHits != null) {
                //System.out.println("Copying hits");
                flags = extraEvent.getMetaData(extraHits).getFlags();
                event.put(hitCollName + suffix, extraHits, TrackerHit.class, flags);
            }
            if (extraTracks != null) {
                for (Track trk : extraTracks) {
                    //System.out.println("Copying track");
                    List<TrackerHit> newHits = event.get(TrackerHit.class, hitCollName + suffix);
                    copyTrack(trk, newHits);
                    //                extraTracksOutput.add(newTrk);
                }

                flags = extraEvent.getMetaData(extraTracks).getFlags();
                event.put(trackCollName + suffix, extraTracks, Track.class, flags);
            }

        }

        //        @SuppressWarnings("deprecation")
        private void copyTrack(Track trk, List<TrackerHit> newHits) {
            //            BaseTrack newTrack = new BaseTrack();
            //            newTrack.setTrackParameters(trk.getTrackParameters(), trk.B());
            //            newTrack.setCovarianceMatrix(trk.getErrorMatrix());
            //            newTrack.setFitSuccess(trk.fitSuccess());
            //            newTrack.setChisq(trk.getChi2());
            //            newTrack.setTrackType(trk.getType());
            //            newTrack.setNDF(trk.getNDF());
            //            newTrack.setReferencePoint(trk.getReferencePoint());
            //            newTrack.setRefPointIsDCA(trk.isReferencePointPCA());

            List<TrackerHit> hitsToAdd = new ArrayList<TrackerHit>();
            for (TrackerHit oldHit : trk.getTrackerHits()) {
                double[] oldHitPos = oldHit.getPosition();
                //System.out.printf("old hit %f %f %f \n", oldHitPos[0], oldHitPos[1], oldHitPos[2]);
                for (TrackerHit newHit : newHits) {
                    double[] newHitPos = newHit.getPosition();

                    //  System.out.printf("new hit %f %f %f \n", newHitPos[0], newHitPos[1], newHitPos[2]);
                    if ((oldHitPos[0] == newHitPos[0]) && (oldHitPos[1] == newHitPos[1]) && (oldHitPos[2] == newHitPos[2])) {
                        hitsToAdd.add(newHit);
                        //    System.out.println("added hit");
                        break;
                    }
                }
            }

            trk.getTrackerHits().clear();
            trk.getTrackerHits().addAll(hitsToAdd);

        }

    }
}
