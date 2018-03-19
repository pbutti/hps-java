package org.hps.recon.tracking.kalman;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import hep.physics.matrix.MatrixOp;
import hep.physics.matrix.SymmetricMatrix;
import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.VecOp;

import org.hps.recon.tracking.MaterialSupervisor.SiStripPlane;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.TrackerHit;

public class KalmanInterface {

    public static SiModule createSiModule(SiStripPlane inputPlane, FieldMap fm) {
        // SiModule(int Layer, Plane p, double stereo, double width, double height, double thickness, FieldMap Bfield) {

        //FIXME
        double stereo = 0;
        HpsSiSensor temp = (HpsSiSensor) (inputPlane.getSensor());
        //if (!temp.isAxial())
        // get stereo angle

        Vec pointOnPlane = new Vec(3, inputPlane.origin().v());
        Vec normalToPlane = new Vec(3, VecOp.cross(inputPlane.getMeasuredCoordinate(), inputPlane.getUnmeasuredCoordinate()).v());

        return new SiModule(temp.getLayerNumber(), new Plane(pointOnPlane, normalToPlane), stereo, inputPlane.getWidth(), inputPlane.getLength(), inputPlane.getThickness(), fm);
    }

    public static void fillMeasurements(List<SiModule> mods, List<TrackerHit> hits1D) {
        Map<Integer, ArrayList<TrackerHit>> stripHits = new HashMap<Integer, ArrayList<TrackerHit>>();

        for (TrackerHit stripHit : hits1D) {
            int lay = ((HpsSiSensor) ((RawTrackerHit) stripHit.getRawHits().get(0)).getDetectorElement()).getLayerNumber();

            ArrayList<TrackerHit> hitsInLayer = null;
            if (stripHits.containsKey(lay)) {
                hitsInLayer = stripHits.get(lay);
            } else {
                hitsInLayer = new ArrayList<TrackerHit>();
            }
            hitsInLayer.add(stripHit);
            stripHits.put(lay, hitsInLayer);
        }

        for (ArrayList<TrackerHit> hitsInLayer : stripHits.values()) {
            hitsInLayer.sort(HitComparator);
        }

        for (SiModule mod : mods) {
            if (!stripHits.containsKey(mod.Layer))
                continue;
            ArrayList<TrackerHit> temp = stripHits.get(mod.Layer);
            for (TrackerHit hit : temp) {
                // hit position
                RotMatrix rm = mod.Rinv;
                Vec hitPos = new Vec(3, hit.getPosition());
                Vec hitPosTransformed = rm.rotate(hitPos);

                // uncertainty on position
                double[][] r = rm.M;
                BasicHep3Matrix rotMat = new BasicHep3Matrix(r[0][0], r[0][1], r[0][2], r[1][0], r[1][1], r[1][2], r[2][0], r[2][1], r[2][2]);
                SymmetricMatrix cov = new SymmetricMatrix(3, hit.getCovMatrix(), true);
                double sigma2 = MatrixOp.mult(MatrixOp.mult(rotMat, cov), MatrixOp.inverse(rotMat)).e(0, 0);
                Measurement m = new Measurement(hitPosTransformed.v[0], Math.sqrt(sigma2), null, 0);

                mod.addMeasurement(m);
            }
        }

    }

    public static Comparator<TrackerHit> HitComparator = new Comparator<TrackerHit>() {

        public int compare(TrackerHit h1, TrackerHit h2) {
            double x1 = h1.getPosition()[0];
            double x2 = h2.getPosition()[0];

            // ascending order?
            if (x1 > x2)
                return 1;
            else if (x1 < x2)
                return -1;
            else
                return 0;
        }
    };

    //    public static FieldMap createFieldMap(org.lcsim.geometry.FieldMap inputMap) {
    //
    //    }
}
