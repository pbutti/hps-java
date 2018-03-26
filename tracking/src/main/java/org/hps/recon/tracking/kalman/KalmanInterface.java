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
import org.lcsim.event.LCIOParameters.ParameterName;
import org.lcsim.event.base.BaseTrack;

public class KalmanInterface {

    public Map<Measurement, TrackerHit> hitMap;

    public KalmanInterface() {

    }

    public BaseTrack createTrack(SeedTrack trk, List<Measurement> measList) {

        BaseTrack newTrack = new BaseTrack();
        double[] oldParams = trk.helixParams().v;
        double[] params = new double[5];
        SquareMatrix oldCov = trk.covariance();
        SymmetricMatrix cov = new SymmetricMatrix(5);

        // convert params
        params[ParameterName.d0.ordinal()] = oldParams[0];
        params[ParameterName.phi0.ordinal()] = oldParams[1];
        params[ParameterName.omega.ordinal()] = oldParams[2];
        params[ParameterName.tanLambda.ordinal()] = oldParams[4];
        params[ParameterName.z0.ordinal()] = oldParams[3];

        // convert cov matrix
        for (int i = 0; i <= 2; i++) {
            for (int j = 0; j <= 2; i++) {
                cov.setElement(i, j, oldCov.M[i][j]);
            }
        }
        for (int i = 3; i <= 4; i++) {
            for (int j = 0; j <= 4; i++) {
                cov.setElement(i, j, oldCov.M[j][i]);
                cov.setElement(j, i, oldCov.M[i][j]);
            }
        }

        newTrack.setTrackParameters(params, trk.B());
        newTrack.setCovarianceMatrix(cov);
        newTrack.setFitSuccess(trk.success);
        for (Measurement meas : measList) {
            TrackerHit hit = hitMap.get(meas);
            if (hit != null)
                newTrack.addHit(hit);
        }
        newTrack.setNDF(measList.size());

        return newTrack;
    }

    public static SiModule createSiModule(SiStripPlane inputPlane, FieldMap fm) {
        // SiModule(int Layer, Plane p, double stereo, double width, double height, double thickness, FieldMap Bfield) {

        Vec pointOnPlane = new Vec(3, inputPlane.origin().v());
        Vec normalToPlane = new Vec(3, VecOp.cross(inputPlane.getMeasuredCoordinate(), inputPlane.getUnmeasuredCoordinate()).v());

        double stereo = 0;
        HpsSiSensor temp = (HpsSiSensor) (inputPlane.getSensor());
        if (!temp.isAxial()) {
            stereo = Math.acos(normalToPlane.v[2]);
        }

        return new SiModule(temp.getLayerNumber(), new Plane(pointOnPlane, normalToPlane), stereo, inputPlane.getWidth(), inputPlane.getLength(), inputPlane.getThickness(), fm);
    }

    public void fillMeasurements(List<SiModule> mods, List<TrackerHit> hits1D) {
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
                hitMap.put(m, hit);
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
