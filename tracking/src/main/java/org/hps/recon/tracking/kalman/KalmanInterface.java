package org.hps.recon.tracking.kalman;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import hep.physics.matrix.MatrixOp;
import hep.physics.matrix.SymmetricMatrix;
import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.Hep3Matrix;
import hep.physics.vec.VecOp;

import org.hps.recon.tracking.MaterialSupervisor.SiStripPlane;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.TrackerHit;
import org.lcsim.event.LCIOParameters.ParameterName;
import org.lcsim.event.base.BaseTrack;

public class KalmanInterface {

    public Map<Measurement, TrackerHit> hitMap;
    public Map<SiModule, HpsSiSensor> moduleMap;
    public static SquareMatrix HpsToKalman;
    public static BasicHep3Matrix HpsToKalmanMatrix;

    public KalmanInterface() {
        hitMap = new HashMap<Measurement, TrackerHit>();
        moduleMap = new HashMap<SiModule, HpsSiSensor>();
        double[][] HpsToKalmanVals = { { 0, 1, 0 }, { 1, 0, 0 }, { 0, 0, -1 } };
        HpsToKalman = new SquareMatrix(3, HpsToKalmanVals);
        HpsToKalmanMatrix = new BasicHep3Matrix();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++)
                HpsToKalmanMatrix.setElement(i, j, HpsToKalmanVals[i][j]);
        }
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

    public SiModule createSiModule(SiStripPlane inputPlane, FieldMap fm) {
        // SiModule(int Layer, Plane p, double stereo, double width, double height, double thickness, FieldMap Bfield) {

        Vec pointOnPlane = new Vec(3, inputPlane.origin().v());
        Vec normalToPlane = new Vec(3, VecOp.cross(inputPlane.getMeasuredCoordinate(), inputPlane.getUnmeasuredCoordinate()).v());

        double stereo = 0;
        HpsSiSensor temp = (HpsSiSensor) (inputPlane.getSensor());
        if (!temp.isAxial()) {
            stereo = Math.acos(normalToPlane.v[0] / normalToPlane.mag());
            if (stereo > (Math.PI / 2))
                stereo = Math.PI - stereo;
        }

        Vec normalToPlaneTransformed = normalToPlane.leftMultiply(HpsToKalman);
        Vec pointOnPlaneTransformed = pointOnPlane.leftMultiply(HpsToKalman);
        if (pointOnPlaneTransformed.v[1] < 0)
            pointOnPlaneTransformed = pointOnPlaneTransformed.scale(-1.0);

        //        System.out.printf("normalToPlaneT: %f %f %f \n", normalToPlaneTransformed.v[0], normalToPlaneTransformed.v[1], normalToPlaneTransformed.v[2]);
        SiModule newMod = new SiModule(temp.getLayerNumber(), new Plane(pointOnPlaneTransformed, normalToPlaneTransformed), stereo, inputPlane.getWidth(), inputPlane.getLength(), inputPlane.getThickness(), fm);
        moduleMap.put(newMod, temp);
        return newMod;
    }

    public void fillMeasurements(List<SiModule> mods, List<TrackerHit> hits1D) {
        Map<HpsSiSensor, ArrayList<TrackerHit>> stripHits = new HashMap<HpsSiSensor, ArrayList<TrackerHit>>();

        for (TrackerHit stripHit : hits1D) {
            HpsSiSensor temp = (HpsSiSensor) ((RawTrackerHit) stripHit.getRawHits().get(0)).getDetectorElement();
            int lay = temp.getLayerNumber();
            System.out.printf("making hit with lay %d \n", lay);
            ArrayList<TrackerHit> hitsInLayer = null;
            if (stripHits.containsKey(lay)) {
                hitsInLayer = stripHits.get(lay);
            } else {
                hitsInLayer = new ArrayList<TrackerHit>();
            }
            hitsInLayer.add(stripHit);
            stripHits.put(temp, hitsInLayer);
        }

        for (SiModule mod : mods) {
            if (!stripHits.containsKey(moduleMap.get(mod)))
                continue;
            ArrayList<TrackerHit> temp = stripHits.get(moduleMap.get(mod));
            for (TrackerHit hit : temp) {
                // hit position
                Vec hitPos = (new Vec(3, hit.getPosition())).leftMultiply(HpsToKalman);
                //System.out.printf("hit getPosition: %f %f %f \n", hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2]);
                Vec hitPosTransformed = mod.toLocal(hitPos);

                // uncertainty on position
                RotMatrix rm = mod.Rinv;

                //rm.print("modRinv");

                double[][] r = rm.M;
                Hep3Matrix rotMat = new BasicHep3Matrix(r[0][0], r[0][1], r[0][2], r[1][0], r[1][1], r[1][2], r[2][0], r[2][1], r[2][2]);
                rotMat = VecOp.mult(rotMat, HpsToKalmanMatrix);
                SymmetricMatrix cov = new SymmetricMatrix(3, hit.getCovMatrix(), true);
                double sigma2 = MatrixOp.mult(MatrixOp.mult(rotMat, cov), MatrixOp.inverse(rotMat)).e(1, 1);
                Measurement m = new Measurement(hitPosTransformed.v[1], Math.sqrt(sigma2), new Vec(0, 0, 0), 0);

                //System.out.printf("rotMat %s \n", rotMat.toString());
                //System.out.printf("cov %s \n", cov.toString());

                mod.addMeasurement(m);
                hitMap.put(m, hit);
            }
        }

    }

    //    public static FieldMap createFieldMap(org.lcsim.geometry.FieldMap inputMap) {
    //
    //    }
}
