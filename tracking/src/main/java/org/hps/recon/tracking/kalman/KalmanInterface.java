package org.hps.recon.tracking.kalman;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

//import hep.physics.matrix.MatrixOp;
import hep.physics.matrix.SymmetricMatrix;
import hep.physics.vec.BasicHep3Matrix;
//import hep.physics.vec.BasicHep3Vector;
//import hep.physics.vec.Hep3Matrix;
import hep.physics.vec.VecOp;

//import org.hps.recon.tracking.CoordinateTransformations;
import org.hps.recon.tracking.MaterialSupervisor.SiStripPlane;
import org.hps.recon.tracking.TrackUtils;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.event.LCIOParameters.ParameterName;
import org.lcsim.event.base.BaseTrack;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;
import org.lcsim.recon.tracking.digitization.sisim.TrackerHitType;

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

    public void clearMaps() {
        hitMap.clear();
        moduleMap.clear();
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
            for (int j = 0; j <= 2; j++) {
                cov.setElement(i, j, oldCov.M[i][j]);
            }
        }
        for (int i = 3; i <= 4; i++) {
            for (int j = 0; j <= 4; j++) {
                cov.setElement(i, j, oldCov.M[j][i]);
                cov.setElement(j, i, oldCov.M[i][j]);
            }
        }

        newTrack.setTrackParameters(params, trk.B());
        newTrack.setCovarianceMatrix(cov);
        newTrack.setFitSuccess(trk.success);
        List<TrackerHit> hitsOnTrack = new ArrayList<TrackerHit>();
        for (Measurement meas : measList) {
            TrackerHit hit = hitMap.get(meas);
            if (hit != null) {
                if (!hitsOnTrack.contains(hit)) {
                    newTrack.addHit(hit);
                    hitsOnTrack.add(hit);
                }
            }
        }
        newTrack.setNDF(measList.size());

        return newTrack;
    }

    public SiModule createSiModule(SiStripPlane inputPlane, FieldMap fm) {
        // SiModule(int Layer, Plane p, double stereo, double width, double height, double thickness, FieldMap Bfield) {

        Vec pointOnPlane = new Vec(3, inputPlane.origin().v());
        Vec normalToPlane = new Vec(3, VecOp.cross(inputPlane.getMeasuredCoordinate(), inputPlane.getUnmeasuredCoordinate()).v());

        Vec normalToPlaneTransformed = normalToPlane.leftMultiply(HpsToKalman);
        Vec pointOnPlaneTransformed = pointOnPlane.leftMultiply(HpsToKalman);

        double stereo = 0;
        HpsSiSensor temp = (HpsSiSensor) (inputPlane.getSensor());
        if (temp.isAxial()) {
            stereo = Math.acos(normalToPlane.v[0] / normalToPlane.mag());
            if (stereo > (Math.PI / 2))
                stereo = Math.PI - stereo;
            //pointOnPlaneTransformed = pointOnPlaneTransformed.scale(-1.0);
            normalToPlaneTransformed = normalToPlaneTransformed.scale(-1.0);
        }

        Plane p = new Plane(pointOnPlaneTransformed, normalToPlaneTransformed);
        //        System.out.printf("normalToPlaneT: %f %f %f \n", normalToPlaneTransformed.v[0], normalToPlaneTransformed.v[1], normalToPlaneTransformed.v[2]);
        SiModule newMod = new SiModule(temp.getLayerNumber(), p, stereo, inputPlane.getWidth(), inputPlane.getLength(), inputPlane.getThickness(), fm);
        moduleMap.put(newMod, temp);
        return newMod;
    }

    public void fillMeasurements(List<SiModule> mods, Track track, RelationalTable hitToStrips, RelationalTable hitToRotated) {
        Map<HpsSiSensor, ArrayList<TrackerHit>> hitsMap = new HashMap<HpsSiSensor, ArrayList<TrackerHit>>();

        List<TrackerHit> hits1D = TrackUtils.getStripHits(track, hitToStrips, hitToRotated);

        //eliminate me
        List<TrackerHit> hits2D = track.getTrackerHits();
        for (TrackerHit hit2D : hits2D) {
            HpsSiSensor temp = ((HpsSiSensor) ((RawTrackerHit) hit2D.getRawHits().get(0)).getDetectorElement());
            int lay = temp.getLayerNumber();
            System.out.printf("2Dhit Layer %d   Position: %f %f %f \n", lay, hit2D.getPosition()[0], hit2D.getPosition()[1], hit2D.getPosition()[2]);

        }

        for (TrackerHit hit1D : hits1D) {
            HpsSiSensor temp = ((HpsSiSensor) ((RawTrackerHit) hit1D.getRawHits().get(0)).getDetectorElement());
            int lay = temp.getLayerNumber();
            //       System.out.printf("hit Layer %d   Position: %f %f %f \n", lay, hit1D.getPosition()[0], hit1D.getPosition()[1], hit1D.getPosition()[2]);

            ArrayList<TrackerHit> hitsInLayer = null;
            if (hitsMap.containsKey(lay)) {
                hitsInLayer = hitsMap.get(lay);
            } else {
                hitsInLayer = new ArrayList<TrackerHit>();
            }
            hitsInLayer.add(hit1D);
            hitsMap.put(temp, hitsInLayer);
        }

        for (SiModule mod : mods) {
            if (!hitsMap.containsKey(moduleMap.get(mod)))
                continue;
            ArrayList<TrackerHit> temp = hitsMap.get(moduleMap.get(mod));
            if (temp == null)
                continue;
            for (int i = 0; i < temp.size(); i++) {
                TrackerHit hit = temp.get(i);

                SiTrackerHitStrip1D local = (new SiTrackerHitStrip1D(hit)).getTransformedHit(TrackerHitType.CoordinateSystem.SENSOR);
                //SiTrackerHitStrip1D global = (new SiTrackerHitStrip1D(hit)).getTransformedHit(TrackerHitType.CoordinateSystem.GLOBAL);
                double umeas = local.getPosition()[0];
                double du = Math.sqrt(local.getCovarianceAsMatrix().diagonal(0));
                if (mod.Layer % 2 == 0)
                    umeas *= -1.0;
                //                System.out.printf("hit local umeas %f du %f \n", umeas, du);
                //                System.out.printf("hit global position: %f %f %f \n", global.getPosition()[0], global.getPosition()[1], global.getPosition()[2]);
                // hit position
                //                System.out.printf("hit raw position: %f %f %f \n", hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2]);
                //                Vec hitPos0 = new Vec(3, CoordinateTransformations.transformVectorToTracking(new BasicHep3Vector(hit.getPosition())).v());
                //                System.out.printf("hit transformed position: %f %f %f \n", hitPos0.v[0], hitPos0.v[1], hitPos0.v[2]);
                //                Vec hitPos = hitPos0.leftMultiply(HpsToKalman);
                //                System.out.printf("hit kalman position: %f %f %f \n", hitPos.v[0], hitPos.v[1], hitPos.v[2]);
                //                Vec hitPos2 = hitPos.dif(mod.p.X());
                //                System.out.printf("hit subtracted position: %f %f %f \n", hitPos2.v[0], hitPos2.v[1], hitPos2.v[2]);
                //                Vec hitPosTransformed = mod.Rinv.rotate(hitPos2);
                //Vec hitPosTransformed = mod.toLocal(hitPos);

                // uncertainty on position
                //RotMatrix rm = mod.Rinv;
                //rm.print("modRinv");

                //                double[][] r = rm.M;
                //                Hep3Matrix rotMat = new BasicHep3Matrix(r[0][0], r[0][1], r[0][2], r[1][0], r[1][1], r[1][2], r[2][0], r[2][1], r[2][2]);
                //                rotMat = VecOp.mult(rotMat, HpsToKalmanMatrix);
                //                SymmetricMatrix cov = new SymmetricMatrix(3, hit.getCovMatrix(), true);
                //                double sigma2 = MatrixOp.mult(MatrixOp.mult(rotMat, cov), MatrixOp.inverse(rotMat)).e(1, 1);
                Measurement m = new Measurement(umeas, du, new Vec(0, 0, 0), 0);

                //System.out.printf("rotMat %s \n", rotMat.toString());
                //System.out.printf("cov %s \n", cov.toString());

                mod.addMeasurement(m);
                hitMap.put(m, hit);

            }
        }

    }

}
