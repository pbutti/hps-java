package org.hps.recon.tracking;

import hep.physics.matrix.SymmetricMatrix;
import hep.physics.vec.Hep3Vector;

import java.util.List;

import org.lcsim.detector.ITransform3D;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.TrackerHit;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;
import org.lcsim.recon.tracking.digitization.sisim.TrackerHitType;

public class HpsSiTrackerHitStrip1D extends SiTrackerHitStrip1D {

    private double timeError;

    public double getTimeError() {
        return timeError;
    }

    public void setTimeError(double input) {
        timeError = input;
    }

    // TODO
    public double calculateTimeError(List<FittedRawTrackerHit> rawHits) {
        return 0;
    }

    public HpsSiTrackerHitStrip1D(Hep3Vector position_vector, SymmetricMatrix covariance_matrix, double energy, double time, List<RawTrackerHit> raw_hits, TrackerHitType decoded_type) {
        super(position_vector, covariance_matrix, energy, time, raw_hits, decoded_type);
    }

    public HpsSiTrackerHitStrip1D(TrackerHit hit) {
        super(hit);
    }

    public HpsSiTrackerHitStrip1D(TrackerHit hit, TrackerHitType.CoordinateSystem coordinate_system) {
        super(hit, coordinate_system);
    }

    public HpsSiTrackerHitStrip1D getTransformedHit(TrackerHitType.CoordinateSystem coordinate_system) {
        return new HpsSiTrackerHitStrip1D(super.getTransformedHit(coordinate_system));
    }

    public HpsSiTrackerHitStrip1D getTransformedHit(ITransform3D global_to_local) {
        return new HpsSiTrackerHitStrip1D(super.getTransformedHit(global_to_local));
    }

}
