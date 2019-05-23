package org.lcsim.detector.tracker.silicon;

import hep.physics.matrix.BasicMatrix;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;

import org.lcsim.detector.IDetectorElement;
import org.lcsim.detector.IRotation3D;
import org.lcsim.detector.ITranslation3D;
import org.lcsim.detector.RotationPassiveXYZ;
import org.lcsim.detector.Transform3D;
import org.lcsim.detector.Translation3D;
import org.lcsim.detector.identifier.IIdentifier;
import org.lcsim.detector.solids.Box;
import org.lcsim.detector.solids.Polygon3D;

/**
 * Description of the layer 0 sensors used by the SVT. This class extends
 * {@link HpsSiSensor} but overrides several properties such as strip pitch and
 * sense transfer efficiency.  It should be noted that this sensor has no
 * intermediate strips and that is why the sense transfer efficiency is 0.
 *
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */
public class HpsThinSiSensor extends HpsSiSensor {

    private final double readoutStripCapacitanceIntercept = 0;
    private final double readoutStripCapacitanceSlope = 0.16; // pf/mm
    private final double senseStripCapacitanceIntercept = 0;
    private final double senseStripCapacitanceSlope = 0.16; // pf/mm

    public HpsThinSiSensor(int sensorid, String name, IDetectorElement parent, String support, IIdentifier id) {
        super(sensorid, name, parent, support, id);
    }

    /** @return The total number of sense strips per sensor. */
    @Override
    public int getNumberOfSenseStrips() {
        return 256;
    }

    /** @return The readout strip pitch in mm. */
    @Override
    public double getReadoutStripPitch() {
        return 0.055; // mm
    }

    /** @return The sense strip pitch in mm. */
    @Override
    public double getSenseStripPitch() {
        return 0.055; // mm
    }

    /**
     * Get the charge transfer efficiency of the sense strips. The thin sensors
     * don't have an intermediate strip so the charge transfer efficiency is
     * zero.
     *
     * @return The charge transfer efficiency of the sense strips.
     */
    @Override
    public double getSenseTransferEfficiency() {
        return 0.0;
    }
    
    /**
     * Setup the geometry and electrical characteristics of an {@link HpsSiSensor}
     */
    @Override
    public void initialize() {

        // Get the solid corresponding to the sensor volume
        final Box sensorSolid = (Box) this.getGeometry().getLogicalVolume().getSolid();
        
        System.out.println("Sensor Solid " + sensorSolid.getFacesNormalTo(new BasicHep3Vector(0, 0, 1)));

        // Get the faces of the solid corresponding to the n and p sides of the sensor
        final Polygon3D pSide = sensorSolid.getFacesNormalTo(new BasicHep3Vector(0, 0, 1)).get(0);
        final Polygon3D nSide = sensorSolid.getFacesNormalTo(new BasicHep3Vector(0, 0, -1)).get(0);
        
        System.out.println("pSide " + pSide);
        System.out.println("nSide" + nSide);

        // p side collects holes.
        this.setBiasSurface(ChargeCarrier.HOLE, pSide);

        // n side collects electrons.
        this.setBiasSurface(ChargeCarrier.ELECTRON, nSide);

        Hep3Vector electrodeShift = new BasicHep3Vector(0,-sensorSolid.getYHalfLength()/2,0);
        // Translate to the outside of the sensor solid in order to setup the electrodes
        final ITranslation3D electrodesPosition1 = new Translation3D(VecOp.add(VecOp.mult(-pSide.getDistance(), pSide.getNormal()),electrodeShift));
        final ITranslation3D electrodesPosition2 = new Translation3D(VecOp.add(VecOp.mult(-pSide.getDistance(), pSide.getNormal()),VecOp.neg(electrodeShift)));
        
        System.out.println("electrodesPosition1 " + electrodesPosition1 + " " + -pSide.getDistance() + " " + pSide.getNormal());
        System.out.println("electrodesPosition2 " + electrodesPosition2);

        // Align the strips with the edge of the sensor.
        final IRotation3D electrodesRotation = new RotationPassiveXYZ(0, 0, 0);
        final Transform3D electrodesTransform1 = new Transform3D(electrodesPosition1, electrodesRotation);
        final Transform3D electrodesTransform2 = new Transform3D(electrodesPosition2, electrodesRotation);
        
        System.out.println("electrodesRotation " + electrodesRotation);
        System.out.println("electrodesTransform1 " + electrodesTransform1);
        System.out.println("electrodesTransform2 " + electrodesTransform2);

        // Set the number of readout and sense electrodes.
        final SiStrips readoutElectrodes1 = new SiStrips(ChargeCarrier.HOLE, getReadoutStripPitch(), this,electrodesTransform1);
        final SiStrips senseElectrodes1 = new SiStrips(ChargeCarrier.HOLE, getSenseStripPitch(),this.getNumberOfSenseStrips(), this, electrodesTransform1);
        final SiStrips readoutElectrodes2 = new SiStrips(ChargeCarrier.HOLE, getReadoutStripPitch(), this,electrodesTransform2);
        final SiStrips senseElectrodes2 = new SiStrips(ChargeCarrier.HOLE, getSenseStripPitch(),this.getNumberOfSenseStrips(), this, electrodesTransform2);

        System.out.println("readoutElectrodes1 " + readoutElectrodes1);
        System.out.println("senseElectrodes1 " + senseElectrodes1);
        System.out.println("readoutElectrodes2 " + readoutElectrodes2);
        System.out.println("senseElectrodes2 " + senseElectrodes2);

        final double readoutCapacitance = this.getStripLength() > this.longSensorLengthThreshold ? this.readoutLongStripCapacitanceSlope
                : this.readoutStripCapacitanceSlope;
        final double senseCapacitance = this.getStripLength() > this.longSensorLengthThreshold ? this.senseLongStripCapacitanceSlope
                : this.senseStripCapacitanceSlope;

        // Set the strip capacitance.
        readoutElectrodes1.setCapacitanceIntercept(this.readoutStripCapacitanceIntercept);
        readoutElectrodes1.setCapacitanceSlope(readoutCapacitance);
        readoutElectrodes2.setCapacitanceIntercept(this.readoutStripCapacitanceIntercept);
        readoutElectrodes2.setCapacitanceSlope(readoutCapacitance);
        senseElectrodes1.setCapacitanceIntercept(this.senseStripCapacitanceIntercept);
        senseElectrodes1.setCapacitanceSlope(senseCapacitance);
        senseElectrodes2.setCapacitanceIntercept(this.senseStripCapacitanceIntercept);
        senseElectrodes2.setCapacitanceSlope(senseCapacitance);

        // Set sense and readout electrodes.
        this.setSenseElectrodes(senseElectrodes1);
        this.setReadoutElectrodes(readoutElectrodes1);
        this.setSenseElectrodes(senseElectrodes2);
        this.setReadoutElectrodes(readoutElectrodes2);

        // Set the charge transfer efficiency of both the sense and readout
        // strips.
        final double[][] transferEfficiencies = {{this.getReadoutTransferEfficiency(), this.getSenseTransferEfficiency()}};
        this.setTransferEfficiencies(ChargeCarrier.HOLE, new BasicMatrix(transferEfficiencies));

    }
}
