package org.lcsim.detector.tracker.silicon;

import java.util.ArrayList;
import java.util.List;

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
import org.lcsim.detector.solids.Point3D;
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

        System.out.println("Half Length Y " + sensorSolid.getYHalfLength());
        Hep3Vector electrodeShift = new BasicHep3Vector(0,-sensorSolid.getYHalfLength()/2,0);
        //Hep3Vector electrodeShift = new BasicHep3Vector(0,0,0);
        // Translate to the outside of the sensor solid in order to setup the electrodes
        final ITranslation3D electrodesPosition1 = new Translation3D(VecOp.add(VecOp.mult(-pSide.getDistance(), pSide.getNormal()),electrodeShift));
        final ITranslation3D electrodesPosition2 = new Translation3D(VecOp.add(VecOp.mult(-pSide.getDistance(), pSide.getNormal()),VecOp.neg(electrodeShift)));
        
        System.out.println("electrodesPosition1 " + electrodesPosition1 + " " + -pSide.getDistance() + " " + pSide.getNormal());
        System.out.println("electrodesPosition2 " + electrodesPosition2);

        // Align the strips with the edge of the sensor.
        final IRotation3D electrodesRotation1 = new RotationPassiveXYZ(0, 0, 0);
        final IRotation3D electrodesRotation2 = new RotationPassiveXYZ(0, 0, 0);
        final Transform3D electrodesTransform1 = new Transform3D(electrodesPosition1, electrodesRotation1);
        final Transform3D electrodesTransform2 = new Transform3D(electrodesPosition2, electrodesRotation2);
        
        System.out.println("electrodesRotation1 " + electrodesRotation1);
        System.out.println("electrodesRotation2 " + electrodesRotation2);
        System.out.println("electrodesTransform1 " + electrodesTransform1);
        System.out.println("electrodesTransform2 " + electrodesTransform2);

        // Set the number of readout and sense electrodes.
        final SiStrips readoutElectrodes1 = new SiStrips(ChargeCarrier.HOLE, getReadoutStripPitch(), this,electrodesTransform1);
        final SiStrips senseElectrodes1 = new SiStrips(ChargeCarrier.HOLE, getSenseStripPitch(),this.getNumberOfSenseStrips(), this, electrodesTransform1);
        final SiStrips readoutElectrodes2 = new SiStrips(ChargeCarrier.HOLE, getReadoutStripPitch(), this,electrodesTransform2);
        final SiStrips senseElectrodes2 = new SiStrips(ChargeCarrier.HOLE, getSenseStripPitch(),this.getNumberOfSenseStrips(), this, electrodesTransform2);

        System.out.println("Geometry1 " + readoutElectrodes1.getGeometry());
        
        List<Point3D> vertices1 = readoutElectrodes1.getGeometry().getVertices();
        Point3D point1_1 = new Point3D(new BasicHep3Vector(vertices1.get(0).x(),vertices1.get(0).y()-sensorSolid.getYHalfLength()/2,vertices1.get(0).z()));
        Point3D point2_1 = new Point3D(new BasicHep3Vector(vertices1.get(1).x(),vertices1.get(1).y()-sensorSolid.getYHalfLength()/2,vertices1.get(1).z()));
        Point3D point3_1 = new Point3D(new BasicHep3Vector(vertices1.get(2).x(),vertices1.get(2).y()+sensorSolid.getYHalfLength()/2,vertices1.get(2).z()));
        Point3D point4_1 = new Point3D(new BasicHep3Vector(vertices1.get(3).x(),vertices1.get(3).y()+sensorSolid.getYHalfLength()/2,vertices1.get(3).z()));
        //System.out.println("Point1 " + point1);
        //System.out.println("Point2 " + point2);
        //System.out.println("Point3 " + point3);
        //System.out.println("Point4 " + point4);
        List<Point3D> newvertices1 = new ArrayList<>();
        newvertices1.add(point1_1);
        newvertices1.add(point2_1);
        newvertices1.add(point3_1);
        newvertices1.add(point4_1);
        
        System.out.println("Geometry2 " + readoutElectrodes2.getGeometry());
        
        List<Point3D> vertices2 = readoutElectrodes2.getGeometry().getVertices();
        Point3D point1_2 = new Point3D(new BasicHep3Vector(vertices2.get(0).x(),vertices2.get(0).y()-sensorSolid.getYHalfLength()/2,vertices2.get(0).z()));
        Point3D point2_2 = new Point3D(new BasicHep3Vector(vertices2.get(1).x(),vertices2.get(1).y()-sensorSolid.getYHalfLength()/2,vertices2.get(1).z()));
        Point3D point3_2 = new Point3D(new BasicHep3Vector(vertices2.get(2).x(),vertices2.get(2).y()+sensorSolid.getYHalfLength()/2,vertices2.get(2).z()));
        Point3D point4_2 = new Point3D(new BasicHep3Vector(vertices2.get(3).x(),vertices2.get(3).y()+sensorSolid.getYHalfLength()/2,vertices2.get(3).z()));
        //System.out.println("Point1 " + point1);
        //System.out.println("Point2 " + point2);
        //System.out.println("Point3 " + point3);
        //System.out.println("Point4 " + point4);
        List<Point3D> newvertices2 = new ArrayList<>();
        newvertices2.add(point1_2);
        newvertices2.add(point2_2);
        newvertices2.add(point3_2);
        newvertices2.add(point4_2);
        
        readoutElectrodes1.setGeometry(new Polygon3D(newvertices1));
        senseElectrodes1.setGeometry(new Polygon3D(newvertices1));
        readoutElectrodes2.setGeometry(new Polygon3D(newvertices2));
        senseElectrodes2.setGeometry(new Polygon3D(newvertices2));
        
        System.out.println("Geometry Redo Readout 1" + readoutElectrodes1.getGeometry());
        System.out.println("Geometry Redo Sense 1" + senseElectrodes1.getGeometry());
        System.out.println("Geometry Redo Readout 2" + readoutElectrodes2.getGeometry());
        System.out.println("Geometry Redo Sense 2" + senseElectrodes2.getGeometry());

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

        /*List<SiStrips> senseElectrodes = new ArrayList<>();
        List<SiStrips> readoutElectrodes = new ArrayList<>();
        senseElectrodes.add(senseElectrodes1);
        senseElectrodes.add(senseElectrodes2);
        readoutElectrodes.add(readoutElectrodes1);
        readoutElectrodes.add(readoutElectrodes2);*/
        
        
        // Set sense and readout electrodes.
        this.setSenseElectrodes(senseElectrodes1.getChargeCarrier(),senseElectrodes1);
        this.setReadoutElectrodes(readoutElectrodes1.getChargeCarrier(),readoutElectrodes1);
        this.setSenseElectrodes(senseElectrodes2.getChargeCarrier(),senseElectrodes2);
        this.setReadoutElectrodes(readoutElectrodes2.getChargeCarrier(),readoutElectrodes2);
        
        System.out.println("Sense Electrodes 1 " + senseElectrodes1 + "  readoutElectrodes 1 " + readoutElectrodes1);
        System.out.println("Sense Electrodes 2 " + senseElectrodes2 + "  readoutElectrodes 2 " + readoutElectrodes2);
        System.out.println("Sense Electrodes " + this.getSenseElectrodes() + "  readoutElectrodes " + this.getReadoutElectrodes());

        // Set the charge transfer efficiency of both the sense and readout
        // strips.
        final double[][] transferEfficiencies = {{this.getReadoutTransferEfficiency(), this.getSenseTransferEfficiency()}};
        this.setTransferEfficiencies(ChargeCarrier.HOLE, new BasicMatrix(transferEfficiencies));

    }
    
    public void setSenseElectrodes(ChargeCarrier carrier, SiSensorElectrodes sense_electrodes)
    {
        _sense_electrodes.put(carrier,sense_electrodes);
    }
    
    public void setReadoutElectrodes(ChargeCarrier carrier, SiSensorElectrodes readout_electrodes)
    {
        _readout_electrodes.put(carrier,readout_electrodes);
    }
}
