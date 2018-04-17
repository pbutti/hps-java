package org.hps.analysis.geometry;

import hep.physics.vec.Hep3Vector;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.data.detectors.DetectorDataResources;
import org.hps.recon.tracking.CoordinateTransformations;
import org.hps.recon.tracking.MaterialSupervisor;
import org.lcsim.detector.ILogicalVolume;
import org.lcsim.detector.IPhysicalVolume;
import org.lcsim.detector.IPhysicalVolumeContainer;
import org.lcsim.detector.material.IMaterial;
import org.lcsim.detector.solids.Box;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiTrackerModule;
import org.lcsim.geometry.Detector;

/**
 *
 * @author Norman Graf
 */
public class SvtGeometryInspector {

    public static void main(String[] args) throws Exception {
        boolean debug = false;
        String detectorName = "HPS-EngRun2015-Nominal-v6-0-fieldmap";
        if (args.length > 0) {
            detectorName = args[0];
        }
        Set<String> availableDetectors = DetectorDataResources.getDetectorNames();
        boolean foundIt = false;
        for (String s : availableDetectors) {
            if (detectorName.equals(s)) {
                foundIt = true;
            }
        }
        if (!foundIt) {
            System.out.println("Detector name " + detectorName + " not found.");
            System.out.println("Please pick from one of the following supported detectors: ");
            for (String s : availableDetectors) {
                System.out.println(s);
            }
            throw new RuntimeException();
        }
        final DatabaseConditionsManager manager = new DatabaseConditionsManager();
        manager.setDetector(detectorName, 5772);
        Detector detector = manager.getCachedConditions(Detector.class, "compact.xml").getCachedData();
        System.out.println(detector.getName());

        MaterialSupervisor materialManager = new MaterialSupervisor();
//        MultipleScattering scattering = new MultipleScattering(materialManager);
        materialManager.buildModel(detector);
        List<MaterialSupervisor.ScatteringDetectorVolume> stripPlanes = materialManager.getMaterialVolumes();
        //TODO replace these lists with a helper class.
        List<String> names = new ArrayList<String>();
        List<Hep3Vector> oList = new ArrayList<Hep3Vector>();
        List<Hep3Vector> uList = new ArrayList<Hep3Vector>();
        List<Hep3Vector> vList = new ArrayList<Hep3Vector>();
        List<Hep3Vector> nList = new ArrayList<Hep3Vector>();
        List<Double> measDim = new ArrayList<Double>();
        List<Double> unmeasDim = new ArrayList<Double>();
        List<Boolean> isAxial = new ArrayList<Boolean>();

        for (MaterialSupervisor.ScatteringDetectorVolume vol : stripPlanes) {
            MaterialSupervisor.SiStripPlane plane = (MaterialSupervisor.SiStripPlane) vol;
            if (debug) {
                System.out.println(plane.getName());
            }
            names.add(plane.getName());
            Hep3Vector oprime = CoordinateTransformations.transformVectorToDetector(plane.origin());
            Hep3Vector nprime = CoordinateTransformations.transformVectorToDetector(plane.normal());
            if (debug) {
                System.out.println(" origin: " + oprime);
            }
            if (debug) {
                System.out.println(" normal: " + nprime);
            }
            if (debug) {
                System.out.println(" Plane is: " + plane.getMeasuredDimension() + " x " + plane.getUnmeasuredDimension());
            }
            HpsSiSensor sensor = (HpsSiSensor) plane.getSensor();

//            if (debug) {
//                System.out.println(SvtUtils.getInstance().isAxial(sensor) ? "axial" : "stereo");
//            }
            Hep3Vector measDir = CoordinateTransformations.transformVectorToDetector(plane.getMeasuredCoordinate());
            if (debug) {
                System.out.println("measured coordinate:    " + measDir);
            }
            Hep3Vector unmeasDir = CoordinateTransformations.transformVectorToDetector(plane.getUnmeasuredCoordinate());
            if (debug) {
                System.out.println("unmeasured coordinate:   " + unmeasDir);
            }
            if (debug) {
                System.out.println("thickness: " + plane.getThickness() + " in X0: " + plane.getThicknessInRL());
            }
            SiTrackerModule module = (SiTrackerModule) plane.getSensor().getGeometry().getDetectorElement().getParent();
            IPhysicalVolume parent = module.getGeometry().getPhysicalVolume();
            IPhysicalVolumeContainer daughters = parent.getLogicalVolume().getDaughters();
            if (debug) {
                System.out.printf(" found %d daughters to SiTrackerModule\n", daughters.size());
            }
            for (IPhysicalVolume daughter : daughters) {
                ILogicalVolume logicalVolume = daughter.getLogicalVolume();
                IMaterial material = logicalVolume.getMaterial();
                //System.out.println(material);
                String name = material.getName();
                double X0 = 10.0 * material.getRadiationLength() / material.getDensity();
                double X0cm = material.getRadiationLengthWithDensity();
                if (debug) {
                    System.out.println("material: " + name + " with X0= " + X0 + " mm " + X0cm);
                }
                Box solid = (Box) logicalVolume.getSolid();
                if (debug) {
                    System.out.printf(" x %f y %f z %f box\n", solid.getXHalfLength(), solid.getYHalfLength(), solid.getZHalfLength());
                }
                double halfThickness = solid.getZHalfLength();
            }
            oList.add(oprime);
            nList.add(nprime);
            uList.add(measDir);
            vList.add(unmeasDir);
            measDim.add(plane.getMeasuredDimension());
            unmeasDim.add(plane.getUnmeasuredDimension());
            isAxial.add(sensor.isAxial());
        }
        DecimalFormat df = new DecimalFormat("###.######");
        for (int i = 0; i < oList.size(); ++i) {
            Hep3Vector o = oList.get(i);
            Hep3Vector n = nList.get(i);
            Hep3Vector u = uList.get(i);
            Hep3Vector v = vList.get(i);
            double l = unmeasDim.get(i) / 2.;
            double h = measDim.get(i) / 2.;
//            String planeType = isAxial.get(i) ? "axial" : "stereo";
//            System.out.println("//" + planeType);
            System.out.println(names.get(i));
            System.out.println(df.format(o.x()) + " , " + df.format(o.y()) + " , " + df.format(o.z()) + " " + df.format(n.x()) + " , " + df.format(n.y()) + " , " + df.format(n.z()));

        }
    }
}
