package org.hps.analysis.geometry;

import hep.physics.vec.Hep3Vector;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.hps.detector.ecal.EcalCrystal;
import org.hps.detector.ecal.HPSEcalDetectorElement;
import org.hps.recon.tracking.CoordinateTransformations;
import org.hps.recon.tracking.MaterialSupervisor;
import org.hps.recon.tracking.MultipleScattering;
import org.hps.recon.tracking.gbl.GBLStripClusterData;
import org.hps.recon.tracking.gbl.GBLTrackData;
import org.lcsim.detector.DetectorElementStore;
import org.lcsim.detector.IDetectorElement;
import org.lcsim.detector.IGeometryInfo;
import org.lcsim.detector.ILogicalVolume;
import org.lcsim.detector.IPhysicalVolume;
import org.lcsim.detector.IPhysicalVolumeContainer;
import org.lcsim.detector.identifier.IExpandedIdentifier;
import org.lcsim.detector.identifier.IIdentifier;
import org.lcsim.detector.identifier.IIdentifierDictionary;
import org.lcsim.detector.identifier.IIdentifierHelper;
import org.lcsim.detector.material.IMaterial;
import org.lcsim.detector.solids.Box;
import org.lcsim.detector.solids.LineSegment3D;
import org.lcsim.detector.solids.Point3D;
import org.lcsim.detector.solids.Trd;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.detector.tracker.silicon.SiTrackerModule;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.event.LCRelation;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.subdetector.HPSEcal3;
import org.lcsim.recon.tracking.digitization.sisim.SiTrackerHitStrip1D;
import org.lcsim.util.Driver;

/**
 *
 * @author Norman A. Graf
 */
public class DrawDetector extends Driver {

    MaterialSupervisor _materialManager = null;
    MultipleScattering _scattering = null;

    private HPSEcalDetectorElement _api = null;
    private HPSEcal3 _ecal = null;

    private final String trackCollectionName = "MatchedTracks";
    private final String track2GblTrackRelationName = "TrackToGBLTrack";
    private final String gblTrack2StripRelationName = "GBLTrackToStripData";

    private boolean _debug = false;

    private boolean _drawDetector = true;
    private boolean _drawEvent = true;

    protected void detectorChanged(Detector det) {
        _materialManager = new MaterialSupervisor();
        _scattering = new MultipleScattering(_materialManager);
        _materialManager.buildModel(det);
        _ecal = (HPSEcal3) det.getSubdetector("Ecal");
        _api = (HPSEcalDetectorElement) _ecal.getDetectorElement();
        drawDetector(det);
    }

    private void drawDetector(Detector det) {
        boolean debug = false;
        List<MaterialSupervisor.ScatteringDetectorVolume> stripPlanes = _materialManager.getMaterialVolumes();
        //TODO replace these lists with a helper class.
        List<Hep3Vector> oList = new ArrayList<Hep3Vector>();
        List<Hep3Vector> uList = new ArrayList<Hep3Vector>();
        List<Hep3Vector> vList = new ArrayList<Hep3Vector>();
        List<Double> measDim = new ArrayList<Double>();
        List<Double> unmeasDim = new ArrayList<Double>();
        List<Boolean> isAxial = new ArrayList<Boolean>();

        for (MaterialSupervisor.ScatteringDetectorVolume vol : stripPlanes) {
            MaterialSupervisor.SiStripPlane plane = (MaterialSupervisor.SiStripPlane) vol;
            if (debug) {
                System.out.println(plane.getName());
            }
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
                System.out.printf("%s found %d daughters to SiTrackerModule\n", this.getClass().getSimpleName(), daughters.size());
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
                    System.out.printf("%s x %f y %f z %f box\n", this.getClass().getSimpleName(), solid.getXHalfLength(), solid.getYHalfLength(), solid.getZHalfLength());
                }
                double halfThickness = solid.getZHalfLength();
            }
            oList.add(oprime);
            uList.add(measDir);
            vList.add(unmeasDir);
            measDim.add(plane.getMeasuredDimension());
            unmeasDim.add(plane.getUnmeasuredDimension());
            isAxial.add(sensor.isAxial());
        }

        // now get ECal
        List<List<Hep3Vector>> calcells = new ArrayList<List<Hep3Vector>>();
        IIdentifierHelper helper = _api.getIdentifierHelper();
        for (EcalCrystal crystal : _api.getCrystals()) {
            IGeometryInfo geom = crystal.getGeometry();
            Trd trd = (Trd) geom.getLogicalVolume().getSolid();
            List<Point3D> vertices = trd.getVertices();
            List<Hep3Vector> verts = new ArrayList<Hep3Vector>();
            for (Point3D p : vertices) {
                verts.add(geom.transformLocalToGlobal(p.getHep3Vector()));
            }
            calcells.add(verts);
            //
//            IExpandedIdentifier expandedId = helper.createExpandedIdentifier();
//            expandedId.setValue(helper.getFieldIndex("ix"), crystal.getX());
//            expandedId.setValue(helper.getFieldIndex("iy"), crystal.getY());
//            expandedId.setValue(helper.getFieldIndex("system"), _api.getSystemID());
//            IIdentifier id = helper.pack(expandedId);
//            if (id.getValue() != crystal.getIdentifier().getValue()) {
//                throw new RuntimeException("Reencoded ID " + id.getValue() + " does not match crystal ID " + crystal.getIdentifier().getValue());
//            }
//
//            IDetectorElementStore deStore = DetectorElementStore.getInstance();
//            if (deStore.find(crystal.getIdentifier()).size() == 0) {
//                throw new RuntimeException("Failed to find crystal ID in store.");
//            }
//
//            if (deStore.find(id).size() == 0) {
//                throw new RuntimeException("Failed to find repacked ID in store.");
//            }
        }
//        System.out.println("calcells has " + calcells.size() + " entries");
//        // test print...
//        List<Hep3Vector> vs0 = calcells.get(0);
//        System.out.println("cell 0 has " + vs0.size() + "points");
//        System.out.println(vs0);
//        List<Hep3Vector> vs1 = calcells.get(1);
//        System.out.println("cell 1 has " + vs1.size() + "points");
//        System.out.println(vs1);

        // now draw
        if (_drawDetector) {
            //TODO extract the following into a separate method
            System.out.println("\n\n\n");
            System.out.println("size(600);");
            System.out.println("import three;");
            System.out.println("currentprojection=orthographic(");
            System.out.println("camera=(-1250,850,-1500),");
            System.out.println("up=(0.374061249845773,1.48257908496856,0.568508408053946),");
            System.out.println("target=(0,0,0),");
            System.out.println("zoom=1);");

            System.out.println("");
            System.out.println("pen fill=gray(0.6)+opacity(0.2);");
            System.out.println("pen edge=black;");
            System.out.println("pen axialEdge=red;");
            System.out.println("pen stereoEdge=blue;");
            System.out.println("// o is the origin point in the plane, u is the measurement direction, v is the unmeasured direction");
            System.out.println("// clockwise normal orientation corresponds to z = u x v");
            System.out.println("//");
            System.out.println("//     +----------------------+  ^");
            System.out.println("//     |                      |  |");
            System.out.println("//     |                      |  u");
            System.out.println("//     |           o          |  |");
            System.out.println("//     |                      |   ");
            System.out.println("//     |                      |   ");
            System.out.println("//     +----------------------+   ");
            System.out.println("//               ------ v ---->   ");
            System.out.println("path3 myPlane(triple o, triple u, triple v)");
            System.out.println("{");
            System.out.println("    path3 P;");
            System.out.println("    P=o+(u+v)--o-u+v--o-u-v--o+u-v--cycle;");
            System.out.println("    return P;");
            System.out.println("}");
            System.out.println(" void drawAxes(triple o, triple u, triple v)");
            System.out.println(" {");
            System.out.println("pen p = red;");
            System.out.println(" draw(Label(\"$u$\",1),(o--o+u),p,Arrow3);");
            System.out.println("p=green;");
            System.out.println(" draw(Label(\"$v$\",1),(o--o+v),p,Arrow3);");
            System.out.println("p=blue;");
            System.out.println(" draw(Label(\"$w$\",1),(o--o+10.*cross(unit(u),unit(v))),p,Arrow3);");
            System.out.println(" }");

            // ECal crystal frustum volume
            // use volumes for hit crystals in event data
            System.out.println("void drawvolume(triple p1, triple p2, triple p3, triple p4, triple p5, triple p6, triple p7, triple p8, pen surfpen)");
            System.out.println("{");
            System.out.println("surface[] surfs;");
            System.out.println("");
            System.out.println("surfs.push( surface(p1--p2--p4--p3--cycle) );");
            System.out.println("surfs.push( surface(p4--p8--p7--p3--cycle) );");
            System.out.println("surfs.push( surface(p3--p7--p5--p1--cycle) );");
            System.out.println("surfs.push( surface(p2--p1--p5--p6--cycle) );");
            System.out.println("surfs.push( surface(p2--p6--p8--p4--cycle) );");
            System.out.println("surfs.push( surface(p6--p5--p7--p8--cycle) );");
            System.out.println("");
            System.out.println("for (int i = 0; i < surfs.length; ++i) {");
            System.out.println("        draw(surfs[i], surfpen);");
            System.out.println("    }");
            System.out.println("");
            System.out.println("}");

            //ECal crystal frustum edges
            // Use edges to draw detector
            // this causes out-of-memory crashes.
            // keep this here, but don't use it for now.
            System.out.println("void drawEdges(triple p1, triple p2, triple p3, triple p4, triple p5, triple p6, triple p7, triple p8, pen edgepen)");
            System.out.println("{");
            System.out.println("draw(p1--p2--p3--p4--cycle, edgepen);");
            System.out.println("draw(p2--p6--p7--p3--cycle, edgepen);");
            System.out.println("draw(p3--p7--p8--p4--cycle, edgepen);");
            System.out.println("draw(p1--p4--p8--p5--cycle, edgepen);");
            System.out.println("draw(p1--p5--p6--p2--cycle, edgepen);");
            System.out.println("draw(p5--p8--p7--p6--cycle, edgepen);");
            System.out.println("}");

            DecimalFormat df = new DecimalFormat("###.######");

            System.out.println("draw(Label(\"$x$\",1),(O--100.*X),Arrow3(HookHead3));");
            System.out.println("draw(Label(\"$y$\",1),(O--100.*Y),Arrow3(HookHead3));");
            System.out.println("draw(Label(\"$z$\",1),(O--1000.*Z),Arrow3(HookHead3));");

            for (int i = 0; i < oList.size(); ++i) {
                Hep3Vector o = oList.get(i);
                Hep3Vector u = uList.get(i);
                Hep3Vector v = vList.get(i);
                double l = unmeasDim.get(i) / 2.;
                double h = measDim.get(i) / 2.;
                //redefine op as lower left corner, viz.
                // o' = o-(l*V+h*u)
                // Hep3Vector op = VecOp.sub(o, VecOp.add(VecOp.mult(l, v), VecOp.mult(h, u)));
                String planeType = isAxial.get(i) ? "axial" : "stereo";
                System.out.println("//" + planeType);
                StringBuffer sb = new StringBuffer("draw(surface( myPlane( ");
                //origin
                sb.append("( " + df.format(o.x()) + " , " + df.format(o.y()) + " , " + df.format(o.z()) + " ), ");
                //u coordinate
                sb.append(h + "*( " + df.format(u.x()) + " , " + df.format(u.y()) + " , " + df.format(u.z()) + " ), ");
                //v coordinate
                sb.append(l + "*( " + df.format(v.x()) + " , " + df.format(v.y()) + " , " + df.format(v.z()) + " )  ) ), fill, " + planeType + "Edge); ");

                System.out.println(sb.toString());
                // sensor axes
                sb = new StringBuffer("drawAxes( ");
                //origin
                sb.append("( " + df.format(o.x()) + " , " + df.format(o.y()) + " , " + df.format(o.z()) + " ), ");
                //u coordinate
                sb.append(h + "*( " + df.format(u.x()) + " , " + df.format(u.y()) + " , " + df.format(u.z()) + " ), ");
                //v coordinate
                sb.append(l + "*( " + df.format(v.x()) + " , " + df.format(v.y()) + " , " + df.format(v.z()) + " ) ); ");

                System.out.println(sb.toString());
            }
            // now the calorimeter...
            StringBuffer sbcal = new StringBuffer();
            for (List<Hep3Vector> vertices : calcells) {
                addcell(sbcal, vertices);
            }
            System.out.println(sbcal.toString());
        }
    }

    @Override
    protected void process(EventHeader event) {
        boolean analyzeTracks = true;
        boolean analyzeStrips = true;
        if (analyzeStrips) {
            List<RawTrackerHit> rawTrackerHits = event.get(RawTrackerHit.class, "RawTrackerHitMaker_RawTrackerHits");
            EventHeader.LCMetaData meta = event.getMetaData(rawTrackerHits);
            // Setup links to SiSensor objects in detector geometry.
            setSensors(meta, (List<RawTrackerHit>) rawTrackerHits);

            List<SiTrackerHitStrip1D> stripHits = event.get(SiTrackerHitStrip1D.class, "StripClusterer_SiTrackerHitStrip1D");
            for (TrackerHit hit : stripHits) {
                SiTrackerHitStrip1D strip = new SiTrackerHitStrip1D(hit);
//                System.out.println("strip sensor id: "+strip.getSensor().getSensorID());
                //draw this strip
                if (_drawEvent) {
                    drawAsy(strip.getHitSegment());
                }
            }
        }
        if (analyzeTracks) {
            // get the tracks
            if (!event.hasCollection(Track.class, trackCollectionName)) {
                if (_debug) {
                    System.out.printf("%s: No tracks in Event %d \n", this.getClass().getSimpleName(), event.getEventNumber());
                }
                return;
            }

            //get the relations to the GBLtracks
            if (!event.hasItem(track2GblTrackRelationName)) {
                System.out.println("Need Relations " + track2GblTrackRelationName);
                return;
            }
            // and strips
            if (!event.hasItem(gblTrack2StripRelationName)) {
                System.out.println("Need Relations " + gblTrack2StripRelationName);
                return;
            }

            List<LCRelation> track2GblTrackRelations = event.get(LCRelation.class, track2GblTrackRelationName);
            //need a map of GBLTrackData keyed on the Generic object from which it created
            Map<GenericObject, GBLTrackData> gblObjMap = new HashMap<GenericObject, GBLTrackData>();
            //need a map of SeedTrack to GBLTrackData keyed on the track object from which it created
            Map<GBLTrackData, Track> gblToSeedMap = new HashMap<GBLTrackData, Track>();

            // loop over the relations
            for (LCRelation relation : track2GblTrackRelations) {
                Track t = (Track) relation.getFrom();
                GenericObject gblTrackObject = (GenericObject) relation.getTo();
                GBLTrackData gblT = new GBLTrackData(gblTrackObject);
                gblObjMap.put(gblTrackObject, gblT);
                gblToSeedMap.put(gblT, t);
            } //end of loop over tracks

            //get the strip hit relations
            List<LCRelation> gblTrack2StripRelations = event.get(LCRelation.class, gblTrack2StripRelationName);

            // need a map of lists of strip data keyed by the gblTrack to which they correspond
            Map<GBLTrackData, List<GBLStripClusterData>> stripsGblMap = new HashMap<GBLTrackData, List<GBLStripClusterData>>();
            for (LCRelation relation : gblTrack2StripRelations) {
                //from GBLTrackData to GBLStripClusterData
                GenericObject gblTrackObject = (GenericObject) relation.getFrom();
                //Let's get the GBLTrackData that corresponds to this object...
                GBLTrackData gblT = gblObjMap.get(gblTrackObject);
                GBLStripClusterData sd = new GBLStripClusterData((GenericObject) relation.getTo());
                if (stripsGblMap.containsKey(gblT)) {
                    stripsGblMap.get(gblT).add(sd);
                } else {
                    stripsGblMap.put(gblT, new ArrayList<GBLStripClusterData>());
                    stripsGblMap.get(gblT).add(sd);
                }
            }
//            System.out.println("found " + stripsGblMap.size() + " GBL Tracks");
            int nTrack = 0;
            for (GBLTrackData t : stripsGblMap.keySet()) {
                List<GBLStripClusterData> hits = stripsGblMap.get(t);
                // Loop over strips
                int n_strips = hits.size();
                if (_debug) {
                    System.out.println("track " + nTrack + " has " + n_strips + " strip hits");
                }
                for (int istrip = 0; istrip != n_strips; ++istrip) {
                    GBLStripClusterData strip = hits.get(istrip);
                    if (_debug) {
                        System.out.println("HpsGblFitter: Processing strip " + istrip + " with id/layer " + strip.getId());
                    }
                    Hep3Vector u = strip.getU();
                    Hep3Vector v = strip.getV();
                    Hep3Vector w = strip.getW();
                    if (_debug) {
                        System.out.println("u: " + u);
                        System.out.println("v: " + v);
                        System.out.println("w: " + w);
                    }
                }// end of loop over strips
                nTrack++;
            }//end of loop over tracks
        }
    }

    private void addcell(StringBuffer sbcal, List<Hep3Vector> vertices) {
        sbcal.append("drawvolume( ");
        for (Hep3Vector v : vertices) {
            sbcal.append("( " + v.x() + ", " + v.y() + ", " + v.z() + "),");
        }
        sbcal.append(" gray(0.2)+opacity(0.1));" + "\n");
    }

    private void drawAsy(LineSegment3D line) {
        Point3D A = line.getStartPoint();
        Point3D B = line.getEndPoint();
        //draw this strip
        drawAsy(A, B);
    }

    private void drawAsy(Point3D A, Point3D B) {
        //draw((0,0,0)--(0,0,1));
        System.out.println("draw ((" + A.x() + "," + A.y() + "," + A.z() + ")--(" + B.x() + "," + B.y() + "," + B.z() + "));");
    }

    private void setSensors(EventHeader.LCMetaData meta, List<RawTrackerHit> hits) {
        // Get the ID dictionary and field information.
        IIdentifierDictionary dict = meta.getIDDecoder().getSubdetector().getDetectorElement().getIdentifierHelper().getIdentifierDictionary();
        int fieldIdx = dict.getFieldIndex("side");
        int sideIdx = dict.getFieldIndex("strip");

        for (RawTrackerHit hit : hits) {
            // The "side" and "strip" fields needs to be stripped from the ID for sensor lookup.
            IExpandedIdentifier expId = dict.unpack(hit.getIdentifier());
            expId.setValue(fieldIdx, 0);
            expId.setValue(sideIdx, 0);
            IIdentifier strippedId = dict.pack(expId);

            // Find the sensor DetectorElement.
            List<IDetectorElement> des = DetectorElementStore.getInstance().find(strippedId);
            if (des == null || des.size() == 0) {
                throw new RuntimeException("Failed to find any DetectorElements with stripped ID <0x" + Long.toHexString(strippedId.getValue()) + ">.");
            } else if (des.size() == 1) {
                hit.setDetectorElement((SiSensor) des.get(0));
            } else {
                // Use first sensor found, which should work unless there are sensors with duplicate IDs.
                for (IDetectorElement de : des) {
                    if (de instanceof SiSensor) {
                        hit.setDetectorElement((SiSensor) de);
                        break;
                    }
                }
            }

            // No sensor was found.
            if (hit.getDetectorElement() == null) {
                throw new RuntimeException("No sensor was found for hit with stripped ID <0x" + Long.toHexString(strippedId.getValue()) + ">.");
            }
        }
    }

    public void setDrawEvent(boolean drawEvent) {
        _drawEvent = drawEvent;
    }

    public void setdrawDetector(boolean drawDetector) {
        _drawDetector = drawDetector;
    }

    public void setDebug(boolean debug) {
        _debug = debug;
    }
}
