package org.hps.recon.tracking.kalman;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.hps.recon.tracking.MaterialSupervisor;
import org.hps.recon.tracking.MaterialSupervisor.ScatteringDetectorVolume;
import org.hps.recon.tracking.MaterialSupervisor.SiStripPlane;
import org.lcsim.event.EventHeader;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;

public class KalmanDriverHPS extends Driver {

    private ArrayList<SiModule> testData;
    private ArrayList<SiStripPlane> detPlanes;
    private MaterialSupervisor _materialManager;
    private FieldMap fm;
    private String fieldMapFileName = "fieldmap/125acm2_3kg_corrected_unfolded_scaled_0.7992.dat";
    private String trackCollectionName = "GBLTracks";
    private KalmanInterface KI;

    public void setMaterialManager(MaterialSupervisor mm) {
        _materialManager = mm;
    }

    public void setTrackCollectionName(String input) {
        trackCollectionName = input;
    }

    public void setFieldMapFilename(String input) {
        fieldMapFileName = input;
    }

    public void setTestData(ArrayList<SiModule> input) {
        testData = input;
    }

    public void constructTestData() {
        double[] location = { 100., 200., 300., 500., 700., 900. };
        double delta = 5.0;
        double[] heights = { 100., 100., 100., 100., 100., 100. };
        double[] widths = { 150., 150., 150., 300., 300., 300. };
        double[] stereoAngle = { 0.1, 0.1, 0.1, 0.05, 0.05, 0.05 };
        double thickness = 0.3;
        Vec tInt = new Vec(0., 1., 0.);

        testData = new ArrayList<SiModule>(12);
        for (int pln = 0; pln < 6; pln++) {
            Vec rInt1 = new Vec(0., location[pln], 0.);

            Plane pInt1 = new Plane(rInt1, tInt);
            SiModule newModule1 = new SiModule(pln, pInt1, 0., widths[pln], heights[pln], thickness, fm);
            testData.add(newModule1);

            Vec rInt2 = new Vec(0., location[pln] + delta, 0.);
            Plane pInt2 = new Plane(rInt2, tInt);
            SiModule newModule2 = new SiModule(pln, pInt2, stereoAngle[pln], widths[pln], heights[pln], thickness, fm);
            testData.add(newModule2);
        }

    }

    @Override
    public void detectorChanged(Detector det) {
        _materialManager = new MaterialSupervisor();
        _materialManager.buildModel(det);
        KI = new KalmanInterface();
        detPlanes = new ArrayList<SiStripPlane>();
        List<ScatteringDetectorVolume> materialVols = ((MaterialSupervisor) (_materialManager)).getMaterialVolumes();
        for (ScatteringDetectorVolume vol : materialVols) {
            detPlanes.add((SiStripPlane) (vol));
        }
        try {
            fm = new FieldMap(fieldMapFileName, "HPS", 0, 0, 0);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void process(EventHeader event) {
        if (!event.hasCollection(Track.class, trackCollectionName)) {
            System.out.println(trackCollectionName + " does not exist; skipping event");
            return;
        }
        List<Track> tracks = event.get(Track.class, trackCollectionName);
        constructTestData();
        //        for (SiModule SiM : testData) {
        //            SiM.print("SiModFromTestData");
        //        }

        ArrayList<SiModule> SiMlist = testSiModuleCreation();
        //        for (SiModule SiM : SiMlist)
        //            SiM.print("SiModFromHPS");

        for (Track trk : tracks) {
            KI.fillMeasurements(SiMlist, trk.getTrackerHits());
        }
        for (SiModule SiM : SiMlist)
            SiM.print("SiModFilled");
        //
        //        SeedTrack testKalmanTrack = new SeedTrack(testData, 0, 1, 12, false);
        //        testKalmanTrack.print("testKalmanTrack");
        //        Track HPStrk = KI.createTrack(testKalmanTrack, getTestMeasurements());
        //        printTrackInfo(HPStrk);
    }

    private ArrayList<SiModule> testSiModuleCreation() {
        ArrayList<SiModule> SiMods = new ArrayList<SiModule>();
        for (SiStripPlane ssp : detPlanes) {
            SiModule SiM = KI.createSiModule(ssp, fm);
            SiMods.add(SiM);
        }
        return SiMods;
    }

    private void printTrackInfo(Track HPStrk) {
        TrackState ts = HPStrk.getTrackStates().get(0);
        double[] params = ts.getParameters();
        System.out.printf("Track hits: %d \n", HPStrk.getTrackerHits().size());
        System.out.printf("      params: %f %f %f %f %f \n", params[0], params[1], params[3], params[4], params[5]);
    }

    public ArrayList<Measurement> getTestMeasurements() {
        ArrayList<Measurement> measList = new ArrayList<Measurement>();
        for (SiModule SiM : testData) {
            measList.addAll(SiM.hits);
        }
        return measList;
    }

}
