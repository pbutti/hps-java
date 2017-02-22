package org.hps.recon.ecal.tdeg;


public class TdegTweakDriver2015 extends TdegTweakDriver {
    TdegData data = new Tdeg2015();
    
    @Override
    double getCorrFactor(int ix, int iy, double timestamp) {
        data.beamEnergy = 1.056;
        return data.getCorrectionFactor(ix, iy, timestamp);
    }

}
