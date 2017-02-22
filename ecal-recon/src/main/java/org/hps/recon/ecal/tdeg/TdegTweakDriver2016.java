package org.hps.recon.ecal.tdeg;


public class TdegTweakDriver2016 extends TdegTweakDriver {
    TdegData data = new Tdeg2016Const();
    
    @Override
    double getCorrFactor(int ix, int iy, double timestamp) {
        data.beamEnergy = 2.306;
        return data.getCorrectionFactor(ix, iy, timestamp);
    }

}
