package org.hps.recon.ecal;

public class TimeDependentGainTweakDriver2015 extends TimeDependentGainTweakDriver {
    TimeDependentGainData data = new TimeDependentGains2015();
    
    @Override
    double getCorrFactor(int ix, int iy, double timestamp) {
        
        return data.getCorrectionFactor(ix, iy, timestamp);
    }

}
