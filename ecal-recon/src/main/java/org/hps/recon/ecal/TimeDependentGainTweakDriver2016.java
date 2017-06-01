package org.hps.recon.ecal;

public class TimeDependentGainTweakDriver2016 extends TimeDependentGainTweakDriver {
    TimeDependentGainData data = new TimeDependentGains2016Constant();
    
    @Override
    double getCorrFactor(int ix, int iy, double timestamp) {
        
        return data.getCorrectionFactor(ix, iy, timestamp);
    }

}
