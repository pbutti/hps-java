package org.hps.recon.tracking;

import java.util.List;

import org.hps.util.Pair;

public class TimedSiliconResolutionModel extends DefaultSiliconResolutionModel {
    @Override
    public Pair<Double, Double> getTime(List<FittedRawTrackerHit> cluster) {
        double time_sum = 0;
        double signal_sum = 0;
        double err_sum = 0;

        for (FittedRawTrackerHit hit : cluster) {

            double signal = hit.getAmp();
            double time = hit.getT0();
            time_sum += time * signal * signal;
            signal_sum += signal * signal;
            err_sum += Math.pow(hit.getT0Err(), 2) * Math.pow(signal, 4);
        }
        return new Pair<Double, Double>(time_sum / signal_sum, Math.sqrt(err_sum) / signal_sum);
    }
}
