package org.hps.recon.tracking;

import java.util.List;
import java.util.ArrayList;

import org.lcsim.event.GenericObject;

/**
 * @author Omar Moreno <omoreno1@ucsc.edu>
 */
public class HitTimeData implements GenericObject {

    List<Double> hitTimeError = new ArrayList<Double>();

    public HitTimeData(double input) {
        hitTimeError.add(input);
    }

    public double getHitTimeError() {
        if (hitTimeError.size() == 1)
            return hitTimeError.get(0);
        return 0;
    }

    public void setHitTimeError(double input) {
        hitTimeError.clear();
        hitTimeError.add(input);
    }

    @Override
    public int getNDouble() {
        return hitTimeError.size();
    }

    /**
     * 
     */
    @Override
    public boolean isFixedSize() {
        return true;
    }

    @Override
    public int getNInt() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int getNFloat() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int getIntVal(int index) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public float getFloatVal(int index) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public double getDoubleVal(int index) {
        // TODO Auto-generated method stub
        return 0;
    }
}
