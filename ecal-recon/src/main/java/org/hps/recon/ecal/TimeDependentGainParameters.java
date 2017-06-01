package org.hps.recon.ecal;

import java.util.HashMap;
import java.util.Map;

/**
 * Object containing the parameters for time dependent ecal gains, as mapped from individual crystals,
 * ranges of crystals, or the entire Ecal.  
 * @author spaul
 *
 */
public class TimeDependentGainParameters {
    
    public TimeDependentGainParameters(){
        setup();
    }
    protected void setup(){
        
    }
    Map<Integer, TimeDependentGainFunc> funcs = new HashMap();
    
    
    double beamEnergy;
    
    private void addCrystal0(int ix, int iy,  TimeDependentGainFunc func){
        //System.out.println(ix + " " +  iy);
        funcs.put(encode(ix,iy),func);
    }
    /**
     * sets the tdeg parameters for an individual crystals with at (ix, iy).  
     * If iy == +- 2, this also sets the parameters for the crystal with iy == +-1 in the same column.
     * Likewise, if iy == +- 4, this also sets the parameters for the crystal with iy == +-5 in the same column.
     * @param ix
     * @param iy
     * @param params
     */
    protected void addCrystal(int ix, int iy, double ... params){
        TimeDependentGainFunc func = new TimeDependentGainFunc(params);
        addCrystal0(ix,iy, func);
        if(iy == 2 || iy == -4){
            addCrystal0(ix, iy-1, func);
        }
        if(iy == -2 || iy == 4){
            addCrystal0(ix, iy+1, func);
        }    
    }
    
    protected void addRange(int ix1, int iy1, int ix2, int iy2, double ...params){
        TimeDependentGainFunc func = new TimeDependentGainFunc(params);
        for(int ix = ix1; ix<=ix2; ix++){
            for(int iy = iy1; iy<=iy2; iy++){
                addCrystal0(ix, iy, func);
            }
        }
    }
    
    protected void setGlobal(double...params){
        TimeDependentGainFunc func = new TimeDependentGainFunc(params);
        for(int ix = -23; ix<= 23; ix ++){
            for (int iy = -5; iy <=5; iy++){
                if(iy*ix == 0)
                    continue;
                addCrystal0(ix, iy, func);
            }
        }
    }
    /**
     * sets the tdeg parameters for all crystals with ix >= ixMin 
     * @param ixMin
     * @param params
     */
    protected void farPositronSide(int ixMin, double ... params){
        TimeDependentGainFunc func = new TimeDependentGainFunc(params);
        for(int ix = ixMin; ix<= 23; ix ++){
            for (int iy = -5; iy <=5; iy++){
                if(iy*ix == 0)
                    continue;
                addCrystal0(ix,iy, func);
            }
        }
    }
    /**
     * sets the tdeg parameters for all crystals with ix <= ixMax
     * @param ixMax
     * @param params
     */
    protected void farElectronSide(int ixMax, double ... params){
        TimeDependentGainFunc func = new TimeDependentGainFunc(params);
        for(int ix = -23; ix<= ixMax; ix ++){
            for (int iy = -5; iy <=5; iy++){
                if(iy*ix == 0)
                    continue;
                addCrystal0(ix,iy, func);
            }
        }
    }
    
    public double getCorrectionFactor(int ix, int iy, double timestamp){
        TimeDependentGainFunc func = funcs.get(encode(ix,iy));
        if(func == null){
            System.err.println("function undefined for crystal (" + ix + ", " + iy + ")");
            return 1;
        }
        return beamEnergy/func.getFittedFunction(timestamp);
    }
    
    private int encode(int ix, int iy){
        return (ix+23)*10+(iy+5) + (iy>0 ? 0 :1);
    }
    
}
