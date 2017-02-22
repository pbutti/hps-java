package org.hps.recon.ecal.tdeg;

import java.util.HashMap;
import java.util.Map;


public class TdegData {
    public TdegData(){
        setup();
    }
    protected void setup(){
        
    }
    Map<Integer, TdegFunc> funcs = new HashMap();
    
    
    double beamEnergy;
    
    private void addCrystal0(int ix, int iy, double ...params){
        //System.out.println(ix + " " +  iy);
        funcs.put(encode(ix,iy), new TdegFunc(params));
    }
    /**
     * sets the tdeg parameters for an individual crystals with at (ix, iy).  
     * If iy == +- 2, this also sets the parameters for the crystal with iy == +-1 in the same column.
     * Likewise, if iy == +- 4, this also sets the parameters for the crystal with iy == +-5 in the same column.
     * @param ix
     * @param iy
     * @param params
     */
    void addCrystal(int ix, int iy, double ... params){
        addCrystal0(ix,iy, params);
        if(iy == 2 || iy == -4){
            addCrystal0(ix, iy-1, params);
        }
        if(iy == -2 || iy == 4){
            addCrystal0(ix, iy+1, params);
        }    
    }
    
    void setGlobal(double...params){
        TdegFunc func = new TdegFunc(params);
        for(int ix = -23; ix<= 23; ix ++){
            for (int iy = -5; iy <=5; iy++){
                if(iy*ix == 0)
                    continue;
                funcs.put(encode(ix,iy), func);
            }
        }
    }
    /**
     * sets the tdeg parameters for all crystals with ix >= ixMin 
     * @param ixMin
     * @param params
     */
    void farPositronSide(int ixMin, double ... params){
        TdegFunc func = new TdegFunc(params);
        for(int ix = ixMin; ix<= 23; ix ++){
            for (int iy = -5; iy <=5; iy++){
                if(iy*ix == 0)
                    continue;
                funcs.put(encode(ix,iy), func);
            }
        }
    }
    /**
     * sets the tdeg parameters for all crystals with ix <= ixMax
     * @param ixMax
     * @param params
     */
    void farElectronSide(int ixMax, double ... params){
        TdegFunc func = new TdegFunc(params);
        for(int ix = -23; ix<= ixMax; ix ++){
            for (int iy = -5; iy <=5; iy++){
                if(iy*ix == 0)
                    continue;
                funcs.put(encode(ix,iy), func);
            }
        }
    }
    
    double getCorrectionFactor(int ix, int iy, double timestamp){
        TdegFunc func = funcs.get(encode(ix,iy));
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
