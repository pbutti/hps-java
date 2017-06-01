package org.hps.recon.ecal;


/**
 * This is the version to use during recon passes, to correct the energy of the individual hits (rather than post-recon corrections of the cluster energy).  
 * @author spaul
 *
 */
public class TimeDependentGainDriver2016 extends TimeDependentGainDriver{
    
    public TimeDependentGainDriver2016(){
        this.tdeg = new TimeDependentGains2016();
    }
    
    
}
