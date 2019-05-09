/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.hps.recon.hodo;

import org.lcsim.event.EventHeader;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;

/**
 *
 * @author rafopar
 */
public class testDriver extends Driver {

    @Override
    public void startOfData() {

        System.out.println("======== Start of Data ======");
    }

    @Override
    public void detectorChanged(Detector detector) {
        System.out.println("======== detectorChanged ======");
    }  
    
    
    @Override
    public void process(EventHeader event) {
        System.out.println("======== process ======");
    }    
    
}
