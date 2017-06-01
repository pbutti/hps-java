package org.hps.recon.ecal;

import java.util.ArrayList;
import java.util.List;

import org.hps.conditions.beam.BeamEnergy.BeamEnergyCollection;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.EventHeader;
import org.lcsim.event.base.BaseCalorimeterHit;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;

public class TimeDependentGainDriver extends Driver{
    private String inputHitsCollectionName;
    private String outputHitsCollectionName;
    
    TimeDependentGainParameters parameterLists[] = {new TimeDependentGainParameters2015(),  new TimeDependentGainParameters2016()};
    
    
    
    @Override
    public void process(EventHeader event){
        List<CalorimeterHit> hits = event.get(CalorimeterHit.class, inputHitsCollectionName);
        int flags = event.getMetaData(hits).getFlags();
        double timestamp = event.getTimeStamp()/1e9; //convert to ns.  
        List<CalorimeterHit> outhits = new ArrayList();
        
        //determine whether to use the parameters for the 2016 range or the 2015 range.  
        TimeDependentGainParameters tdeg = timestamp > 1457140000 ? parameterLists[1] : parameterLists[0];
        
        for(CalorimeterHit inhit : hits){
            double corrFactor = tdeg.getCorrectionFactor(inhit.getIdentifierFieldValue("ix"), inhit.getIdentifierFieldValue("iy"), timestamp);
            BaseCalorimeterHit corrhit = new BaseCalorimeterHit(
                    inhit.getRawEnergy()*corrFactor, 
                    inhit.getCorrectedEnergy()*corrFactor,
                    inhit.getEnergyError()*corrFactor,
                    inhit.getTime(),
                    inhit.getCellID(), 
                    inhit.getPositionVec(), 
                    inhit.getType(),
                    inhit.getMetaData()
            );
            outhits.add(corrhit);
        }
        event.put(outputHitsCollectionName, outhits, CalorimeterHit.class, flags);
    }
    
    public void setInputHitsCollectionName(String val){
        this.inputHitsCollectionName = val;
    }
    public void setOutputHitsCollectionName(String val){
        this.outputHitsCollectionName = val;
    }
}
