package org.hps.recon.ecal.tdeg;

import java.util.ArrayList;
import java.util.List;

import org.hps.conditions.beam.BeamEnergy.BeamEnergyCollection;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.EventHeader;
import org.lcsim.event.base.BaseCalorimeterHit;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
/**
 * This is the version to use during recon passes, to correct the energy of the individual hits (rather than post-recon corrections of the cluster energy).  
 * @author spaul
 *
 */
public class TdegGainDriver extends Driver{
    String inputHitsCollectionName;
    String outputHitsCollectionName;
    
    TdegData tdeg = new Tdeg2016();
    
    @Override
    public void detectorChanged(Detector det){
       BeamEnergyCollection beamEnergyCollection = this.getConditionsManager().getCachedConditions(BeamEnergyCollection.class, "beam_energies").getCachedData(); 
        tdeg.beamEnergy = beamEnergyCollection.get(0).getBeamEnergy();
    }
    
    @Override
    public void process(EventHeader event){
        List<CalorimeterHit> hits = event.get(CalorimeterHit.class, inputHitsCollectionName);
        int flags = event.getMetaData(hits).getFlags();
        double timestamp = event.getTimeStamp();
        List<CalorimeterHit> outhits = new ArrayList();
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
}
