package org.hps.recon.filtering;

import java.util.ArrayList;
import java.util.List;


//===> import org.hps.conditions.deprecated.SvtUtils;
import org.lcsim.detector.identifier.IIdentifier;
import org.lcsim.detector.tracker.silicon.ChargeCarrier;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.detector.tracker.silicon.SiSensor;
import org.lcsim.detector.tracker.silicon.SiSensorElectrodes;
import org.lcsim.detector.tracker.silicon.SiTrackerIdentifierHelper;
import org.lcsim.event.EventHeader;
import org.lcsim.event.RawTrackerHit;

/**
 * 
 * @author Sebouh Paul
 */
public class EventRangeFilter extends EventReconFilter{
  
    
    boolean inBadRange(int run, int event){
       for(EventRange range : badRanges){
           if(run == range.run && event >= range.min && event <= range.max)
               return false;
       }
       return true;
    }
    
    @Override
    public void process(EventHeader event){
        incrementEventProcessed();
        if(inBadRange(event.getRunNumber(), event.getEventNumber()))
            skipEvent();

        
        incrementEventPassed();
    }
    
    ArrayList<EventRange> badRanges = new ArrayList();

    
    public void setBadRanges(String s){
        String split[] = s.split("[\n\r ,]+");
        int n = split.length/3;
        for(int i = 0; i< n; i++){
            int run = Integer.parseInt(split[i*3]);
            long min = Long.parseLong(split[i*3]+1);
            long max = Long.parseLong(split[i*3]+2);
            badRanges.add(new EventRange(run, min, max));
        }
    }
    class EventRange {
        int run; long min, max;
        EventRange(int run, long min, long max){
            this.run = run; this.min = min; this.max = max;
        }
    }
}
