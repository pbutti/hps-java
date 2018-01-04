package org.hps.recon.filtering;

import java.util.ArrayList;
import java.util.Scanner;

import org.lcsim.event.EventHeader;

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
        
        Scanner scanner = new Scanner(s);
        scanner.useDelimiter("[\n\r\t ,]+");
        while (scanner.hasNextLong()){
            int run = scanner.nextInt();
            long min = scanner.nextLong();
            long max = scanner.nextLong();
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
