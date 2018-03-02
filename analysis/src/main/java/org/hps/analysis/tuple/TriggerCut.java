package org.hps.analysis.tuple;

import org.hps.record.triggerbank.AbstractIntData;
import org.hps.record.triggerbank.TIData;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;

class TriggerCut extends TupleCut {

    public final boolean passes(EventHeader event) {
        String triggerType = getParameter("triggerType");        
        TIData triggerData = getTriggerData(event);
        if (triggerData != null) {
            if (triggerType.contentEquals("") || triggerType.contentEquals("all")) {
                return true;
            }
            if (triggerData.isSingle0Trigger() && triggerType.contentEquals("singles0")) {
                return true;
            }
            if (triggerData.isSingle1Trigger() && triggerType.contentEquals("singles1")) {
                return true;
            }
            if (triggerData.isPair0Trigger() && triggerType.contentEquals("pairs0")) {
                return true;
            }
            if (triggerData.isPair1Trigger() && triggerType.contentEquals("pairs1")) {
                return true;
            }
        }
        return false;
    }
    
    private TIData getTriggerData(EventHeader event) {
        TIData triggerData = null;
        if (event.hasCollection(GenericObject.class, "TriggerBank")) {
            for (GenericObject data : event.get(GenericObject.class, "TriggerBank")) {
                if (AbstractIntData.getTag(data) == TIData.BANK_TAG) {
                    triggerData = new TIData(data);
                }
            }
        }
        return triggerData;
    }
}
