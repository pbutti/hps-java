/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.hps.evio;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.conditions.hodoscope.HodoscopeChannel;
import org.hps.conditions.hodoscope.HodoscopeChannel.GeometryId;
import org.hps.conditions.hodoscope.HodoscopeChannel.HodoscopeChannelCollection;
import org.hps.recon.ecal.FADCGenericHit;
import org.hps.recon.ecal.HitExtraData;
import org.hps.recon.ecal.HitExtraData.Mode7Data;
import org.jlab.coda.jevio.BaseStructure;
import org.jlab.coda.jevio.BaseStructureHeader;
import org.jlab.coda.jevio.CompositeData;
import org.jlab.coda.jevio.EvioEvent;
import org.jlab.coda.jevio.EvioException;
import org.lcsim.detector.IDetectorElementContainer;
import org.lcsim.detector.identifier.IIdentifierHelper;
import org.lcsim.detector.identifier.Identifier;
import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.event.RawCalorimeterHit;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.base.BaseLCRelation;
import org.lcsim.event.base.BaseRawCalorimeterHit;
import org.lcsim.event.base.BaseRawTrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.Subdetector;
import org.lcsim.lcio.LCIOConstants;

/**
 *
 * This class is similar to the ECalEvioReader.  It was copied and modified
 * to work with Hodoscope data.
 * 
 * @author rafopar
 */
public class HodoEvioReader extends EvioReader {

    private int bankTag = 0;
    private Class<?> hitClass = BaseRawCalorimeterHit.class;

    private HodoscopeChannelCollection hodoChannels = null;
    private IIdentifierHelper helper = null;

    private static final String readoutName = "HodoHits";
    private static final String subdetectorName = "Hodoscope";

    private Subdetector subDetector;

    private static final String genericHitCollectionName = "FADCGenericHits";
    private List<FADCGenericHit> genericHits;

    private static final String extraDataRelationsName = "HodoReadoutExtraDataRelations";
    private List<LCRelation> extraDataRelations;

    private static final String extraDataCollectionName = "HodoReadoutExtraData";
    private List<HitExtraData> extraDataList;

    private final Map<List<Integer>, Integer> genericHitCount = new HashMap<List<Integer>, Integer>();

    private static final Logger LOGGER = Logger.getLogger(HodoEvioReader.class.getPackage().getName());

    private int topBankTag, botBankTag;

    public HodoEvioReader(int topBankTag, int botBankTag) {
        this.topBankTag = topBankTag;
        this.botBankTag = botBankTag;
        this.hitCollectionName = "HodoReadoutHits";
    }

    public void setTopBankTag(int topBankTag) {
        this.topBankTag = topBankTag;
    }

    public void setBotBankTag(int botBankTag) {
        this.botBankTag = botBankTag;
    }

    @Override
    public boolean makeHits(EvioEvent event, EventHeader lcsimEvent) {        
        boolean foundHits = false;
        List<Object> hits = new ArrayList<Object>();
        genericHits = new ArrayList<FADCGenericHit>();
        extraDataList = new ArrayList<HitExtraData>();
        extraDataRelations = new ArrayList<LCRelation>();
        int flags = 0;
        for (BaseStructure bank : event.getChildrenList()) {
            BaseStructureHeader header = bank.getHeader();
            int crateBankTag = header.getTag();
            int crate = 0;
            if (crateBankTag == topBankTag) {
                crate = 1;
            } else if (crateBankTag == botBankTag) {
                crate = 2;
            }
            if (crateBankTag == topBankTag || crateBankTag == botBankTag) {
                foundHits = true;
                if (bank.getChildCount() > 0) {
                    if (debug) {
                        System.out.println("Hodo bank tag: " + header.getTag() + "; childCount: " + bank.getChildCount());
                    }
                    try {                                               
                        for (BaseStructure slotBank : bank.getChildrenList()) {
                            if (slotBank.getCompositeData() != null) { //skip SSP and TI banks, if any
                                for (CompositeData cdata : slotBank.getCompositeData()) {
//                            CompositeData cdata = slotBank.getCompositeData();
                                    // For some reason, if the ECal is readout earlier, the index get gets some non 0 value.
                                    cdata.index(0);
//                                    System.out.println("cdata index is " + cdata.index());                                    
                                    if (slotBank.getHeader().getTag() != bankTag) {
                                        bankTag = slotBank.getHeader().getTag();
                                        LOGGER.info(String.format("Hodo format tag: 0x%x\n", bankTag));
                                    }
                                    switch (slotBank.getHeader().getTag()) {
                                        case EventConstants.FADC_MODE1_BANK_TAG:
                                            hits.addAll(makeWindowHits(cdata, crate));
                                            hitClass = RawTrackerHit.class;
                                            flags = 0;
                                            break;
                                        case EventConstants.FADC_PULSE_BANK_TAG:
                                            hits.addAll(makePulseHits(cdata, crate));
                                            hitClass = RawTrackerHit.class;
                                            flags = 0;
                                            break;
                                        case EventConstants.FADC_PULSE_INTEGRAL_BANK_TAG:
                                            hits.addAll(makeIntegralHitsMode3(cdata, crate));
                                            hitClass = RawCalorimeterHit.class;
                                            flags = (1 << LCIOConstants.RCHBIT_TIME); //store timestamp
                                            break;
                                        case EventConstants.FADC_PULSE_INTEGRAL_HIGHRESTDC_BANK_TAG:
                                            hits.addAll(makeIntegralHitsMode7(cdata, crate));
                                            hitClass = RawCalorimeterHit.class;
                                            flags = (1 << LCIOConstants.RCHBIT_TIME); //store timestamp
                                            break;
                                        default:
                                            throw new RuntimeException("Unsupported Hodo format - bank tag " + slotBank.getHeader().getTag());
                                    }
                                }
                            }
                        }
                    } catch (EvioException e) {
                        throw new RuntimeException(e);
                    }
                }
            }
        }

        lcsimEvent.put(this.hitCollectionName, hits, hitClass, flags, readoutName);
        lcsimEvent.put(genericHitCollectionName, genericHits, FADCGenericHit.class, 0);
        if (!extraDataList.isEmpty()) {
            lcsimEvent.put(extraDataCollectionName, extraDataList, Mode7Data.class, 0);
            lcsimEvent.put(extraDataRelationsName, extraDataRelations, LCRelation.class, 0);
        }

        return foundHits;
    }

    private BaseRawTrackerHit makeHodoRawHit(int time, long id, CompositeData cdata, int nSamples) {

        short[] adcValues = new short[nSamples];
        for (int i = 0; i < nSamples; i++) {
            adcValues[i] = cdata.getShort();
            //System.out.println("ADC["+i+"] = " + adcValues[i]);
        }
        IDetectorElementContainer srch = subDetector.getDetectorElement().findDetectorElement(new Identifier(id));
        if (srch.size() == 0) {
            throw new RuntimeException("No detector element was found for hit ID: " + helper.unpack(new Identifier(id)));
        }
        return new BaseRawTrackerHit( 
                time,
                id,
                adcValues,
                new ArrayList<SimTrackerHit>(),
                srch.get(0)
                );
    }

    private static FADCGenericHit makeGenericRawHit(int mode, int crate, short slot, short channel, CompositeData cdata, int nSamples) {
        int[] adcValues = new int[nSamples];
        for (int i = 0; i < nSamples; i++) {
            adcValues[i] = cdata.getShort();
        }
        return new FADCGenericHit(mode, crate, slot, channel, adcValues);
    }

    private List<BaseRawTrackerHit> makeWindowHits(CompositeData cdata, int crate) {
        List<BaseRawTrackerHit> hits = new ArrayList<BaseRawTrackerHit>();
//        if (debug) {
//            int n = cdata.getNValues().size();
//            for (int i = 0; i < n; i++) {
//                System.out.println("cdata.N[" + i + "]=" + cdata.getNValues().get(i));
//            }
//            int ni = cdata.getItems().size();
//            for (int i = 0; i < ni; i++) {
//                System.out.println("cdata.type[" + i + "]=" + cdata.getTypes().get(i));
//            }
//        }   
        while (cdata.index() + 1 < cdata.getItems().size()) {
            short slot = cdata.getByte();
            int trigger = cdata.getInt();
            long timestamp = cdata.getLong();
            int nchannels = cdata.getNValue();            
            if (debug) {
                System.out.println("slot#=" + slot + "; trigger=" + trigger + "; timestamp=" + timestamp + "; nchannels=" + nchannels);
            }
            for (int j = 0; j < nchannels; j++) {
                short channel = cdata.getByte();
                int nSamples = cdata.getNValue();

                if (debug) {
                    System.out.println("  channel=" + channel + "; nSamples=" + nSamples);
                }

                Long id = daqToGeometryId(crate, slot, channel);

                if (debug) {
                    System.out.println("The long id is: " + id);
                }

                if (id == null) {
                    FADCGenericHit hit = makeGenericRawHit(EventConstants.HODO_RAW_MODE, crate, slot, channel, cdata, nSamples);
                    processUnrecognizedChannel(hit);
                } else {
                    BaseRawTrackerHit hit = makeHodoRawHit(0, id, cdata, nSamples);
                    hits.add(hit);
                }
            }
        }
        return hits;
    }

    private Long daqToGeometryId(int crate, short slot, short channel) {

//        HodoscopeChannel hodoChannel = hodoConditions.getChannels().findChannel(crate, 10, channel);
        HodoscopeChannel hodoChannel = hodoChannels.findChannel(crate, slot, channel);

        if (hodoChannel == null) {
            return null;
        }
        int ix = hodoChannel.getIX();
        int iy = hodoChannel.getIY();
        int ilayer = hodoChannel.getLayer();
        //int ihole = hodoChannel.getHole();
        
        /* The 'hole' ID value is set to zero to match with the pixel detector elements in the geometry. --JM */
        int ihole = 0; 
        GeometryId geometryId = new GeometryId(helper, new int[]{subDetector.getSystemID(), ix, iy, ilayer, ihole});
        Long id = geometryId.encode();

        return id;
    }

    private List<BaseRawTrackerHit> makePulseHits(CompositeData cdata, int crate) {
        List<BaseRawTrackerHit> hits = new ArrayList<BaseRawTrackerHit>();
        if (debug) {
            int n = cdata.getNValues().size();
            for (int i = 0; i < n; i++) {
                System.out.println("cdata.N[" + i + "]=" + cdata.getNValues().get(i));
            }
            int ni = cdata.getItems().size();
            for (int i = 0; i < ni; i++) {
                System.out.println("cdata.type[" + i + "]=" + cdata.getTypes().get(i));
            }
        }
        while (cdata.index() + 1 < cdata.getItems().size()) { //the +1 is a hack because sometimes an extra byte gets read (padding)
            short slot = cdata.getByte();
            int trigger = cdata.getInt();
            long timestamp = cdata.getLong();
            int nchannels = cdata.getNValue();
            if (debug) {
                System.out.println("slot#=" + slot + "; trigger=" + trigger + "; timestamp=" + timestamp + "; nchannels=" + nchannels);
            }
            for (int j = 0; j < nchannels; j++) {
                short channel = cdata.getByte();
                int npulses = cdata.getNValue();
                if (debug) {
                    System.out.println("  channel=" + channel + "; npulses=" + npulses);
                }
                Long id = daqToGeometryId(crate, slot, channel);
                for (int k = 0; k < npulses; k++) {
                    short pulseNum = cdata.getByte();
                    int sampleCount = cdata.getNValue();

                    if (id == null) {
                        FADCGenericHit hit = makeGenericRawHit(EventConstants.HODO_PULSE_MODE, crate, slot, channel, cdata, sampleCount);
                        processUnrecognizedChannel(hit);
                    } else {
                        BaseRawTrackerHit hit = makeHodoRawHit(pulseNum, id, cdata, sampleCount);
                        hits.add(hit);
                    }
                }
            }
        }
        return hits;
    }

    private List<RawCalorimeterHit> makeIntegralHitsMode3(CompositeData cdata, int crate) {
        List<RawCalorimeterHit> hits = new ArrayList<RawCalorimeterHit>();
        if (debug) {
            int n = cdata.getNValues().size();
            for (int i = 0; i < n; i++) {
                System.out.println("cdata.N[" + i + "]=" + cdata.getNValues().get(i));
            }
            int ni = cdata.getItems().size();
            for (int i = 0; i < ni; i++) {
                System.out.println("cdata.type[" + i + "]=" + cdata.getTypes().get(i));
            }
        }
        while (cdata.index() + 1 < cdata.getItems().size()) { //the +1 is a hack because sometimes an extra byte gets read (padding)
            short slot = cdata.getByte();
            int trigger = cdata.getInt();
            long timestamp = cdata.getLong();
            int nchannels = cdata.getNValue();
            if (debug) {
                System.out.println("slot#=" + slot + "; trigger=" + trigger + "; timestamp=" + timestamp + "; nchannels=" + nchannels);
            }
            for (int j = 0; j < nchannels; j++) {
                short channel = cdata.getByte();
                int npulses = cdata.getNValue();
                if (debug) {
                    System.out.println("  channel=" + channel + "; npulses=" + npulses);
                }
                Long id = daqToGeometryId(crate, slot, channel);
                for (int k = 0; k < npulses; k++) {
                    short pulseTime = cdata.getShort();
                    int pulseIntegral = cdata.getInt();
                    if (debug) {
                        System.out.println("    pulseTime=" + pulseTime + "; pulseIntegral=" + pulseIntegral);
                    }
                    if (id == null) {
                        int[] data = {pulseIntegral, pulseTime};
                        processUnrecognizedChannel(new FADCGenericHit(EventConstants.HODO_PULSE_INTEGRAL_MODE, crate, slot, channel, data));
                    } else {
                        hits.add(new BaseRawCalorimeterHit(id, pulseIntegral, pulseTime));
                    }
                }
            }
        }
        return hits;
    }

    private List<RawCalorimeterHit> makeIntegralHitsMode7(CompositeData cdata, int crate) {
        List<RawCalorimeterHit> hits = new ArrayList<RawCalorimeterHit>();
        if (debug) {
            int n = cdata.getNValues().size();
            for (int i = 0; i < n; i++) {
                System.out.println("cdata.N[" + i + "]=" + cdata.getNValues().get(i));
            }
            int ni = cdata.getItems().size();
            for (int i = 0; i < ni; i++) {
                System.out.println("cdata.type[" + i + "]=" + cdata.getTypes().get(i));
            }
        }
        while (cdata.index() + 1 < cdata.getItems().size()) { //the +1 is a hack because sometimes an extra byte gets read (padding)
            short slot = cdata.getByte();
            int trigger = cdata.getInt();
            long timestamp = cdata.getLong();
            int nchannels = cdata.getNValue();
            if (debug) {
                System.out.println("slot#=" + slot + "; trigger=" + trigger + "; timestamp=" + timestamp + "; nchannels=" + nchannels);
            }
            for (int j = 0; j < nchannels; j++) {
                short channel = cdata.getByte();
                int npulses = cdata.getNValue();
                if (debug) {
                    System.out.println("  channel=" + channel + "; npulses=" + npulses);
                }
                Long id = daqToGeometryId(crate, slot, channel);
                for (int k = 0; k < npulses; k++) {
                    short pulseTime = cdata.getShort();
                    int pulseIntegral = cdata.getInt();
                    short amplLow = cdata.getShort();
                    short amplHigh = cdata.getShort();
                    if (debug) {
                        System.out.println("    pulseTime=" + pulseTime + "; pulseIntegral=" + pulseIntegral + "; amplLow=" + amplLow + "; amplHigh=" + amplHigh);
                    }
                    if (id == null) {
                        int[] data = {pulseIntegral, pulseTime, amplLow, amplHigh};
                        processUnrecognizedChannel(new FADCGenericHit(EventConstants.HODO_PULSE_INTEGRAL_HIGHRESTDC_MODE, crate, slot, channel, data));
                    } else {
                        RawCalorimeterHit hit = new BaseRawCalorimeterHit(id, pulseIntegral, pulseTime);
                        hits.add(hit);
                        Mode7Data extraData = new Mode7Data(amplLow, amplHigh);
                        extraDataList.add(extraData);
                        extraDataRelations.add(new BaseLCRelation(hit, extraData));
                    }
                }
            }
        }
        return hits;
    }

    private void processUnrecognizedChannel(FADCGenericHit hit) {
        genericHits.add(hit);

        List<Integer> channelAddress = Arrays.asList(hit.getCrate(), hit.getSlot(), hit.getChannel());
        Integer count = genericHitCount.get(channelAddress);
        if (count == null) {
            count = 0;
        }
        count++;
        genericHitCount.put(channelAddress, count);

        // Lowered the log level on these.  Otherwise they print too much. --JM
        if (count < 10) {
            LOGGER.finer(String.format("Crate %d, slot %d, channel %d not found in map", hit.getCrate(), hit.getSlot(), hit.getChannel()));
        } else if (count == 10) {
            LOGGER.fine(String.format("Crate %d, slot %d, channel %d not found in map: silencing further warnings for this channel", hit.getCrate(), hit.getSlot(), hit.getChannel()));
        }
    }

    void initialize() {

        DatabaseConditionsManager mgr = DatabaseConditionsManager.getInstance();
        Detector det = mgr.getDetectorObject();
        
        /**
         * The HodoEvioReader is always active, but we might be processing older data
         * with no Hodoscope subdetector in the compact detector description or no
         * Hodoscope data in the conditions database for the run.  Instead of crashing here, 
         * just check if this is the case and print a warning if either occurs.  No Hodoscope 
         * data should actually be processed by the job, or if there is some configuration
         * issue, crashing later with a null pointer exception during processing is acceptable.
         * Maybe a better way to do this would be configuring the builder only to use
         * the Hodoscope data reader if there is a Hodoscope detector in the compact
         * description. --JM
         */
        
        if (det.getSubdetectorNames().contains(subdetectorName)) {
            subDetector = DatabaseConditionsManager.getInstance().getDetectorObject().getSubdetector(subdetectorName);
            helper = subDetector.getDetectorElement().getIdentifierHelper();
        } else {
            LOGGER.warning("No Hodoscope subdetector was found in the detector description!");
        }

        if (mgr.hasConditionsRecord("hodo_channels")) {
            hodoChannels = mgr.getCachedConditions(HodoscopeChannelCollection.class, "hodo_channels").getCachedData();
        } else {
            LOGGER.warning("No HodoscopeChannelCollection found for run " + mgr.getRun() + " in the conditions database!");
        }
    }
}
