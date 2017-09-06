package org.hps.recon.tracking;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.lcsim.fit.helicaltrack.HelicalTrackFit;
import org.lcsim.fit.helicaltrack.HelicalTrackHit;
import org.lcsim.fit.helicaltrack.HelicalTrackFitter.FitStatus;
import org.lcsim.recon.tracking.seedtracker.FastCheck;
import org.lcsim.recon.tracking.seedtracker.HelixFitter;
import org.lcsim.recon.tracking.seedtracker.HitManager;
import org.lcsim.recon.tracking.seedtracker.Sector;
import org.lcsim.recon.tracking.seedtracker.SeedCandidate;
import org.lcsim.recon.tracking.seedtracker.SeedLayer;
import org.lcsim.recon.tracking.seedtracker.SeedStrategy;
import org.lcsim.recon.tracking.seedtracker.SortHits;
import org.lcsim.recon.tracking.seedtracker.SortLayers;

public class ConfirmerExtender extends org.lcsim.recon.tracking.seedtracker.ConfirmerExtender {

    public ConfirmerExtender(HitManager hitmanager, HelixFitter helixfitter) {
        super(hitmanager, helixfitter);
    }

    public boolean Confirm(SeedCandidate seed, SeedStrategy strategy, double bfield, FastCheck checker, List<HelicalTrackHit> hitcol) {
        //  Create a list to hold the confirmed / extended seeds
        _result = new ArrayList<SeedCandidate>();

        //  Establish which layers are to be checked
        seed.setUncheckedLayers(strategy.getLayers(SeedLayer.SeedType.Confirm));

        //  Process the seed
        doTask(seed, Task.CONFIRM, strategy, bfield, checker, hitcol);

        //  Return true if we found at least one confirming seed candidate
        return _result.size() > 0;
    }

    /**
     * Try to extend a seed using a specified strategy.  The strategy specifies
     * the layers to use in extending the seed as well as the minimum number of
     * layers required for a seed to be a track candidate.  Any track candidates
     * found as a result of the extend operation are merged with the list of
     * track candidates found so far to eliminate inferior fits when a pair of
     * track candidates shares more than one hit.
     *
     * @param seed seed to be extended
     * @param strategy strategy to use
     * @param bfield magnetic field
     * @param foundseeds list of track candidates found so far
     */

    public void Extend(SeedCandidate seed, SeedStrategy strategy, double bfield, List<SeedCandidate> foundseeds, FastCheck checker, List<HelicalTrackHit> hitcol) {

        //  Initialize the list of extended seed candidates to those already found
        _result = foundseeds;

        //  Establish which layers are to be checked
        seed.setUncheckedLayers(strategy.getLayers(SeedLayer.SeedType.Extend));

        //  Extend the seed and return
        doTask(seed, Task.EXTEND, strategy, bfield, checker, hitcol);
        return;
    }

    /**
     * Perform the confirm or extend step.
     *
     * @param inputseed seed to be confirmed/extended
     * @param task confirm or extend (enum)
     * @param strategy strategy to use
     * @param bfield magnetic field
     */
    private void doTask(SeedCandidate inputseed, Task task, SeedStrategy strategy, double bfield, FastCheck checker, List<HelicalTrackHit> hitcol) {
        //  Initialize the counter for the number of fits performed on this seed
        _nfit = 0;

        if (this._applySectorBinning)
            checker.setDoSectorBinCheck(this._hmanager.getSectorManager());

        //  Calculate the minimum number of hits to succeed, retrieve the chisq cuts
        int minhits = strategy.getMinHits();
        if (task == Task.CONFIRM)
            minhits = strategy.getMinConfirm() + 3;
        double badhitchisq = strategy.getBadHitChisq();
        double maxchisq = strategy.getMaxChisq();

        //  Create a LIFO queue of seeds to be searched for a confirmation/extension
        // hit (note that a LIFO queue is used to minimize memory usage)
        LinkedList<SeedCandidate> seedlist = new LinkedList<SeedCandidate>();

        //  The bestseed is a SeedCandidate the meets the requirements for becoming
        //  a track, shares at least one hit with the inputseed, and has been deemed
        //  the best such seed by the track merging criteria
        //
        //  Initialize the best seed to null
        SeedCandidate bestseed = null;

        //  If we have already found track candidates, check for duplicates
        //  that share hits with the seed, finding the best such duplicate candidate.
        for (SeedCandidate trkcand : _result) {
            if (_merger.isDuplicate(inputseed, trkcand))
                bestseed = findBestCandidate(trkcand, bestseed, strategy);
        }

        //  Create a map between the SeedLayers to be checked and a list of hits on the layer to check
        Map<SeedLayer, List<HelicalTrackHit>> hitmap = new HashMap<SeedLayer, List<HelicalTrackHit>>();

        //  Loop over the layers to be checked
        for (SeedLayer lyr : inputseed.getUncheckedLayers()) {

            //  Create a list of hits to check on this layer
            List<HelicalTrackHit> hitlist = new ArrayList<HelicalTrackHit>();

            if (this._applySectorBinning) {
                //  Loop over the sectors on this layer to collect hits to check
                for (Sector sector : _hmanager.getSectors(lyr)) {

                    //  If there are no hits, skip this sector
                    if (sector.Hits().isEmpty())
                        continue;

                    //  See if this sector is consistent with this seed
                    if (!checker.CheckSector(inputseed, sector))
                        continue;

                    //  Add the hits for this sector
                    hitlist.addAll(sector.Hits());
                }
            } else {

                int layer = lyr.getLayer();
                for (HelicalTrackHit hth : hitcol) {
                    if (inputseed.getHits().get(0).z() * hth.z() <= 0)
                        continue;
                    if (hth.Layer() == layer)
                        hitlist.add(hth);
                }
            }

            //  Save the list of hits in the hitmap
            if (!hitlist.isEmpty())
                hitmap.put(lyr, hitlist);
        }

        //  Create a list of layers that have hits to check
        List<SeedLayer> lyrlist = new ArrayList<SeedLayer>();
        lyrlist.addAll(hitmap.keySet());

        //  Sort the layers in order of increasing number of hits
        SortLayers lyrsort = new SortLayers(hitmap);
        Collections.sort(lyrlist, lyrsort);

        //  Store the layers to check in the seed
        inputseed.setUncheckedLayers(lyrlist);

        //  Start with the input seed
        seedlist.add(inputseed);

        //  Keep looping until we have fully processed all seed candidates
        while (seedlist.size() > 0) {

            //  If we have exceeded the maximum number of fits, print warning and stop processing seed candidates
            if (_nfit > _maxfit) {
                System.out.println("Maximum number of fits exceeded in " + task.toString() + " step");
                if (bestseed == null) {
                    System.out.println("No track candidates are associated with the seed hits");
                } else {
                    System.out.println("Track candidate with " + bestseed.getHits().size() + " hits and chisq of " + bestseed.getHelix().chisqtot() + " associated with the seed hits");
                }
                break;
            }

            //  Pull the last seed off the queue (use a LIFO queue to minimize queue length)
            SeedCandidate seed = seedlist.poll();

            //  Check if there are enough unchecked layers to meet the minimum number of hits
            int lyrsleft = seed.getUncheckedLayers().size();
            int possiblehits = lyrsleft + seed.getHits().size();
            if (possiblehits < minhits)
                continue;

            //  If there is a best fit candidate, see if there is still a chance of beating it
            if (bestseed != null) {

                //  If the maximimum hits we can achieve is >1 fewer than the best fit, skip this candidate
                int besthits = bestseed.getHits().size();
                if (possiblehits < besthits - 1)
                    continue;

                //  If the maximum hits we can achieve equals the best fit, skip if we have a worse chi2
                double chisq = seed.getHelix().chisqtot();
                double bestchisq = seed.getHelix().chisqtot();
                if ((possiblehits == besthits) && chisq > bestchisq)
                    continue;

                //  If the maximum hits we can achieve is 1 fewer than the best fit, skip if the bad hit criteria can't be met
                if ((possiblehits == besthits - 1) && (chisq > bestchisq - badhitchisq))
                    continue;
            }

            //  See if there are any layers left for confirm/extend
            if (lyrsleft == 0) {

                //  Take final action on this seed
                if (task == Task.CONFIRM) {

                    //  No more layers and min hit requirement is met, seed is confirmed
                    _result.add(seed);

                } else if (task == Task.EXTEND) {

                    //  Merge the seed into the list of extended seeds
                    boolean merged = _merger.Merge(_result, seed, strategy);

                    //  If the seed survived the merge, make it our new best candidate
                    if (merged)
                        bestseed = findBestCandidate(seed, bestseed, strategy);
                }

                //  Done with this seed
                continue;
            }

            //  Pull the next layer off the queue
            SeedLayer lyr = seed.getNextLayer();
            HelicalTrackFit helix = seed.getHelix();

            //  Retrieve the chisq for the last fit and initialize the best fit chisq for this layer
            double oldchisq = helix.chisqtot();
            double oldcirclechisq = helix.chisq()[0];
            double chisqbest = 1.e99;

            //  Get the list of hits to check for this layer and sort them by x-y distance from current helix
            List<HelicalTrackHit> hitlist = hitmap.get(lyr);
            SortHits comp = new SortHits(helix);
            Collections.sort(hitlist, comp);

            //  Loop over the sorted hits in this layer
            for (HelicalTrackHit hit : hitlist) {

                //  Make a test seed including the new hit
                SeedCandidate test = new SeedCandidate(seed);
                test.addHit(hit);

                //  Check that this hit is potentially viable
                if (!checker.CheckHitSeed(hit, seed)) {
                    if (_diag != null)
                        _diag.fireCheckHitFailed(hit, test);
                    continue;
                }

                //  Fit the test seed
                boolean success = _fitter.FitCandidate(test, strategy);
                _nfit++;

                //  Check if the fit was successful
                if (success) {

                    //  Success - attach the fit to the test seed
                    HelicalTrackFit newhelix = _fitter.getHelix();
                    test.setHelix(newhelix);

                    //  Add the seed to the LIFO queue of seed candidates and update the best chisq
                    seedlist.addFirst(test);
                    chisqbest = Math.min(chisqbest, newhelix.chisqtot());

                } else {

                    //  Stop checking hits in this layer if circle chisq increase is too big
                    if (_fitter.getFitStatus() != FitStatus.CircleFitFailed) {
                        double circlechisq = _fitter.getCircleFit().chisq();
                        if (circlechisq > oldcirclechisq + maxchisq)
                            break;
                    }
                }
            }

            //  Finished checking hits in the current layer.  If all the fit trials for
            //  this layer are potentially bad hits, include the starting seed (less the
            //  current layer, which was popped off the layer queue) in the seed list.
            if (chisqbest - oldchisq > strategy.getBadHitChisq())
                seedlist.addFirst(seed);
        }

        //  Finished looping over the seeds in the LIFO candidate queue - we are done!
        return;
    }
}
