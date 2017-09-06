package org.hps.recon.tracking;

import java.util.List;

import org.lcsim.fit.helicaltrack.HelicalTrackHit;
import org.lcsim.recon.tracking.seedtracker.FastCheck;
import org.lcsim.recon.tracking.seedtracker.HelixFitter;
import org.lcsim.recon.tracking.seedtracker.HitManager;
import org.lcsim.recon.tracking.seedtracker.SeedCandidate;
import org.lcsim.recon.tracking.seedtracker.SeedStrategy;

public class SeedTrackFinder extends org.lcsim.recon.tracking.seedtracker.SeedTrackFinder {

    public SeedTrackFinder(HitManager hitmanager, HelixFitter helixfitter) {
        super(hitmanager, helixfitter);

        _confirmer = new ConfirmerExtender(_hitmanager, _helixfitter);

    }

    public boolean FindTracks(SeedStrategy strategy, double bfield, FastCheck checker, List<HelicalTrackHit> hitcol) {
        //  Loop over the first seed layer
        for (HelicalTrackHit hit1 : hitcol) {

            //  Loop over the second seed layer and check that we have a hit pair consistent with our strategy
            for (HelicalTrackHit hit2 : hitcol) {

                //  Check if the pair of hits is consistent with the current strategy
                if (!checker.TwoPointCircleCheck(hit1, hit2, null)) {
                    continue;
                }
                if (hit1.z() * hit2.z() <= 0)
                    continue;

                //  Loop over the third seed layer and check that we have a hit triplet consistent with our strategy
                for (HelicalTrackHit hit3 : hitcol) {

                    if (hit1.z() * hit3.z() <= 0)
                        continue;

                    //  Form a seed candidate from the seed hits
                    SeedCandidate seed = new SeedCandidate(strategy, bfield);
                    seed.addHit(hit1);
                    seed.addHit(hit2);
                    seed.addHit(hit3);

                    //  Check if the triplet of hits is consistent with the current strategy
                    if (!checker.ThreePointHelixCheck(hit1, hit2, hit3)) {
                        continue;
                    }

                    //  Form a seed candidate from the seed hits
                    //  See if we can fit a helix to this seed candidate
                    boolean success = _helixfitter.FitCandidate(seed, strategy);
                    if (!success)
                        continue;

                    seed.setHelix(_helixfitter.getHelix());
                    success = ((ConfirmerExtender) _confirmer).Confirm(seed, strategy, bfield, checker, hitcol);
                    if (!success)
                        continue;

                    //  Try to extend each confirmed seed candidates to make a track candidate
                    List<SeedCandidate> confirmedlist = _confirmer.getResult();
                    for (SeedCandidate confirmedseed : confirmedlist) {

                        //  See if we can extend this seed candidate
                        ((ConfirmerExtender) _confirmer).Extend(confirmedseed, strategy, bfield, _trackseeds, checker, hitcol);
                    }
                }
            }
        }

        return _trackseeds.size() > 0;
    }
}
