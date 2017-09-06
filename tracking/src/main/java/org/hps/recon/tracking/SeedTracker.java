/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.hps.recon.tracking;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.lcsim.event.EventHeader;
import org.lcsim.event.MCParticle;
import org.lcsim.fit.helicaltrack.HelicalTrackHit;
import org.lcsim.recon.tracking.seedtracker.SeedCandidate;
import org.lcsim.recon.tracking.seedtracker.SeedStrategy;

/**
 * Class extending lcsim version to allow extra flexibility
 * 
 * @author Per Hansson Adrian <phansson@slac.stanford.edu>
 */
public class SeedTracker extends org.lcsim.recon.tracking.seedtracker.SeedTracker {

    private int _iterativeConfirmedFits = 0;
    private boolean doIterativeHelix = false;
    private boolean debug;
    private boolean doSectoring = false;

    public SeedTracker(List<SeedStrategy> strategylist) {
        // use base class only if this constructor is called!
        super(strategylist);
    }

    public void setIterativeHelix(boolean value) {
        doIterativeHelix = value;
    }

    private void initialize(List<SeedStrategy> strategylist, boolean useHPSMaterialManager, boolean includeMS) {

        // Explicitly only replace the objects that might change to avoid getting the lcsim versions

        //  Instantiate the material manager for HPS,  the helix fitter and seed track finder as tey depends on the material manager
        if (useHPSMaterialManager) {
            MaterialSupervisor materialSupervisor = new MaterialSupervisor(includeMS);
            materialSupervisor.setDebug(true);
            _materialmanager = materialSupervisor;
            _helixfitter = new HelixFitter(materialSupervisor, doIterativeHelix);
        } else {
            MaterialManager materialmanager = new MaterialManager(includeMS);
            _materialmanager = materialmanager; //mess around with types here...
            _helixfitter = new HelixFitter(materialmanager, doIterativeHelix);
        }

        //  Instantiate the helix finder since it depends on the material manager
        _finder = new SeedTrackFinder(_hitmanager, _helixfitter);

    }

    public SeedTracker(List<SeedStrategy> strategylist, boolean includeMS) {
        super(strategylist);
        initialize(strategylist, true, includeMS);
    }

    public SeedTracker(List<SeedStrategy> strategylist, boolean useHPSMaterialManager, boolean includeMS) {
        super(strategylist);
        initialize(strategylist, useHPSMaterialManager, includeMS);
    }

    public SeedTracker(List<SeedStrategy> strategylist, boolean useHPSMaterialManager, boolean includeMS, boolean doIterative) {
        super(strategylist);
        setIterativeHelix(doIterative);
        initialize(strategylist, useHPSMaterialManager, includeMS);
    }

    public void setIterativeConfirmed(int maxfits) {
        this._iterativeConfirmedFits = maxfits;
        super.setIterativeConfirmed(maxfits);
    }

    /**
    * Set to enable debug output
    * 
    * @param debug switch
    */
    @Override
    public void setDebug(boolean debug) {
        super.setDebug(debug);
        _materialmanager.setDebug(debug);
        _helixfitter.setDebug(debug);
        this.debug = debug;
    }

    /**
    * Set to enable the sectoring to use the sector bins in checking for consistent hits.
    *
    * @param applySectorBinning apply sector binning switch
    */
    public void setApplySectorBinning(boolean applySectorBinning) {
        _finder.setApplySectorBinning(applySectorBinning);
        _finder.getConfirmer().setApplySectorBinning(applySectorBinning);
        doSectoring = false;
    }

    public void setSubdetectorName(String subdetectorName) {
        ((MaterialSupervisor) this._materialmanager).setSubdetectorName(subdetectorName);
    }

    @Override
    protected void process(EventHeader event) {
        // Pass the event to the diagnostics package
        if (_diag != null)
            _diag.setEvent(event);

        // Get the hit collection from the event
        List<HelicalTrackHit> hitcol = event.get(HelicalTrackHit.class, _inputCol);

        // Sort the hits for this event
        _hitmanager.OrganizeHits(hitcol);

        // Make sure that we have cleared the list of track seeds in the finder
        _finder.clearTrackSeedList();

        // Loop over strategies and perform track finding
        for (SeedStrategy strategy : _strategylist) {

            // Set the strategy for the diagnostics
            if (_diag != null)
                _diag.fireStrategyChanged(strategy);

            FastCheck checker = new FastCheck(strategy, _bfield, _diag);
            // Perform track finding under this strategy
            if (doSectoring)
                ((SeedTrackFinder) _finder).FindTracks(strategy, _bfield, checker);
            else {
                ((SeedTrackFinder) _finder).FindTracks(strategy, _bfield, checker, hitcol);
            }
        }

        // Get the list of final list of SeedCandidates
        List<SeedCandidate> trackseeds = _finder.getTrackSeeds();
        ((HelixFitter) _helixfitter).setIterative(true);

        if (_iterativeConfirmedFits > 0) {
            // Iteratively re-fit tracks to take into account helix and hit position correlations
            if (this.debug)
                System.out.printf("%s: Iteratively improve %d seeds\n", this.getClass().getSimpleName(), trackseeds.size());
            List<SeedCandidate> seedsToRemove = new ArrayList<SeedCandidate>();
            for (SeedCandidate seed : trackseeds) {
                SeedStrategy strategy = seed.getSeedStrategy();
                boolean success = false;
                for (int iterFit = 0; iterFit < _iterativeConfirmedFits; ++iterFit) {
                    success = _helixfitter.FitCandidate(seed, strategy);
                }
                if (!success) {
                    seedsToRemove.add(seed);
                } else {
                    if (this.debug)
                        System.out.printf("%s: done iterating, this seed will be added to event:\n%s\n", this.getClass().getSimpleName(), seed.toString());
                }
            }
            for (SeedCandidate badseed : seedsToRemove) {
                trackseeds.remove(badseed);
            }
        }

        // Make tracks from the final list of track seeds
        _maketracks.Process(event, trackseeds, _bfield);

        // Save the MC Particles that have been seeded / confirmed if diagnostics are enabled
        if (_diag != null) {
            Set<MCParticle> seededmcpset = _finder.getSeededMCParticles();
            List<MCParticle> seededmcp = new ArrayList<MCParticle>(seededmcpset);
            event.put("SeededMCParticles", seededmcp, MCParticle.class, 0);
            Set<MCParticle> confirmedmcpset = _finder.getConfirmedMCParticles();
            List<MCParticle> confirmedmcp = new ArrayList<MCParticle>(confirmedmcpset);
            event.put("ConfirmedMCParticles", confirmedmcp, MCParticle.class, 0);
        }

        // Clear the list of track seeds accumulated in the track finder
        _finder.clearTrackSeedList();

        return;
    }

}
