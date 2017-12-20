package org.hps.analysis.geometry;

import java.util.ArrayList;
import java.util.List;
import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.conditions.svt.SvtAlignmentConstant;
import org.lcsim.geometry.compact.converter.MilleParameter;

/**
 *
 * @author Norman A. Graf
 */
public class QuerySvtDatabase {

    public static void main(String[] args) {
        String detectorName = "HPS-EngRun2015-Nominal-v6-0-fieldmap";
        int runNum = 5772;
        if (args.length > 1) {
            detectorName = args[0];
        }
        if (args.length > 2) {
            runNum = Integer.getInteger(args[1]);
        }
        try {
            testSvtAlignment(detectorName, runNum);
        } catch (Exception e) {
            e.printStackTrace();

        }
    }

    public static void testSvtAlignment(String detectorName, int runNum) throws Exception {

        final DatabaseConditionsManager manager = new DatabaseConditionsManager();
        manager.setDetector(detectorName, runNum);

        final List<MilleParameter> milleParameters = new ArrayList<MilleParameter>();

        final SvtAlignmentConstant.SvtAlignmentConstantCollection alignmentConstants = manager.getCachedConditions(
                SvtAlignmentConstant.SvtAlignmentConstantCollection.class, "svt_alignments").getCachedData();

        for (final SvtAlignmentConstant constant : alignmentConstants) {
            final MilleParameter p = new MilleParameter(constant.getParameter(), constant.getValue(), 0.0);
            if (constant.getValue() != 0.0) {
                System.out.println(p);
            }
            milleParameters.add(p);
        }
    }

}
