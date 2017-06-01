package org.hps.recon.ecal;

public class TimeDependentGainFunc {
    
    private static final int iSTART = 0;
    private static final int iEND = 1;
    private static final int iA = 2;
    private static final int iB = 3;
    private static final int iC = 4; 
    public TimeDependentGainFunc(double[] params){
        this.values = params;
    }
    private double[] values;
    
    public double getFittedFunction(double timestamp){
        //System.out.println(timestamp);
        for(int i = 0; i<= values.length/5; i++){
            /*if(i == values.length/5){
                
                String message = "timestamp out of range:" + timestamp;
                message += "\nacceptable ranges are";
                for(int j = 0; j<values.length/5; j++){
                    message+= "\n" + values[j*5+iSTART] + " - " + values[j*5+iEND];
                }
                throw new RuntimeException(message);
            }*/
            double start = values[i*5+iSTART];
            double end = values[i*5+iEND];
            if(timestamp < end && timestamp > start){
                double A = values[i*5+iA];
                double B = values[i*5+iB];
                double C = values[i*5+iC];
                return A-B*Math.exp(-(timestamp-start)/C);
            }
        }
        return 1;
    }

}
