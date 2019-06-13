package org.hps.readout.triggerstudies;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import hep.aida.IHistogram;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.IHistogram3D;

public class PlotToTextModule {
    /**
     * Converts the argument histogram to a text file formatted for
     * import into Wolfram Mathematica as a list of points. Supported
     * object types are {@link hep.aida.IHistogram1D IHistogram1D},
     * {@link hep.aida.IHistogram2D IHistogram2D}, and {@link
     * hep.aida.IHistogram3D IHistogram3D}.
     * @param plot - The plot to convert.
     * @return Returns a {@link java.util.String} object that can be
     * read into Mathematica as a list of points.
     * @throws IllegalArgumentException Occurs if the argument object
     * type is not supported.
     */
    public static final String aidaToMathematica(IHistogram plot) throws IllegalArgumentException {
        if(plot instanceof IHistogram1D) {
            return aidaToMathematica((IHistogram1D) plot);
        } else if(plot instanceof IHistogram2D) {
            return aidaToMathematica((IHistogram2D) plot);
        } else if(plot instanceof IHistogram3D) {
            return aidaToMathematica((IHistogram3D) plot);
        } else {
            throw new IllegalArgumentException("Unrecognized histogram type \"" + plot.getClass().getSimpleName() + "\".");
        }
    }
    
    /**
     * Converts an AIDA histogram to a text file formatted for import
     * into Wolfram Mathematica as a list of points. Note that this
     * method does not handle writing the output file.
     * @param plot - The plot to convert.
     * @return Returns a {@link java.util.String} object that can be
     * read into Mathematica as a list of points.
     */
    public static final String aidaToMathematica(IHistogram1D plot) {
        StringBuffer outputBuffer = new StringBuffer();
        
        int bins = plot.axis().bins();
        
        for(int x = 0; x < bins; x++) {
            outputBuffer.append(String.format("%f,%f%n",plot.axis().binCenter(x), plot.binHeight(x)));
        }
        
        return outputBuffer.toString();
    }
    
    /**
     * Converts an AIDA histogram to a text file formatted for import
     * into Wolfram Mathematica as a list of points. Note that this
     * method does not handle writing the output file.
     * @param plot - The plot to convert.
     * @return Returns a {@link java.util.String} object that can be
     * read into Mathematica as a list of points.
     */
    public static final String aidaToMathematica(IHistogram2D plot) {
        StringBuffer outputBuffer = new StringBuffer();
        
        int xBins = plot.xAxis().bins();
        int yBins = plot.yAxis().bins();
        
        for(int x = 0; x < xBins; x++) {
            for(int y = 0; y < yBins; y++) {
                outputBuffer.append(String.format("%f,%f,%f%n", plot.xAxis().binCenter(x), plot.yAxis().binCenter(y), plot.binHeight(x, y)));
            }
        }
        
        return outputBuffer.toString();
    }
    
    /**
     * Converts an AIDA histogram to a text file formatted for import
     * into Wolfram Mathematica as a list of points. Note that this
     * method does not handle writing the output file.
     * @param plot - The plot to convert.
     * @return Returns a {@link java.util.String} object that can be
     * read into Mathematica as a list of points.
     */
    public static final String aidaToMathematica(IHistogram3D plot) {
        StringBuffer outputBuffer = new StringBuffer();
        
        int xBins = plot.xAxis().bins();
        int yBins = plot.yAxis().bins();
        int zBins = plot.zAxis().bins();
        
        for(int x = 0; x < xBins; x++) {
            for(int y = 0; y < yBins; y++) {
                for(int z = 0; z < zBins; y++) {
                    outputBuffer.append(String.format("%f,%f,%f,%f%n",
                            plot.xAxis().binCenter(x), plot.yAxis().binCenter(y), plot.zAxis().binCenter(z), plot.binHeight(x, y, z)));
                }
            }
        }
        
        return outputBuffer.toString();
    }
    
    /**
     * Converts an array of coordinates to a text file formatted for
     * import into Wolfram Mathematica as a list of points and then
     * writes it to the designated filepath.
     * @param coors - The coordinates to convert.
     * @param outputFile - The file to which the output should be
     * written.
     * @throws IOException Occurs if there is an issue with writing
     * the file.
     */
    public static final void writeCoordinateSet(Coordinate[] coors, File outputFile) throws IOException {
        writeFile(coordinateSetToMathematica(coors), outputFile);
    }
    
    /**
     * Converts an AIDA histogram to a text file formatted for import
     * into Wolfram Mathematica as a list of points and then writes
     * it to the designated filepath.
     * @param plot - The plot to convert.
     * @param outputFile - The file to which the output should be
     * written.
     * @throws IOException Occurs if there is an issue with writing
     * the file.
     */
    public static final void writePlot(IHistogram plot, File outputFile) throws IOException {
        writeFile(aidaToMathematica(plot), outputFile);
    }
    
    /**
     * Converts an array of coordinates to text for outputting to a
     * file.
     * @param coors - The coordinates.
     * @return Returns the list of coordinates in a text form that is
     * compatible with Mathematica.
     */
    private static final String coordinateSetToMathematica(Coordinate[] coors) {
        StringBuffer outputBuffer = new StringBuffer();
        for(Coordinate coor : coors) {
            outputBuffer.append(String.format("%f,%f%n", coor.getX(), coor.getY()));
        }
        return outputBuffer.toString();
    }
    
    /**
     * Writes the output text to the specified file.
     * @param outputText - The text to output.
     * @param outputFile - The file to output it to.
     * @throws IOException Occurs if there is an issue with writing
     * the file.
     */
    private static final void writeFile(String outputText, File outputFile) throws IOException {
        FileWriter writer = new FileWriter(outputFile);
        writer.write(outputText);
        writer.close();
    }
}