package org.hps.readout.triggerstudies;

public class Coordinate implements Comparable<Coordinate> {
    private final double x;
    private final double y;
    
    public Coordinate(double x, double y) {
        this.x = x;
        this.y = y;
    }
    
    public double getX() { return x; }
    
    public double getY() { return y; }
    
    @Override
    public String toString() {
        return String.format("(%f, %f)", x, y);
    }
    
    @Override
    public boolean equals(Object o) {
        if(!(o instanceof Coordinate)) { return false; }
        else {
            Coordinate coor = (Coordinate) o;
            return (x == coor.x && y == coor.y);
        }
    }
    
    @Override
    public int compareTo(Coordinate coor) {
        if(x == coor.x) {
            return Double.compare(y, coor.y);
        } else {
            return Double.compare(x, coor.x);
        }
    }
}