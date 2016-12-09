/*
 * Copyright © Alexandre Boulch, Renaud Marlet, Ecole Nationale des Ponts et Chaussees
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this
 * software and associated documentation files (the "Software"), to deal in the Software
 * without restriction, including without limitation the rights to use, copy, modify,
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following
 * conditions:

 * The above copyright notice and this permission notice shall be included in all copies
 * or substantial portions of the Software.

 * The Software is provided "as is", without warranty of any kind, express or implied,
 * including but not limited to the warranties of merchantability, fitness for a
 * particular purpose and noninfringement. In no event shall the authors or copyright
 * holders be liable for any claim, damages or other liability, whether in an action of
 * contract, tort or otherwise, arising from, out of or in connection with the software
 * or the use or other dealings in the Software.
 */


import java.util.ArrayList;
import java.util.Random;

import primitives_compare.Surface;
import primitives_compare.surface_comparator_KS;
import primitives_compare.surface_comparator_MW;

class Point{
	double x,y,z;

	public Point(){x = y = z = 0;}
	public Point(double _x, double _y, double _z){
		x = _x;
		y = _y;
		z = _z;
	}

	public double dot( Point pt){
		return x*pt.x+y*pt.y+z*pt.z;
	}

	public double norm(){
		return Math.sqrt(x*x+y*y+z*z);
	}

	public Point add(Point pt){
		return new Point(x+pt.x, y+pt.y, z+pt.z);
	}

	public Point minus(Point pt){
		return new Point(x-pt.x, y-pt.y, z-pt.z);
	}

	public Point mult(double v){
		return new Point(x*v, y*v, z*v);
	}
}


class Plane extends Surface<Point>{
	Point point;
	Point normal;

	public double distance(Point p){
		return Math.abs(normal.dot(p.minus(point)));
	}
}



public class Primitive_Fusion_example{


	public static void main(String[] args){
		System.out.println("Primitive Fusion Example");

		Random rand = new Random();
		
		ArrayList<Point> P1 = new ArrayList<Point>();
		ArrayList<Point> P2 = new ArrayList<Point>();

		double dist_plane =0.01;
		int nbr_pts = 1000;
		double noise = 0.03;
		
		for(int i=0; i<nbr_pts; i++){

			Point pt1 = new Point(rand.nextDouble(), rand.nextDouble(), 0);
			Point pt2 = new Point(rand.nextDouble(), rand.nextDouble(), dist_plane);
			
			//box Muller method for normal variable
			double U = rand.nextDouble() ;
			double V = rand.nextDouble();
			double dir1 = Math.sqrt(-2*Math.log(U))*Math.cos(2*3.14159265*V)*noise;
			U = rand.nextDouble() ;
			V = rand.nextDouble() ;
			double dir2 = Math.sqrt(-2*Math.log(U))*Math.cos(2*3.14159265*V)*noise;

			pt1.z += dir1;
			pt2.z += dir2;
			
			P1.add(pt1);
			P2.add(pt2);
			
		}
		
		
		Plane pl1 = new Plane();
		Plane pl2 = new Plane();
		pl1.point = new Point(0,0,0);
		pl2.point = new Point(0,0,dist_plane);
		pl1.normal = pl2.normal = new Point(0,0,1);
		

		System.out.println("SCMW");
		surface_comparator_MW<Point, Plane, Plane, Surface<Point> > scmw = 
				new surface_comparator_MW<Point, Plane, Plane, Surface<Point> >(pl1,P1,pl2,P2);
		System.out.print("MW: ");
		System.out.println(scmw.apply_MW_test(0.01, 0.05));

		System.out.println("SCKS");
		surface_comparator_KS<Point, Plane, Plane, Surface<Point> > scks = 
				new surface_comparator_KS<Point, Plane, Plane, Surface<Point> >(pl1,P1,pl2,P2);
		System.out.print("KS: ");
		System.out.println(scks.apply_KS_test(0.01, 0.05));


	}

}