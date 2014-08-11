/*
	$Revision: 1.1 $	$Date: 2009/11/04 23:54:15 $

	Definitions for 3D geometrical structures
*/

typedef struct Point {
	double x;
	double y;
	double z;
} Point;

typedef struct Box {
	double x0;
	double x1;
	double y0;
	double y1;
	double z0;
	double z1;
} Box;
