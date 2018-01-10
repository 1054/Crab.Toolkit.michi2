#include "EPSPlot.h"
#include <math.h>

void drawClipping(EPSPlot &plot, float x, float y, float radius) {
	plot.setLineCap(0);
	plot.setLineJoin(1);
	const int points = 1025;
	const int rounds = 17;
	plot.setLineWidth(0.2f);
	float *linex = new float[points];
	float *liney = new float[points];
	for (int round = 0; (round < rounds); round++) {
		float tempradius = (float)(radius*sqrt(2)*round/rounds);
		for (int point = 0; (point < points); point++) {
			float angle = (float)(3.1415926535897932384*2*point/(points-1));
			float wobble = (float)((sin(angle*7-5)+1)*radius/2);
			linex[point] = x+(float)sin(angle)*(tempradius + wobble);
			liney[point] = y+(float)cos(angle)*(tempradius + wobble);
		}
		float shade = 1-(float)round/rounds;
		plot.setRGBColor((float)pow(shade, 1), (float)pow(shade, 2), (float)pow(shade, 4));
		plot.drawLines(linex, liney, points, 0.02f, x-radius, y-radius, x+radius, y+radius);
	}
	plot.drawText(x, y, "Clipping", 0.5f, -0.2f, 15);
	plot.drawText(x, y, "40x40", 0.5f, 1.2f, 15);
	plot.setLineJoin(0);
	plot.setGray(0);
	plot.drawBox(x-radius, y-radius, x+radius, y+radius);
	delete[] linex;
	delete[] liney;
}

void drawLines(EPSPlot &plot, float x, float y) {	
	plot.setGray(0);
	plot.setLineCap(1);
	plot.setFontSize(12);
	plot.drawText(x, y, "Lines", 1, 1);
	plot.setFontSize(6);
	plot.setLineWidth(0);
	plot.drawLegend(x, y-=15, "0pt");
	plot.setLineWidth(0.05f);
	plot.drawLegend(x, y-=6, "0.05pt");
	plot.setLineWidth(0.1f);
	plot.drawLegend(x, y-=6, "0.1pt");
	plot.setLineWidth(0.2f);
	plot.drawLegend(x, y-=6, "0.2pt");
	plot.setLineWidth(0.25f);
	plot.drawLegend(x, y-=6, "0.25pt");
	plot.setLineWidth(0.5f);
	plot.drawLegend(x, y-=6, "0.5pt");
	plot.setLineWidth(0.75f);
	plot.drawLegend(x, y-=6, "0.75pt");
	plot.setLineWidth(1.0f);
	plot.drawLegend(x, y-=6, "1.0pt");
	plot.setLineWidth(2.0f);
	plot.setLineCap(1);
	plot.drawLegend(x, y-=6, "(2pt) Round cap");
	plot.setLineCap(0);
	plot.drawLegend(x, y-=6, "Butt cap");
	plot.setLineCap(2);
	plot.drawLegend(x, y-=6, "Extended butt cap");
	plot.setLineCap(0);
	plot.setDash(1,2);
	plot.drawLegend(x, y-=6, "Dash [1, 2]");
	plot.setDash(2,1);
	plot.drawLegend(x, y-=6, "Dash [2, 1]");
	plot.setDash(4,1,1,1);
	plot.drawLegend(x, y-=6, "Dash [4, 1, 1, 1]");
}

void drawColors(EPSPlot &plot, float x, float y) {	
	plot.setGray(0);
	plot.setLineCap(0);
	plot.setNoDash();
	plot.setFontSize(12);
	plot.drawText(x, y, "RGB", 1, 1);
	plot.setFontSize(6);
	y -= 15;
	plot.setLineWidth(10);
	plot.setGray(0);
	plot.drawLine(x-2, y+3, x-2, y-39);
	plot.setLineWidth(2);
	plot.setRGBColor(1,0,0);
	plot.drawLegend(x, y, "(1,0,0)");
	plot.setRGBColor(0,1,0);
	plot.drawLegend(x, y-=6, "(0,1,0)");
	plot.setRGBColor(0,0,1);
	plot.drawLegend(x, y-=6, "(0,0,1)");
	plot.setRGBColor(0,1,1);
	plot.drawLegend(x, y-=6, "(0,1,1)");
	plot.setRGBColor(1,0,1);
	plot.drawLegend(x, y-=6, "(1,0,1)");
	plot.setRGBColor(1,1,0);
	plot.drawLegend(x, y-=6, "(1,1,0)");
	plot.setRGBColor(1,0.5f,0);
	plot.drawLegend(x, y-=6, "(1,0.5,0)");
}

void drawGrays(EPSPlot &plot, float x, float y) {	
	plot.setGray(0);
	plot.setLineCap(0);
	plot.setNoDash();
	plot.setFontSize(12);
	plot.drawText(x, y, "Gray", 1, 1);
	plot.setFontSize(6);
	plot.setGray(0);
	y -= 15;
	plot.setLineWidth(10);
	plot.drawLine(x-2, y+3, x-2, y-27);
	plot.setLineWidth(2);
	plot.drawLegend(x, y, "(0)");
	plot.setGray(0.25f);
	plot.drawLegend(x, y-=6, "(0.25)");
	plot.setGray(0.5f);
	plot.drawLegend(x, y-=6, "(0.5)");
	plot.setGray(0.75f);
	plot.drawLegend(x, y-=6, "(0.75)");
	plot.setGray(1);
	plot.drawLegend(x, y-=6, "(1)");
}

void drawJoint(EPSPlot &plot, float x, float y) {
	float linex[3] = {x, x-5, x-10};
	float liney[3] = {y-3, y+3, y-3};
	plot.drawLines(linex, liney, 3);
}

void drawJoins(EPSPlot &plot, float x, float y) {	
	plot.setGray(0);
	plot.setLineCap(0);
	plot.setLineWidth(2.5);
	plot.setNoDash();
	plot.setFontSize(12);
	plot.drawText(x, y, "Joins", 1, 1);
	plot.setFontSize(6);
	plot.setLineJoin(0);
	drawJoint(plot, x, y-=15);
	plot.drawText(x-12, y, "Miter join", 1, 0.5);
	plot.setLineJoin(1);
	drawJoint(plot, x, y-=8);
	plot.drawText(x-12, y, "Round join", 1, 0.5);
	plot.setLineJoin(2);
	drawJoint(plot, x, y-=8);
	plot.drawText(x-12, y, "Bevel join", 1, 0.5);
}

void drawFonts(EPSPlot &plot, float x, float y) {	
	plot.setGray(0);
	plot.setFontSize(12);
	plot.drawText(x, y, "Sizes", 1, 1);
	plot.setFontSize(4);
	plot.drawText(x, y-=15, "4pt font", 1);
	plot.setFontSize(6);
	plot.drawText(x, y-=6, "6pt font", 1);
	plot.setFontSize(8);
	plot.drawText(x, y-=8, "8pt font", 1);
	plot.setFontSize(10);
	plot.drawText(x, y-=10, "10pt font", 1);
	plot.setFontSize(12);
	plot.drawText(x, y-=12, "12pt font", 1);
}

void main() {
	EPSPlot plot("testpage.eps", 35, 0, 202, 101);
	drawClipping(plot, 58, 22, 20);
	drawLines(plot, 150, 100);
	drawJoins(plot, 200, 100);
	drawColors(plot, 110, 100);
	drawGrays(plot, 70, 100);
	drawFonts(plot, 200, 57);
}