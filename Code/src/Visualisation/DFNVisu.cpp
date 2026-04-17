/*
 * DFNVisu.cpp
 *
 *  Created on: Jul 24, 2013
 *      Author: roubinet
 */


#include "DFNVisu.h"
#include <GL/gl.h>
#include <GL/glut.h>
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

DFNVisu::DFNVisu(){}
DFNVisu::DFNVisu(NetworkMeshes net_mesh_){
	net_mesh=net_mesh_;
}

DFNVisuSeg::DFNVisuSeg(){}
DFNVisuSeg::DFNVisuSeg(MapSegmentTimesCum seg_results_,MapSegmentVelocities segment_velocities_,Domain domain_,double snapshot_time_){
	seg_results=seg_results_;
	segment_velocities=segment_velocities_;
	domain=domain_;
	snapshot_time=snapshot_time_;
}

void DFNVisu::display_DFN(){
	// 1. window setting
	glClear(GL_COLOR_BUFFER_BIT);	// clear the window with current clearing color
	glTranslatef(-1,-1,-1);	// change coordinate reference (translate it from the center to the bottom left of the domain)
	double Lx=net_mesh.domain.domain_size_x(),Ly=net_mesh.domain.domain_size_y();
	glScalef(2./Lx,2./Ly,1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	// 2. value display
	// draw a polygon for each cell with a color representing its value
	glPolygonMode(GL_FRONT,GL_FILL);	// polygon with a color
	// white polygon for the back
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	SetColor(-1.);
	glBegin(GL_POLYGON);
	double factor=0.5*Lx;
	glVertex2f(net_mesh.domain.min_pt.i+factor,net_mesh.domain.min_pt.j+factor);
	glVertex2f(net_mesh.domain.min_pt.i+factor,net_mesh.domain.max_pt.j+factor);
	glVertex2f(net_mesh.domain.max_pt.i+factor,net_mesh.domain.max_pt.j+factor);
	glVertex2f(net_mesh.domain.max_pt.i+factor,net_mesh.domain.min_pt.j+factor);
	glEnd();
	// plot each fracture
	vector<FractureMesh> meshes=net_mesh.meshes;
	double factor_ap=1;
	for (int i=0;i<meshes.size();i++){
		glLineWidth(meshes[i].aperture*factor_ap*100);
		glLineWidth(1.0);
		glBegin(GL_LINES);
		SetColor(0.);
		glVertex2f(meshes[i].p_ori.p.x()+factor,meshes[i].p_ori.p.y()+factor);
		glVertex2f(meshes[i].p_tar.p.x()+factor,meshes[i].p_tar.p.y()+factor);
		glEnd();
	}

	// polygon contour
	//double contour_width=1;
	double contour_width=10;
	SetColor(0.6);
	glLineWidth(contour_width);
	glBegin(GL_LINES);
	glVertex2f(0,0);
	glVertex2f(0,Ly);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(Lx,0);
	glVertex2f(Lx,Ly);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(0,0);
	glVertex2f(Lx,0);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(0,Ly);
	glVertex2f(Lx,Ly);
	glEnd();

	glFlush();
}

void DFNVisuSeg::display_DFN_Seg(){
	// 1. window setting
	glClear(GL_COLOR_BUFFER_BIT);	// clear the window with current clearing color
	glTranslatef(-1,-1,-1);	// change coordinate reference (translate it from the center to the bottom left of the domain)
	double Lx=domain.domain_size_x(),Ly=domain.domain_size_y();
	glScalef(2./Lx,2./Ly,1.);	// scale window as its spatial discretization (from [-1 1] to [0 N]) (Ny-1 to avoid ghost display)
	double factor=0.5*Lx;
	double factor_fracture_aperture=100;	// factor to increase the fracture aperture
	// 2. value display
	// draw a polygon for each cell with a color representing its value
	glPolygonMode(GL_FRONT,GL_FILL);	// polygon with a color
	// white polygon for the back
	glEnable(GL_LINE_SMOOTH);
			glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glBegin(GL_POLYGON);
	SetColor(-1.);
	glVertex2i(domain.min_pt.i+factor,domain.min_pt.j+factor);
	glVertex2i(domain.min_pt.i+factor,domain.max_pt.j+factor);
	glVertex2i(domain.max_pt.i+factor,domain.max_pt.j+factor);
	glVertex2i(domain.max_pt.i+factor,domain.min_pt.j+factor);
	glEnd();
	// plot each fracture
	double part_prop;map<double,double>::iterator it3;bool value_found;
	for (MapSegmentTimesCum::iterator it1=seg_results.begin();it1!=seg_results.end();it1++){
		part_prop=-1;value_found=false;
		for (map<double,double>::iterator it2=it1->second.begin();it2!=it1->second.end();it2++){
			if (it2!=it1->second.end()){
				it3=it2;it3++;
				if (it2->first<snapshot_time&&snapshot_time<it3->first){
					part_prop=it2->second;value_found=true;
				}
			}
		}
		if (!value_found){
			it3=it1->second.end();it3--;
			part_prop=it3->second;
		}
		glLineWidth(segment_velocities[it1->first]*factor_fracture_aperture);
		glBegin(GL_LINES);
		//SetColorBis(part_prop);
		SetColor(part_prop);
		cout << "part_prop=" << part_prop << endl;
		glVertex2f(it1->first.first.i+factor,it1->first.first.j+factor);
		glVertex2f(it1->first.second.i+factor,it1->first.second.j+factor);
		glEnd();
	}
	glFlush();
}

// set a color corresponding to x value for 0<x<1
void SetColor(double x){

	double r,g,b;

	//if((x>=0)&&(x<0.25)) {r=0; g=0; b=1;}
	if((x>=0)&&(x<0.05)) {r=0; g=0; b=1;}
	else if((x>0.05)&&(x<=0.25)) {r=0; g=1; b=0;}
	else if((x>0.25)&&(x<=0.5)) {r=0.33; g=0; b=0.66;}
	else if((x>0.5)&&(x<=0.75)) {r=0.66; g=0; b=0.33;}
	else if((x>0.75)&&(x<=1)) {r=1; g=0; b=0;}
	else if (x<0){r=1; g=1; b=1;}
	else{
		cout << "WARNING in SetColor (VisualizationUtilitary.cpp): color not defined" << endl;
		cout << "x = " << x << endl;
/* flowing through different levels of  purple between blue and red if we wanted to do something very smooth we would do it for each discretization according to the code below. Which we might not have time for
	if((x>=0)&&(x<0.125)) {r=0; g=0; b=1;}
	else if((x>0.125)&&(x<=0.375)) {r=0.25; g=0; b=0.75;}
	else if((x>0.375)&&(x<=0.625)) {r=0.5; g=0; b=0.5;}
	else if((x>0.625)&&(x<=0.875)) {r=0.75; g=0; b=0.25;}
	else if((x>0.875)&&(x<=1)) {r=1; g=0; b=0;}
	else if (x<0){r=1; g=1; b=1;}
	else{
		cout << "WARNING in SetColor (VisualizationUtilitary.cpp): color not defined" << endl;
		cout << "x = " << x << endl;

		*/
/*
	if((x>=0)&&(x<0.125)) {r=0; g=0; b=0.5+4*x;}
	else if((x>0.125)&&(x<=0.375)) {r=0; g=4*x-0.5; b=1;}
	else if((x>0.375)&&(x<=0.625)) {r=4*x-1.5; g=1; b=2.5-4*x;}
	else if((x>0.625)&&(x<=0.875)) {r=1; g=3.5-4*x; b=0;}
	else if((x>0.875)&&(x<=1)) {r=4.5-4*x; g=0; b=0;}
	else if (x<0){r=1; g=1; b=1;}
	else{
		cout << "WARNING in SetColor (VisualizationUtilitary.cpp): color not defined" << endl;
		cout << "x = " << x << endl;
		*/
	}

	glColor3f(r,g,b);


/*
	typedef unsigned char uint8;

	struct Colour
	{
	    uint8 red;
	    uint8 green;
	    uint8 blue;

	    Colour(uint8 r, uint8 g, uint8 b) : red(r), green(g), blue(b) {}
	};

	Colour Source(0, 0, 255/insert x of the previous segment);  // blue
	Colour Target(255/x, 0, 0);  // red

	Colour MyColour = Source;

	const int NumSteps = 100;

	for (int i = 0; i <= NumSteps; i++)
	{
	  MyColour.red   = Source.red   + (((Target.red   - Source.red)   * i) / NumSteps);
	  MyColour.green = Source.green + (((Target.green - Source.green) * i) / NumSteps);
	  MyColour.blue  = Source.blue  + (((Target.blue  - Source.blue)  * i) / NumSteps);

	  // Do something, like update the display
	 // DoSomethingWithColour(MyColour);

	  glColor3f(MyColour.red,MyColour.green,MyColour.blue);
	}
*/



}

void SetColorBis(double x){
	double r,g,b;
	if((x>=0)&&(x<0.33)) {r=0; g=0; b=1;}
	else if((x>=0.33)&&(x<=0.66)) {r=0; g=1; b=0;}
	else if((x>=0.66)&&(x<=1)) {r=1; g=0; b=0;}
	else{
		cout << "WARNING in SetColor (VisualizationUtilitary.cpp): color not defined" << endl;
		cout << "x = " << x << endl;
	}
	glColor3f(r,g,b);
}

// set a color corresponding to x value for x_min<x<x_max
void SetColor(double x, double x_min, double x_max){
	SetColor((x-x_min)/(x_max-x_min));
}
