/*
 * DisplayResults.cpp
 *
 *  Created on: Jul 24, 2013
 *      Author: roubinet
 */


#include "DisplayResults.h"
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

using namespace std;

DisplayResults::DisplayResults() {}
DisplayResults::~DisplayResults() {}

DFNVisu* current_results;
DFNVisuSeg* current_results_seg;

void glutReshapeFunc(int width, int height){
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
}

void display_DFN(){current_results->display_DFN();}
void display_DFN_Seg(){current_results_seg->display_DFN_Seg();}

void SetupRC(){
	glClearColor(1,1,1,1);	// initialize window color to white
}

void keyboardFunc(unsigned char key, int x, int y){
	//current_results->keyboardFunc(key);
	glutPostRedisplay();
}

DisplayResults::DisplayResults(int argc,char **argv,DFNVisu results_){
	// class members affectation
	results = results_;
	current_results = &results_;

	int width=500;
	int height=500;

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(width,height);
	glutInitWindowPosition(200,200);
	glutInit(&argc, argv);
	glutCreateWindow("Discrete Fracture Network");
	glutKeyboardFunc(keyboardFunc);
	glutReshapeFunc(width,height);
	glutDisplayFunc(display_DFN);
	glutMainLoop();
}

DisplayResults::DisplayResults(int argc,char **argv,DFNVisuSeg results_seg_){
	// class members affectation
	results_seg=results_seg_;
	current_results_seg = &results_seg_;

	int width=500;
	int height=500;

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(width,height);
	glutInitWindowPosition(200,200);
	glutInit(&argc, argv);
	glutCreateWindow("Discrete Fracture Network");
	glutKeyboardFunc(keyboardFunc);
	glutReshapeFunc(width,height);
	glutDisplayFunc(display_DFN_Seg);
	glutMainLoop();
}
