/*
 * DisplayResults.h
 *
 *  Created on: Jul 24, 2013
 *      Author: roubinet
 */

#ifndef DISPLAYRESULTS_H_
#define DISPLAYRESULTS_H_


#include <iostream>
#include "DFNVisu.h"


class DisplayResults {
private:
	DFNVisu results;
	DFNVisuSeg results_seg;
public:
	DisplayResults();
	virtual ~DisplayResults();
	DisplayResults(int,char **,DFNVisu);
	DisplayResults(int,char **,DFNVisuSeg);
};


#endif /* DISPLAYRESULTS_H_ */
