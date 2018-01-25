#include <iostream>
#include <vector>

#include <omp.h>

#define STB_IMAGE_IMPLEMENTATION

#include "PhysicsEngine.h"
#include "GraphicsEngine.h"
#include "Definitions.h"

#include "ConstitutiveRelation.h"

int main (int argc, char ** argv)
{
	// physics engine initialization ------------------------------------------
	PhysicsEngine thePhysicsEngine;
	thePhysicsEngine.initializeWorld_Random();
	// graphics engine initialization -----------------------------------------
	GraphicsEngine theGraphicsEngine;
	// run simulation
	theGraphicsEngine.runVisualization(&thePhysicsEngine);

	return(0);
}

