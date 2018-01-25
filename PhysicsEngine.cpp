#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
PhysicsEngine::PhysicsEngine()
{
};
// ----------------------------------------------------------------------------
PhysicsEngine::~PhysicsEngine()
{
	//delete all objects created by Factory classes
	for(unsigned int index = 0; index < allParticle.size(); index++)
        delete allParticle[index];

	for(unsigned int index = 0; index < allWall.size(); index++)
        delete allWall[index];
}
// ----------------------------------------------------------------------------
void PhysicsEngine::reportConsole(std::string sDescription)
{
	std::string strConsole = "";
	strConsole += sDescription;
	strConsole += "\ttime: " + Script(d_Time,4);
	strConsole += "\tRuntime: " + Script(d_Runtime_Total,3);

	glm::dvec3 d3Momentum = glm::dvec3(0.0,0.0,0.0);
	for(unsigned int index = 0; index < allParticle.size(); index++)
	{// momentum
		Particle_CC *thisP = allParticle[index];

		d3Momentum += thisP->d_Mass * thisP->d3_Velocity;
	}
	strConsole += "\tMomentum(y): " + Script(d3Momentum.y, 4);

	strConsole += "\n";

	std::cout << strConsole;

//	std::string strFileName = str_LogFile;
//	std::ofstream OutputFile(strFileName.c_str(), std::ios_base::app);
//
//	OutputFile << strConsole;
//
//	OutputFile.close();
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
