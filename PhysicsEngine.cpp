#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
PhysicsEngine::PhysicsEngine()
{
	// create pool for PPC links
	v_Link_PPC.resize(_MAX_PARTICLE * _MAX_PPC_PER_PARTICLE);
	for(unsigned int index = 0; index < v_Link_PPC.size(); index++)
	{
		Link_PPC_CC *newLink = new Link_PPC_CC;

		v_Link_PPC[index] = newLink;
	}
	// create pool for PWC links
	v_Link_PWC.resize(_MAX_PARTICLE * _MAX_PWC_PER_PARTICLE);
	for(unsigned int index = 0; index < v_Link_PWC.size(); index++)
	{
		Link_PWC_CC *newLink = new Link_PWC_CC;

		v_Link_PWC[index] = newLink;
	}
};
// ----------------------------------------------------------------------------
PhysicsEngine::~PhysicsEngine()
{
	// delete all objects created by Factory classes

	// particles
	for(unsigned int index = 0; index < v_Particle.size(); index++)
        delete v_Particle[index];
	// walls
	for(unsigned int index = 0; index < v_Wall.size(); index++)
        delete v_Wall[index];
	// PPC links
	for(unsigned int index = 0; index < v_Link_PPC.size(); index++)
		delete v_Link_PPC[index];
	// PWC links
	for(unsigned int index = 0; index < v_Link_PWC.size(); index++)
		delete v_Link_PWC[index];
}
// ----------------------------------------------------------------------------
void PhysicsEngine::reportConsole(std::string sDescription)
{
	std::string strConsole = "";
	strConsole += sDescription;
	strConsole += "\ttime: " + Script(d_Time,4);
	strConsole += "\tRuntime: " + Script(d_Runtime_Total,3);

	glm::dvec3 d3Momentum = glm::dvec3(0.0,0.0,0.0);
	for(unsigned int index = 0; index < v_Particle.size(); index++)
	{// momentum
		Particle_CC *thisP = v_Particle[index];

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
