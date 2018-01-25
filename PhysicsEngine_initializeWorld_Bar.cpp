#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_Random(void)
{
	// size of world
	d3_Length_World = glm::dvec3(0.100, 0.100, 0.100);

	if(true)
	{// random particles ------------------------------------------------------
		std::vector<Particle_CC *> thisParticleDomain;
		for(unsigned int index_P = 0; index_P < 10; index_P++)
		{// assign material point initial values
			Particle_CC *thisP = new Particle_CC;
			thisParticleDomain[index_P] = thisP;

			thisP->i_ID = 1;

			thisP->d_Radius = 0.005;
			thisP->d_Volume = 4.0/3.0 * _PI * glm::pow(thisP->d_Radius, 3.0);

			thisP->d_Mass  = 2700.0 * thisP->d_Volume;

			thisP->d_ElasticModulus = 70.0e9;
			thisP->d_PoissonRatio = 0.25;

			double dx = _RANDOM(0.0,0.1,100);
			double dy = _RANDOM(0.0,0.1,100);
			double dz = 0.0;
			thisP->d3_Position = glm::dvec3(dx, dy, dz);
			thisP->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);
		}
		for(unsigned int index_P = 0; index_P < thisParticleDomain.size(); index_P++)
		{// send to MP vectors
			Particle_CC *thisP = thisParticleDomain[index_P];
			// all MPs
			allParticle.push_back(thisP);
		}
	}

	double dMass_Domain = 0.0;
	for(unsigned int index_P = 0; index_P < allParticle.size(); index_P++)
	{// calculate debug values
		dMass_Domain += allParticle[index_P]->d_Mass;
	}

	a_Runtime.fill(0.0);
	d_DampingCoefficient = 0.00;

	d_TimeIncrement_Maximum = 5.0e-9;
	d_TimeEnd = 1.0;
	d_TimeConsole_Interval = 5.0e-5;

	std::string sDescription = "";
	{
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time (&rawtime);
		timeinfo = localtime(&rawtime);

		strftime(buffer,80,"%d-%m-%Y %H:%M:%S",timeinfo);
		std::string strTime(buffer);

		sDescription += "-------------------------------------------------------------\n";
		sDescription += "-------------------------------------------------------------\n";
		sDescription += "Process started on: " + strTime + "\n";
		sDescription += "-------------------------------------------------------------\n";
//		sDescription += "Number of threads: " + Script(_MAX_N_THREADS) + "\n";
		sDescription += "Time increment: " + Script(d_TimeIncrement_Maximum, 6) + "\n";
		sDescription += "Particle count: " + Script(allParticle.size()) + "\n";
		sDescription += "Mass: " + Script(dMass_Domain,6) + "\n";
		sDescription += "-------------------------------------------------------------\n";

		sDescription += "Global Damping: " + Script(d_DampingCoefficient, 3) + "\n";
	}

	d_TimeConsole_Last = 0.0;
	this->reportConsole(sDescription);
}
// ----------------------------------------------------------------------------
