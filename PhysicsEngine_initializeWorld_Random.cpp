#include "Particle_CC.h"
#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
void PhysicsEngine::initializeWorld_Random(void)
{
	// size of world
	d3_Size_World = glm::dvec3(0.100, 0.100, 0.100);

	if(true)
	{// walls -----------------------------------------------------------------
		std::vector<Wall_CC *> thisWallDomain;

		// bottom wall
		Wall_CC *newW_Bottom = new Wall_CC();

		newW_Bottom->i_ID = 1;

		newW_Bottom->d_ElasticModulus = 200.0e9;
		newW_Bottom->d_PoissonRatio = 0.3;

		newW_Bottom->d3_UnitNormal = glm::dvec3(0.0, 1.0, 0.0);

		newW_Bottom->d3_Position = glm::dvec3(0.5*d3_Size_World.x, 0.0, 0.5*d3_Size_World.z);
		newW_Bottom->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);

		v_Wall.push_back(newW_Bottom);

		// top wall
		Wall_CC *newW_Top = new Wall_CC();

		newW_Top->i_ID = 0;

		newW_Top->d_ElasticModulus = 200.0e9;
		newW_Top->d_PoissonRatio = 0.3;

		newW_Top->d3_UnitNormal = glm::dvec3(0.0, -1.0, 0.0);

		newW_Top->d3_Position = glm::dvec3(0.5*d3_Size_World.x, d3_Size_World.y, 0.5*d3_Size_World.z);
		newW_Top->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);

		v_Wall.push_back(newW_Top);

		// side walls
		for(double dAngle = 0.0; dAngle < 2.0*_PI; dAngle += 0.5*_PI)
		{
			Wall_CC *newW = new Wall_CC();

			newW->i_ID = 0;

			newW->d_ElasticModulus = 200.0e9;
			newW->d_PoissonRatio = 0.3;

			newW->d3_UnitNormal = -glm::normalize(glm::dvec3(glm::cos(dAngle), 0.0, glm::sin(dAngle)));

			newW->d3_Position = 0.5*d3_Size_World + 0.5*d3_Size_World * glm::dvec3(glm::cos(dAngle), 0.0, glm::sin(dAngle));
			newW->d3_Velocity = glm::dvec3(0.0, 0.0, 0.0);

			v_Wall.push_back(newW);
		}
	}

	if(true)
	{// particles -------------------------------------------------------------
		std::vector<Particle_CC *> thisParticleDomain;
		for(unsigned int index_P = 0; index_P < 2; index_P++)
		{// assign material point initial values
			// trial, create a random particle
			double dR = 0.005;
			double dx = 0.5*d3_Size_World.x + index_P*0.2*dR;//_RANDOM(1.0*dR,d3_Size_World.x-1.0*dR,100);
			double dy = 1.0*dR + index_P*2.1*dR;//_RANDOM(1.0*dR,d3_Size_World.y-1.0*dR,100);
			double dz = 0.5*d3_Size_World.z;//_RANDOM(1.0*dR,d3_Size_World.z-1.0*dR,100);

			bool bContact = false;
			for(int index_P_Other = 0; index_P_Other < thisParticleDomain.size(); index_P_Other++)
			{
				Particle_CC *otherP = thisParticleDomain[index_P_Other];

				double dDistance = glm::length(glm::dvec3(dx,dy,dz) - otherP->d3_Position);

				if(dDistance < dR + otherP->d_Radius)
					bContact = true;
			}

			if(bContact == true)
				continue;

			Particle_CC *newP = new Particle_CC();
			thisParticleDomain.push_back(newP);

			newP->i_ID = index_P;

			newP->d_Radius = dR;
			newP->d_Volume = 4.0/3.0 * _PI * glm::pow(newP->d_Radius, 3.0);

			newP->d_Mass  = 2700.0 * newP->d_Volume;
			newP->d_MomentInertia = 0.4*(newP->d_Mass)*glm::pow(dR,2.0);

			newP->d_ElasticModulus = 70.0e9;
			newP->d_PoissonRatio = 0.25;

			newP->d_DampingCoefficient = 1.0e9 * newP->d_Volume;

			newP->d3_Position = glm::dvec3(dx, dy, dz);
			newP->d3_Velocity = glm::dvec3(0., 0.0, 0.0);
			newP->d3_Force_External = newP->d_Mass * glm::dvec3(0.0,-10.0,0.0);

			v_Particle.push_back(newP);
		}
		if(v_Particle.size() > _MAX_PARTICLE)
		{
			std::cout << "PhysicsEngine::initializeWorld, exceeding allowable particle count." << std::endl;
		}
	}

	double dMass_Domain = 0.0;
	for(unsigned int index_P = 0; index_P < v_Particle.size(); index_P++)
	{// calculate debug values
		dMass_Domain += v_Particle[index_P]->d_Mass;
	}

	a_Runtime.fill(0.0);
	d_DampingCoefficient = 0.00;

	d_TimeIncrement_Maximum = 1.0e-6;
	d_TimeEnd = 1.0;
	d_TimeConsole_Interval = 1.0e-2;

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
		sDescription += "Number of threads: " + Script(_MAX_N_THREADS) + "\n";
		sDescription += "Time increment: " + Script(d_TimeIncrement_Maximum, 6) + "\n";
		sDescription += "Particle count: " + Script(v_Particle.size()) + "\n";
		sDescription += "Mass: " + Script(dMass_Domain,6) + "\n";
		sDescription += "-------------------------------------------------------------\n";

		sDescription += "Global Damping: " + Script(d_DampingCoefficient, 3) + "\n";
	}

	d_TimeConsole_Last = 0.0;
	this->reportConsole(sDescription);
}
// ----------------------------------------------------------------------------
