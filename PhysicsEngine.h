#ifndef PHYSICSENGINE_H
#define PHYSICSENGINE_H

#include <math.h>
#include <vector>
#include <algorithm>
#include <time.h>

#include <omp.h>

#include "Definitions.h"
////#include ".\include\glm\glm.hpp" // windows
//#include "./include/glm/glm.hpp"//sina, glm is a column major implementation

#define _MAX_N_THREADS 1

#define _MAX_PARTICLE 1000
#define _MAX_PPC_PER_PARTICLE 6
#define _MAX_PWC_PER_PARTICLE 3

#include "Particle_CC.h"
#include "Wall_CC.h"

#include "Link_PPC_CC.h"
#include "Link_PWC_CC.h"

#include "ConstitutiveRelation.h"
#include "TimeLine.h"


class PhysicsEngine
{
	public:
		PhysicsEngine();
		virtual ~PhysicsEngine();

		void	initializeWorld_Random(void);

		glm::dvec3 d3_Size_World = glm::dvec3(0.0, 0.0, 0.0);

		// MPM ----------------------------------------------------------------
		int		runSimulation(double dTimeIncrement_Total);

		// methods to communicate with outside -------------------------------
		double getTime_Runtime(void) {return(d_Runtime_Total);}
		double getTime_Current(void) {return(d_Time);}
		double getTime_End(void) {return(d_TimeEnd);}
		double getTime_Increment(void) {return(d_TimeIncrement_Maximum);}
		double getTime_ConsoleInterval(void) {return(d_TimeConsole_Interval);}

		// graphics interface -------------------------------------------------
		unsigned int	getCount_Particles(void) {return(v_Particle.size());}
		unsigned int	getCount_Walls(void) {return(v_Wall.size());}

		std::vector<Particle_CC *>	getParticles(void) {return(v_Particle);}
		std::vector<Wall_CC *>		getWalls(void) {return(v_Wall);}
	protected:
		double d_DampingCoefficient = 0.0;

		// counters
		unsigned long int i_Particle_Count = 0;
		unsigned long int i_Wall_Count = 0;
		unsigned long int i_Link_PPC_Count = 0;
		unsigned long int i_Link_PWC_Count = 0;

		// particles and walls
		std::vector<Particle_CC *> v_Particle;
		std::vector<Wall_CC *> v_Wall;
		// links
		std::vector<Link_PPC_CC *> v_Link_PPC;
		std::vector<Link_PWC_CC *> v_Link_PWC;

		// time related variables
		double d_Time = 0.0; // current simulation time
		double d_TimeIncrement_Maximum = 1.0e-3;
		double d_TimeEnd = 10.0;//5.0e-4;
		int i_TimeCycle = 0;

		double d_Runtime_Total = 0.0;
		std::array<double, 8> a_Runtime;

		double d_TimeConsole_Interval = 500.0*d_TimeIncrement_Maximum;
		double d_TimeConsole_Last = -1.0e12; // before creation

		std::string str_LogFile = "log_.txt";
		void reportConsole(std::string sDescription = "");
	private:
};

#endif // PHYSICSENGINE_H
