#ifndef WALL_CC_H
#define WALL_CC_H

#include "Definitions.h"

//-----------------------------------------------
class Wall_CC
{
	public:
		Wall_CC() {;}
		virtual ~Wall_CC() {;}

		int		i_ID;

		double	d_ElasticModulus = 0.0;
		double	d_PoissonRatio = 0.0;

		glm::dvec3	d3_UnitNormal; // outward normal, pointing toward air
		glm::dvec3	d3_Position;
		glm::dvec3	d3_Velocity;
	protected:
	private:
};
//-----------------------------------------------
#endif
