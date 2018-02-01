#ifndef LINK_PPC_CC_H
#define LINK_PPC_CC_H

#include <iostream>

#include "Definitions.h"

#include "Particle_CC.h"

class Link_PPC_CC
{
	public:
		Link_PPC_CC() {;}
		virtual ~Link_PPC_CC() {;}

		unsigned int index_P1 = 0;
		unsigned int index_P2 = 0;

		double d_DampingCoefficient = 0.0;

		double d_Deformation_Normal = 0.0;
		double d_Deformation_Tangent = 0.0;

		double d_DeformationRate_Normal = 0.0;

		glm::dvec3 d3_Relative_Position = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_Relative_Velocity = glm::dvec3(0.0,0.0,0.0);

		glm::dvec3 d3_UnitNormal = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_UnitTangent = glm::dvec3(0.0,0.0,0.0);

		glm::dvec3 d3_Force_Normal_P1 = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_Force_Normal_P2 = glm::dvec3(0.0,0.0,0.0);

		glm::dvec3 d3_Force_Tangent_P1 = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_Force_Tangent_P2 = glm::dvec3(0.0,0.0,0.0);
	protected:

	private:
};

#endif // LINK_PPC_CC_H
