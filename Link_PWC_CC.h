#ifndef LINK_PWC_CC_H
#define LINK_PWC_CC_H

#include <iostream>

#include "Definitions.h"

class Link_PWC_CC
{
	public:
		Link_PWC_CC() {;}
		virtual ~Link_PWC_CC() {;}

		unsigned int index_P = 0;
		unsigned int index_W = 0;

		double d_DampingCoefficient = 0.0;

		double d_Deformation_Normal = 0.0;
		double d_Deformation_Tangent = 0.0;

		double d_DeformationRate_Normal = 0.0;

		glm::dvec3 d3_Relative_Position = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_Relative_Velocity = glm::dvec3(0.0,0.0,0.0);

		glm::dvec3 d3_UnitNormal = glm::dvec3(0.0,0.0,0.0);// will be equal to the unit normal of the wall
		glm::dvec3 d3_UnitTangent = glm::dvec3(0.0,0.0,0.0);// will be in the direction of the tangent of the velocity of the particle

		glm::dvec3 d3_Force_Normal_P = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_Force_Normal_W = glm::dvec3(0.0,0.0,0.0);

		glm::dvec3 d3_Force_Tangent_P = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3 d3_Force_Tangent_W = glm::dvec3(0.0,0.0,0.0);
	protected:

	private:
};

#endif // LINK_PWC_CC_H
