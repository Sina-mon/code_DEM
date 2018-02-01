#ifndef PARTICLE_CC_H
#define PARTICLE_CC_H

#include <iostream>

#include "Definitions.h"

#define	_LINPOS	1	//0000000001
#define	_LINVEL	2	//0000000010
#define	_LINACC	4	//0000000100
#define	_LINJER	8	//0000001000
#define	_ANGPOS	16	//0000010000
#define	_ANGVEL	32	//0000100000
#define	_ANGACC	64	//0001000000
#define	_ANGJER	128	//0010000000
#define	_MATERI	256	//0100000000
#define	_RADIUS	512	//1000000000

//-----------------------------------------------
class Particle_CC
{
	public:
		Particle_CC() {;}
		virtual ~Particle_CC() {;}

		bool	b_Marked = false;

		int		i_ID;
		double	d_Radius = 0.0;
		double	d_Volume = 0.0;
		double	d_Mass = 0.0;
		double	d_MomentInertia = 0.0;

		double	d_ElasticModulus = 0.0;
		double	d_PoissonRatio = 0.0;
		double	d_DampingCoefficient = 0.0;

		glm::dvec3	d3_Position = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3	d3_Velocity = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3	d3_Acceleration = glm::dvec3(0.0,0.0,0.0);

		glm::dvec3	d3_AngularPosition = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3	d3_AngularVelocity = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3	d3_AngularAcceleration = glm::dvec3(0.0,0.0,0.0);

		glm::dvec3	d3_Force = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3	d3_Force_External = glm::dvec3(0.0,0.0,0.0);

		glm::dvec3	d3_Moment = glm::dvec3(0.0,0.0,0.0);
		glm::dvec3	d3_Moment_External = glm::dvec3(0.0,0.0,0.0);

		std::string		Script			(int) const;
	protected:
	private:
};
//-----------------------------------------------
#endif
