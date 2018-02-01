#include "PhysicsEngine.h"

// ----------------------------------------------------------------------------
int PhysicsEngine::runSimulation(double dTimeIncrement_Total)
{
	omp_set_num_threads(_MAX_N_THREADS);

	clock_t clockCurrent_Total;
	clock_t clockCurrent;

	double dRuntime_MP = 0.0;
	double dRuntime_Block = 0.0;

	double dTimeIncrement_Accumulated = 0.0;
	while(dTimeIncrement_Accumulated < dTimeIncrement_Total)
	{
		// default time-increment value
		double dTimeIncrement = d_TimeIncrement_Maximum;
		// unless reaching the end, where the remaining increment is less than the maximum
		if(d_TimeIncrement_Maximum > (dTimeIncrement_Total-dTimeIncrement_Accumulated))
			dTimeIncrement = dTimeIncrement_Total - dTimeIncrement_Accumulated;
		// calculate increments accumulated
		dTimeIncrement_Accumulated += dTimeIncrement;

		clockCurrent_Total = clock();
		dRuntime_MP = omp_get_wtime();
		#pragma omp parallel
		{
			int iThread_Count = omp_get_num_threads();
			int	iThread_This = omp_get_thread_num();

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// reset particle variables: force
			#pragma omp for
			for(unsigned int index_P = 0; index_P < v_Particle.size(); index_P++)
			{
				Particle_CC *thisP = v_Particle[index_P];

				thisP->b_Marked = false;

				thisP->d3_Force = thisP->d3_Force_External;
			}

			// reset PPC counter
			i_Link_PPC_Count = 0;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// find PP contacts
			// sina, this will need locks
			#pragma omp for
			for(unsigned int index_P1 = 0; index_P1 < v_Particle.size(); index_P1++)
			{
				Particle_CC *thisP1 = v_Particle[index_P1];

				for(unsigned int index_P2 = index_P1 + 1; index_P2 < v_Particle.size(); index_P2++)
				{
					Particle_CC *thisP2 = v_Particle[index_P2];

					double dRadius_P1 = thisP1->d_Radius;
					double dRadius_P2 = thisP2->d_Radius;

					glm::dvec3 d3Position_P1 = thisP1->d3_Position;
					glm::dvec3 d3Position_P2 = thisP2->d3_Position;

					glm::dvec3 d3Velocity_P1 = thisP1->d3_Velocity;
					glm::dvec3 d3Velocity_P2 = thisP2->d3_Velocity;

					glm::dvec3 d3Relative_Position = d3Position_P2 - d3Position_P1;
					glm::dvec3 d3Relative_Velocity = d3Velocity_P2 - d3Velocity_P1;

					double dDistance = glm::length(d3Relative_Position);

					if(dDistance < dRadius_P2 + dRadius_P1)
					{
						thisP1->b_Marked = true;
						thisP2->b_Marked = true;

						v_Link_PPC[i_Link_PPC_Count]->index_P1 = index_P1;
						v_Link_PPC[i_Link_PPC_Count]->index_P2 = index_P2;

						i_Link_PPC_Count++;
					}
				}
			}

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// calculate PP contacts
			#pragma omp for
			for(unsigned int index_PPC = 0; index_PPC < i_Link_PPC_Count; index_PPC++)
			{
				// initial values
				unsigned int index_P1 = v_Link_PPC[index_PPC]->index_P1;
				unsigned int index_P2 = v_Link_PPC[index_PPC]->index_P2;

				Particle_CC *thisP1 = v_Particle[index_P1];
				Particle_CC *thisP2 = v_Particle[index_P2];

				double dRadius_P1 = thisP1->d_Radius;
				double dRadius_P2 = thisP2->d_Radius;

				double dE_P1 = thisP1->d_ElasticModulus;
				double dE_P2 = thisP2->d_ElasticModulus;

				double dP_P1 = thisP1->d_PoissonRatio;
				double dP_P2 = thisP2->d_PoissonRatio;

				glm::dvec3 d3Position_P1 = thisP1->d3_Position;
				glm::dvec3 d3Position_P2 = thisP2->d3_Position;

				glm::dvec3 d3Velocity_P1 = thisP1->d3_Velocity;
				glm::dvec3 d3Velocity_P2 = thisP2->d3_Velocity;

				glm::dvec3 d3Relative_Position = d3Position_P2 - d3Position_P1;
				glm::dvec3 d3Relative_Velocity = d3Velocity_P2 - d3Velocity_P1;

				// calculated values
				//http://en.wikipedia.org/wiki/Contact_mechanics
				double dDistance = glm::length(d3Relative_Position);

				glm::dvec3 d3UnitNormal = glm::normalize(d3Relative_Position);

				double dDeformation_Normal = (dRadius_P1 + dRadius_P2) - dDistance;

				// normal contact force
				//hertzian parameters
				double dE_Hertzian = 0.0;
				dE_Hertzian += (1.0-dP_P1*dP_P1)/dE_P1;
				dE_Hertzian += (1.0-dP_P2*dP_P2)/dE_P2;
				dE_Hertzian = 1.0/dE_Hertzian;

				double dRadius_Hertzian = 0.0;
				dRadius_Hertzian = 1.0/dRadius_P1 + 1.0/dRadius_P2;
				dRadius_Hertzian = 1.0/dRadius_Hertzian;

				double dStiffness_Hertzian = 0.0;
				dStiffness_Hertzian = 4.0/3.0 * sqrt(dRadius_Hertzian) * dE_Hertzian;

				glm::dvec3 d3Force_Contact_Normal = glm::dvec3(0.0, 0.0, 0.0);
				d3Force_Contact_Normal = d3UnitNormal;
				d3Force_Contact_Normal *= (dStiffness_Hertzian * pow(dDeformation_Normal,1.5));

				thisP1->d3_Force -= d3Force_Contact_Normal;
				thisP2->d3_Force += d3Force_Contact_Normal;

				// normal damping force
				double dC_P1 = thisP1->d_DampingCoefficient;
				double dC_P2 = thisP2->d_DampingCoefficient;
				double dC = 0.5*(dC_P1 + dC_P2);

				glm::dvec3 d3Force_Damping_Normal = glm::dvec3(0.0, 0.0, 0.0);

				d3Force_Damping_Normal = d3UnitNormal;
				d3Force_Damping_Normal *= glm::dot(d3Relative_Velocity, d3UnitNormal);
				d3Force_Damping_Normal *= dC;

				thisP1->d3_Force += d3Force_Damping_Normal;
				thisP2->d3_Force -= d3Force_Damping_Normal;
			}

			// reset PWC counter
			i_Link_PWC_Count = 0;

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// find PW contacts
			// sina, this will need locks
			#pragma omp for
			for(unsigned int index_P = 0; index_P < v_Particle.size(); index_P++)
			{
				Particle_CC *thisP = v_Particle[index_P];

				for(unsigned int index_W = 0; index_W < v_Wall.size(); index_W++)
				{
					Wall_CC *thisW = v_Wall[index_W];

					double dRadius_P = thisP->d_Radius;

					glm::dvec3 d3Position_P = thisP->d3_Position;
					glm::dvec3 d3Position_W = thisW->d3_Position;

					glm::dvec3 d3Velocity_P = thisP->d3_Velocity;
					glm::dvec3 d3Velocity_W = thisW->d3_Velocity;

					glm::dvec3 d3Relative_Position = d3Position_P - d3Position_W;
					glm::dvec3 d3Relative_Velocity = d3Velocity_P - d3Velocity_W;

					glm::dvec3 d3UnitNormal = thisW->d3_UnitNormal;

					double dNormalDistance = glm::length(glm::dot(d3Relative_Position, d3UnitNormal));

					if(dNormalDistance < dRadius_P)
					{
						thisP->b_Marked = true;

						v_Link_PWC[i_Link_PWC_Count]->index_P = index_P;
						v_Link_PWC[i_Link_PWC_Count]->index_W = index_W;

						i_Link_PWC_Count++;
					}
				}
			}

			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// calculate PW contacts
			#pragma omp for
			for(unsigned int index_PWC = 0; index_PWC < i_Link_PWC_Count; index_PWC++)
			{
				// initial values
				unsigned int index_P = v_Link_PWC[index_PWC]->index_P;
				unsigned int index_W = v_Link_PWC[index_PWC]->index_W;

				Particle_CC *thisP = v_Particle[index_P];
				Wall_CC		*thisW = v_Wall[index_W];

				glm::dvec3 d3Position_P = thisP->d3_Position;
				glm::dvec3 d3Position_W = thisW->d3_Position;

				glm::dvec3 d3Velocity_P = thisP->d3_Velocity;
				glm::dvec3 d3Velocity_W = thisW->d3_Velocity;

				glm::dvec3 d3Relative_Position = d3Position_P - d3Position_W;
				glm::dvec3 d3Relative_Velocity = d3Velocity_P - d3Velocity_W;

				glm::dvec3 d3UnitNormal = thisW->d3_UnitNormal;

				double dNormalDistance = glm::length(glm::dot(d3Relative_Position, d3UnitNormal));

				double dDeformation_Normal = thisP->d_Radius - dNormalDistance;

				//normal contact force
				//http://en.wikipedia.org/wiki/Contact_mechanics
				double dForce	= 0.0;
				double dStiffness_Hertzian = 0.0;

				double dE_P = thisP->d_ElasticModulus;
				double dE_W = thisW->d_ElasticModulus;

				double dP_P = thisP->d_PoissonRatio;
				double dP_W = thisW->d_PoissonRatio;

				double dE_Hertzian = 0.0;
				dE_Hertzian += (1.0-dP_P*dP_P)/dE_P;
				dE_Hertzian += (1.0-dP_W*dP_W)/dE_W;
				dE_Hertzian = 1.0/dE_Hertzian;

				double dRadius_Hertzian = thisP->d_Radius;

				dStiffness_Hertzian = 4.0/3.0 * sqrt(dRadius_Hertzian) * dE_Hertzian;

				glm::dvec3 d3Force_Contact_Normal = glm::dvec3(0.0, 0.0, 0.0);
				d3Force_Contact_Normal = d3UnitNormal;
				d3Force_Contact_Normal *= (dStiffness_Hertzian * pow(dDeformation_Normal,1.5));

				thisP->d3_Force += d3Force_Contact_Normal;
				thisW->d3_Force -= d3Force_Contact_Normal;

				// normal damping force
				double dC = thisP->d_DampingCoefficient;

				glm::dvec3 d3Force_Damping_Normal = glm::dvec3(0.0, 0.0, 0.0);
				d3Force_Damping_Normal = d3UnitNormal;
				d3Force_Damping_Normal *= -glm::dot(d3Relative_Velocity, d3UnitNormal);
				d3Force_Damping_Normal *= dC;

				thisP->d3_Force += d3Force_Damping_Normal;
			}
/*
			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// PW contact
			#pragma omp for
			for(unsigned int index_P = 0; index_P < v_Particle.size(); index_P++)
			{
				Particle_CC *thisP = v_Particle[index_P];

				for(unsigned int index_W = 0; index_W < v_Wall.size(); index_W++)
				{
					Wall_CC *thisW = v_Wall[index_W];

					glm::dvec3 d3Relative_Position = thisP->d3_Position - thisW->d3_Position;
					glm::dvec3 d3Relative_Velocity = thisP->d3_Velocity - thisW->d3_Velocity;

					double dNormalDistance = glm::length(glm::dot(d3Relative_Position, thisW->d3_UnitNormal));

					bool bContact = false;
					if(dNormalDistance < thisP->d_Radius)
					{
						bContact = true;
						thisP->b_Marked = true;
					}
					else
					{
						thisP->b_Marked = false;
					}

					if(bContact == true)
					{// reaction force
						glm::dvec3 d3UnitNormal = thisW->d3_UnitNormal;
						double dOverlap = thisP->d_Radius - dNormalDistance;

						//calculate normal forces
						//http://en.wikipedia.org/wiki/Contact_mechanics
						double dForce	= 0.0;
						double dStiffness_Hertzian = 0.0;

						double dE_P = thisP->d_ElasticModulus;
						double dE_W = thisW->d_ElasticModulus;

						double dP_P = thisP->d_PoissonRatio;
						double dP_W = thisW->d_PoissonRatio;

						//hertzian modulus
						double dE_Hertzian = 0.0;
						dE_Hertzian += (1.0-dP_P*dP_P)/dE_P;
						dE_Hertzian += (1.0-dP_W*dP_W)/dE_W;
						dE_Hertzian = 1.0/dE_Hertzian;

						//hertzian radius
						double dRadius_Hertzian = thisP->d_Radius;

						//hertzian stiffness
						dStiffness_Hertzian = 4.0/3.0 * sqrt(dRadius_Hertzian) * dE_Hertzian;
						dForce	= dStiffness_Hertzian * pow(dOverlap,1.5);

						// force
						glm::dvec3 d3Force = d3UnitNormal * dForce;

						thisP->d3_Force += d3Force;
					}
					if(bContact == true)
					{// damping force
						glm::dvec3 d3UnitNormal = thisW->d3_UnitNormal;
						double dC = thisP->d_DampingCoefficient;

						glm::dvec3 d3Force = glm::dvec3(0.0, 0.0, 0.0);

						d3Force = d3UnitNormal;
						d3Force *= -glm::dot(d3Relative_Velocity, d3UnitNormal);
						d3Force *= (dC * (thisP->d_Volume));

						thisP->d3_Force += d3Force;
					}
				}
			}
*/
/*
			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// PP contact
			#pragma omp for
			for(unsigned int index_P = 0; index_P < v_Particle.size(); index_P++)
			{
				Particle_CC *thisP = v_Particle[index_P];

				for(unsigned int index_P_Other = index_P + 1; index_P_Other < v_Particle.size(); index_P_Other++)
				{
					Particle_CC *thisP_Other = v_Particle[index_P_Other];

					glm::dvec3 d3Relative_Position = thisP_Other->d3_Position - thisP->d3_Position;
					glm::dvec3 d3Relative_Velocity = thisP_Other->d3_Velocity - thisP->d3_Velocity;

					double dDistance = glm::length(thisP->d3_Position - thisP_Other->d3_Position);

					bool bContact = false;
					if(dDistance < thisP->d_Radius + thisP_Other->d_Radius)
					{
						bContact = true;
						thisP->b_Marked = true;
					}
					else
					{
						thisP->b_Marked = false;
					}

					if(bContact == true)
					{
						double dR = thisP->d_Radius;
						double dE = thisP->d_ElasticModulus;
						double dP = thisP->d_PoissonRatio;

						double dR_Other = thisP_Other->d_Radius;
						double dE_Other = thisP_Other->d_ElasticModulus;
						double dP_Other = thisP_Other->d_PoissonRatio;

						glm::dvec3 d3UnitNormal = glm::normalize(d3Relative_Position);
						double dOverlap = (dR + dR_Other) - glm::length(d3Relative_Position);

						//calculate normal forces
						//http://en.wikipedia.org/wiki/Contact_mechanics
						//hertzian parameters
						double dE_Hertzian = 0.0;
						dE_Hertzian += (1.0-dP*dP)/dE;
						dE_Hertzian += (1.0-dP_Other*dP_Other)/dE_Other;
						dE_Hertzian = 1.0/dE_Hertzian;

						double dRadius_Hertzian = 0.0;
						dRadius_Hertzian = 1.0/dR + 1.0/dR_Other;
						dRadius_Hertzian = 1.0/dRadius_Hertzian;

						double dStiffness_Hertzian = 0.0;
						dStiffness_Hertzian = 4.0/3.0 * sqrt(dRadius_Hertzian) * dE_Hertzian;

						glm::dvec3 d3Force = glm::dvec3(0.0, 0.0, 0.0);
						d3Force = d3UnitNormal;
						d3Force	*= (dStiffness_Hertzian * pow(dOverlap,1.5));

						thisP->d3_Force			-= d3Force;
						thisP_Other->d3_Force	+= d3Force;
					}
					if(bContact == true)
					{// damping force
						glm::dvec3 d3UnitNormal = glm::normalize(d3Relative_Position);

						double dC = thisP->d_DampingCoefficient;
						double dC_Other = thisP_Other->d_DampingCoefficient;

						glm::dvec3 d3Force = glm::dvec3(0.0, 0.0, 0.0);

						d3Force = d3UnitNormal;
						d3Force *= glm::dot(d3Relative_Velocity, d3UnitNormal);

						thisP->d3_Force += d3Force*(dC*(thisP->d_Volume));
						thisP_Other->d3_Force -= d3Force*(dC_Other*(thisP_Other->d_Volume));
					}
				}
			}
*/
			#pragma omp barrier
			dRuntime_Block = omp_get_wtime();
			// calculate accelerations, velocities, updated positions
			#pragma omp for
			for(unsigned int index_P = 0; index_P < v_Particle.size(); index_P++)
			{
				Particle_CC *thisP = v_Particle[index_P];

				thisP->d3_Acceleration = thisP->d3_Force * (1.0/thisP->d_Mass);

				// update velocity and position
				thisP->d3_Velocity += thisP->d3_Acceleration * dTimeIncrement;
				thisP->d3_Position += thisP->d3_Velocity * dTimeIncrement;
			}
		}

		d_Runtime_Total += omp_get_wtime() - dRuntime_MP;
		//report to console ---------------------------------------------------
		if(d_Time - d_TimeConsole_Last > d_TimeConsole_Interval)
		{
			d_TimeConsole_Last = d_Time;
			this->reportConsole();
		}

		d_Time += dTimeIncrement;
		i_TimeCycle++;
	}

	if(d_Time < d_TimeEnd)
		return(0);
	else
	{
		d_TimeConsole_Last = d_Time;
		this->reportConsole();
		std::cout << "Analysis complete ..." << std::endl;
		return(1);
	}
}
// ----------------------------------------------------------------------------
