#include "GraphicsEngine.h"

// ----------------------------------------------------------------------------
void GraphicsEngine::drawGame(void)
{
	glm::vec3 f3Size_World = mpm_PhysicsEngine->d3_Size_World;

	if(true)
	{// create shadow map
		gl_Shadow_Texture->bindRenderTarget();
//		gl_Canvas_Texture->bindRenderTarget(i_ScreenWidth, i_ScreenHeight);
		gl_ShadowProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		GLuint transformationShadowLocation = gl_ShadowProgram.getUniformLocation("transformationShadowMatrix");

		std::vector<Particle_CC *> vParticles = mpm_PhysicsEngine->getParticles();

		for(int index_MP = 0; index_MP < vParticles.size(); index_MP++)
		{
			Particle_CC *thisP = vParticles[index_MP];
			// particle position
			float fSize = thisP->d_Radius;
			Transformation glTransformation(thisP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection() * glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// draw
			gl_Particle_Mesh->Draw();
		}

		gl_ShadowProgram.unuse();
	}

	if(true)
	{// draw all objects to the canvas (not the screen)
		v_Canvas_Texture[(int)enum_Canvas::MAIN]->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Diffuse_Texture->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		// particles
		std::vector<Particle_CC *> vParticles = mpm_PhysicsEngine->getParticles();
		for(int index_P = 0; index_P < vParticles.size(); index_P++)
		{
			Particle_CC *thisP = vParticles[index_P];

			// position, and size
			float fSize = thisP->d_Radius;
			glm::vec3 f3Position = thisP->d3_Position;

			// color
			glm::vec4 f4objectColor = _RED;
			if(thisP->b_Marked == true)
				f4objectColor = _GREEN;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(f3Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}

		// walls
		std::vector<Wall_CC *> vWalls = mpm_PhysicsEngine->getWalls();
		for(int index_W = 0; index_W < vWalls.size(); index_W++)
		{
			Wall_CC *thisW = vWalls[index_W];

			// position, and size
			//glm::vec3 f3Size = 0.1f*f3Size_World * glm::vec3(thisW->d3_UnitNormal + 0.05*glm::dvec3(1.0,1.0,1.0));
			glm::vec3 f3Size = 0.1f*f3Size_World * glm::vec3(glm::dvec3(1.0,1.0,1.0) - glm::abs(thisW->d3_UnitNormal));
			glm::vec3 f3Position = thisW->d3_Position;

			// color
			glm::vec4 f4objectColor = _GRAY;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(f3Position, glm::vec3(0.0, 0.0, 0.0), f3Size);
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}
		gl_BasicProgram.unuse();
	}

	if(true)
	{// bind the screen for final output
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0); // drawing to the window

		// which gl program to use
		gl_FinalProgram.use();

		// clear
		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		float fScreenRatio = 1.0 / (int)enum_Canvas::COUNT;
		glm::vec2 f2ScreenRatio = glm::vec2((float)i_ScreenHeight/i_ScreenWidth,1.0);

		if((int)enum_Canvas::MAIN < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::MAIN]->bindTextureUnit(0);
			// set viewport
			glm::vec2 f2PositionRatio = glm::vec2((float)(int)enum_Canvas::MAIN/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			gl_Canvas_Mesh->Draw();
		}

		// undind texture
		gl_FinalProgram.unuse();
	}

	SDL_GL_SwapWindow(p_Window);
}
// ----------------------------------------------------------------------------
