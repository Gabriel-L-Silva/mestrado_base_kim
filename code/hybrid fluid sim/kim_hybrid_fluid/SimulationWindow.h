#ifndef __SimulationWindow_h
#define __SimulationWindow_h

#include "graphics/GLWindow.h"
#include "graphics/GLBuffer.h"
#include "math/Matrix4x4.h"
#include "GLDrawSpheres.h"
//#include <jet/jet.h>
//#include <jet/box2.h>
//#include <jet/flip_solver2.h>
//#include <jet/rigid_body_collider2.h>
//#include <jet/sphere2.h>
//#include <jet/volume_particle_emitter2.h>

class SimulationWindow final: public cg::GLWindow
{
public:
	SimulationWindow(int width, int height);

	~SimulationWindow() {
		/*if (_positions)
			delete _positions;*/
	}

private:
	bool _paused{ true };
	bool _enableColorMap{ false };
	double _radius;
	float _particleSpacing{ 0.005f };
	cg::vec4f _particleColor{cg::Color::blue};
	cg::GLSL::Program _program;
	/*cg::GLBuffer<jet::Vector2D>* _positions;
	cg::GLBuffer<jet::Vector2D>* _velocities;*/
	GLuint _vao;

	//jet::Frame _frame;
	//jet::FlipSolver2Ptr _solver;
	//jet::Surface2Ptr _surface; // used to emit the particles
	//jet::VolumeParticleEmitter2Ptr _emitter; // particle emitter

	void initialize() override;
	void gui() override;
	void render() override;
	void terminate() override;

	void resetSimulation();
	void surfaceOptions();
};


#endif // __SimulationWindow_h
