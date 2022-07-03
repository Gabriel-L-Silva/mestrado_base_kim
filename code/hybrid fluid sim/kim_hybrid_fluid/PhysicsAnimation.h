#ifndef __PhysicsAnimation_h
#define __PhysicsAnimation_h

#include <iostream>
#include "math/Real.h"

namespace cg
{

/**
* Frame struct that provides relevant information about a simulation frame.
* 
* This struct provides relevant information about a simulation frame, such as
* the frame index and the frame time.
*/
struct Frame final
{
  int index = 0; ///< Frame index.

  double timeIntervalInSeconds = 1.0 / 60.0; ///< Frame time interval.

  /** Default constructor. */
  Frame();

  /** Constructs frame using \p index and optionally its time interval. */
  Frame(int, double);

  /** \returns the frame time interval. */
  double timeInSeconds() const;

  /** Increments the frame index. */
  void advance();

  /** Adds \p delta units to the frame index. */
  void advance(unsigned int delta);

  /** Prefix increment operator does the same as Frame::advance. */
  Frame& operator++(); // prefix

  /** Postfix increment operator does the same as Frame::advance. */
  Frame operator++(int i); // postfix

}; // Frame

/**
* Abstract base class for physics based animations.
* 
* This class implements the top-level interface for physics based animations.
* For underlying simulations, fixed and adaptive subtimestepping are available
* as when performing time-intergration this practice may produce better
* results. When inheriting from this class, subclasses must provide a way to
* compute the next simulation state by overriding the
* PhysicsAnimation::onAdvanceTimeStep method. Initialization routines should
* go in PhysicsAnimation::initialize. This class performs the time-integrations
* by successively advancing frames via PhysicsAnimation::advanceFrame.
*/
class PhysicsAnimation
{
public:
  /** Default constructor. */
  PhysicsAnimation();

  /** Default destructor. */
  virtual ~PhysicsAnimation();

  /** \returns the frame. */
  Frame frame() const;

  /** Sets the given frame. */
  void setFrame(const Frame& frame);

  /** Advances the simulation state. */
  void advanceFrame(const Frame& frame);

  /** \returns the number of fixed sub-timesteps. */
  auto numberOfSubTimeSteps() const { return _numberOfFixedSubTimeSteps; }

  /** Sets the number of fixed sub-timesteps. */
  void setNumberOfSubTimeSteps(size_t subTimeSteps)
  {
    _numberOfFixedSubTimeSteps = math::max<size_t>(subTimeSteps, 1ULL);
  }

  /** \returns \c true if using fixed sub-timestepping, \c false otherwise. */
  bool isUsingFixedSubTimeSteps() const { return _usingFixedSubTimeSteps; }

  /**
  * \brief Sets the sub-timestepping approach.
  * 
  * If \p enable is \c true uses fixed sub-timestepping, otherwise uses
  * adaptive sub-timestepping.
  */
  void setIsUsingFixedSubTimeSteps(bool enable) { _usingFixedSubTimeSteps = enable; }

  /** \returns the current simulation time. */
  auto currentTime() const { return _currentTime; }

protected:
  /**
  * Returns the required number of sub-timesteps for given time interval.
  * 
  * The required number of sub-timesteps can be different depending on the
  * model used. Override this method to implement your own logic for your
  * model specific sub-timestepping for a given interval.
  * 
  * \param[in] timeInterval The time interval in seconds.
  * \return The required number of sub-timesteps.
  */
  virtual size_t numberOfSubTimeSteps(double timeInterval) const;

  /**
  * Called at the first frame to initialize the physics state.
  * 
  * Subclasses can override this method to setup their initial state.
  */
  virtual void initialize() = 0;

  /**
  * Called when a single time-step should be advanced.
  * 
  * When PhysicsAnimation::advanceFrame is called, this class will internally
  * subdivide a frame into sub-steps if needed. Each sub-step, or time-step
  * is then taken to move forward in time. This function is called for each
  * time-step, and a subclass that inherits PhysicsAnimation class should
  * implement this function for its own physics model.
  * 
  * \param[in] timeInterval The time in seconds.
  */
  virtual void onAdvanceTimeStep(double timeInterval) = 0;

private:
  /** Simulation frame. */
  Frame _frame;
  /** Indicates if is using fixed sub-timestepping. */
  bool _usingFixedSubTimeSteps = true;
  /** Number of fixed sub-timestepping. */
  size_t _numberOfFixedSubTimeSteps = 1ULL;
  /** Current simulation time. */
  double _currentTime = 0.0;

  /**
  * Called by PhysicsAnimation::advanceFrame to subdivide the time-step and
  * then invoke PhysicsAnimation::onAdvanceTimeStep to advance the simulation
  * state.
  */
  void advanceTimeStep(double timeInterval);

}; // PhysicsAnimation

} // end namespace cg

#endif // __PhysicsAnimation_h
