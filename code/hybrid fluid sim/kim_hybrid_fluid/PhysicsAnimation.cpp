#include "PhysicsAnimation.h"
#include "utils/Stopwatch.h"

namespace cg
{

#ifdef _DEBUG
#define debug(...) fprintf(stdout, __VA_ARGS__)
#else
#define debug(...) 
#endif

Frame::Frame() {
  // do nothing
}

Frame::Frame(int index, double timeIntervalInSeconds = 1.0 / 60.0) :
  index(index), timeIntervalInSeconds(timeIntervalInSeconds)
{ /* do nothing */ }

Frame& Frame::operator++()
{
  advance();
  return *this;
}

Frame Frame::operator++(int i) // postfix
{
  Frame ret = *this;
  advance();
  return ret;
}

double Frame::timeInSeconds() const
{
  return index * timeIntervalInSeconds;
}

void Frame::advance()
{
  ++index;
}

void Frame::advance(unsigned int delta)
{
  index += delta;
}

PhysicsAnimation::PhysicsAnimation()
{
  _frame.index = -1;
}

PhysicsAnimation::~PhysicsAnimation()
{
  // do nothing
}

Frame
PhysicsAnimation::frame() const
{
  return _frame;
}

void
PhysicsAnimation::setFrame(const Frame& frame)
{
  _frame = frame;
}

size_t
PhysicsAnimation::numberOfSubTimeSteps(double timeInterval) const
{
  return _numberOfFixedSubTimeSteps;
}

void
PhysicsAnimation::advanceFrame(const Frame& frame)
{
  if (frame.index > _frame.index)
  {
    if (_frame.index < 0)
      initialize();

    int numberOfFrames = frame.index - _frame.index;
    for (auto i = 0; i < numberOfFrames; ++i) {
      advanceTimeStep(frame.timeIntervalInSeconds);
    }

    _frame = frame;
  }
}

void
PhysicsAnimation::advanceTimeStep(double timeInterval)
{
  _currentTime = _frame.timeInSeconds();

  if (_usingFixedSubTimeSteps)
  {
    debug("Using fixed sub-timesteps: %llu\n" , _numberOfFixedSubTimeSteps);

    const auto actualTimeInterval = timeInterval /
      static_cast<double>(_numberOfFixedSubTimeSteps);

    for (size_t i = 0; i < _numberOfFixedSubTimeSteps; ++i)
    {
      Stopwatch s;
      s.start();
      onAdvanceTimeStep(actualTimeInterval);
      debug("[INFO] End onAdvanceTimeStep: %lld ms\n", s.lap());

      _currentTime += actualTimeInterval;
    }
  }
  else
  {
    debug("Using adaptive sub-timesteps\n");

    // Adaptive time-stepping
    auto remainingTime = timeInterval;
    while (math::isPositive(remainingTime))
    {
      Stopwatch s;

      auto numSteps = numberOfSubTimeSteps(remainingTime);
      auto actualTimeInterval = remainingTime / static_cast<double>(numSteps);
      
      debug("Number of remaining sub-timesteps: %llu\n", numSteps);

      s.start();
      onAdvanceTimeStep(actualTimeInterval);

      debug("[INFO] End onAdvanceTimeStep: %lld ms\n", s.lap());

      remainingTime -= actualTimeInterval;
      _currentTime += actualTimeInterval;
    }
  }
}


} // end namespace cg
