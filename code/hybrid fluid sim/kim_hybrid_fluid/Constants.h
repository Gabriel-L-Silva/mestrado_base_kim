#ifndef __Constants_h
#define __Constants_h

/**
* Constants namespace
*/
namespace cg::constants
{

//! No direction.
constexpr int directionNone = 0;

//! Left direction.
constexpr int directionLeft = 1 << 0;

//! Right direction.
constexpr int directionRight = 1 << 1;

//! Down direction.
constexpr int directionDown = 1 << 2;

//! Up direction.
constexpr int directionUp = 1 << 3;

//! Back direction.
constexpr int directionBack = 1 << 4;

//! Front direction.
constexpr int directionFront = 1 << 5;

//! All directions.
constexpr int directionAll = directionLeft | directionRight |
  directionDown | directionUp | directionBack | directionFront;

} // end namespace cg::constants

#endif // __Constants_h
