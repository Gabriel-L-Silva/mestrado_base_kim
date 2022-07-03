#ifndef __KTXHelper_h
#define __KTXHelper_h

#include <cinttypes>

namespace cg
{ // begin namespace cg

namespace ktx
{ // begin namespace ktx

uint32_t load(const char* filename, uint32_t tex = 0, uint32_t unitIdx = 0);

} // end namespace ktx

} // end namespace cg

#endif // __KTXHelper_h
