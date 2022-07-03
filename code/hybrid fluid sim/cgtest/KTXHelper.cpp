#include "KTXHelper.h"
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#else
#define NOMINMAX
#include <GL/gl3w.h>
#endif
#include <cstdio>
#include <cstring>

namespace cg
{ // begin namespace cg

namespace ktx
{ // begin namespace ktx

namespace file
{ // begin namespace file

struct header
{
  unsigned char identifier[12];
  uint32_t endianness;
  uint32_t gltype;
  uint32_t gltypesize;
  uint32_t glformat;
  uint32_t glinternalformat;
  uint32_t glbaseinternalformat;
  uint32_t pixelwidth;
  uint32_t pixelheight;
  uint32_t pixeldepth;
  uint32_t arrayelements;
  uint32_t faces;
  uint32_t miplevels;
  uint32_t keypairbytes;
};

static const unsigned char identifier[] =
{
  0xAB, 0x4B, 0x54, 0x58, 0x20, 0x31, 0x31, 0xBB, 0x0D, 0x0A, 0x1A, 0x0A
};

static const auto
swap32(const uint32_t u32)
{
  union
  {
    uint32_t u32;
    unsigned char u8[4];
  } a, b;

  a.u32 = u32;
  b.u8[0] = a.u8[3];
  b.u8[1] = a.u8[2];
  b.u8[2] = a.u8[1];
  b.u8[3] = a.u8[0];
  return b.u32;
}

static const auto
swap16(const uint16_t u16)
{
  union
  {
    uint16_t u16;
    unsigned char u8[2];
  } a, b;

  a.u16 = u16;
  b.u8[0] = a.u8[1];
  b.u8[1] = a.u8[0];
  return b.u16;
}

static auto
computeStride(const header& h, uint32_t width, uint32_t pad = 4)
{
  uint32_t channels{};

  switch (h.glbaseinternalformat)
  {
    case GL_RED:
      channels = 1;
      break;
    case GL_RG:
      channels = 2;
      break;
    case GL_BGR:
    case GL_RGB:
      channels = 3;
      break;
    case GL_BGRA:
    case GL_RGBA:
      channels = 4;
      break;
  }

  auto stride = h.gltypesize * channels * width;

  stride = (stride + (pad - 1)) & ~(pad - 1);
  return stride;
}

inline static auto
computeFaceSize(const header& h)
{
  return computeStride(h, h.pixelwidth) * h.pixelheight;
}

inline uint32_t
read(FILE* f, uint32_t tex, uint32_t unitIdx)
{
  header h;

  if (fread(&h, sizeof(h), 1, f) != 1)
    return 0;
  if (memcmp(h.identifier, identifier, sizeof(identifier)) != 0)
    return 0;
  if (h.endianness == 0x01020304)
  {
    h.endianness = swap32(h.endianness);
    h.gltype = swap32(h.gltype);
    h.gltypesize = swap32(h.gltypesize);
    h.glformat = swap32(h.glformat);
    h.glinternalformat = swap32(h.glinternalformat);
    h.glbaseinternalformat = swap32(h.glbaseinternalformat);
    h.pixelwidth = swap32(h.pixelwidth);
    h.pixelheight = swap32(h.pixelheight);
    h.pixeldepth = swap32(h.pixeldepth);
    h.arrayelements = swap32(h.arrayelements);
    h.faces = swap32(h.faces);
    h.miplevels = swap32(h.miplevels);
    h.keypairbytes = swap32(h.keypairbytes);
  }
  else if (h.endianness != 0x04030201)
    return 0;

  GLenum target;

  if (h.pixelheight == 0)
    target = h.arrayelements == 0 ? GL_TEXTURE_1D : GL_TEXTURE_1D_ARRAY;
  else if (h.pixeldepth == 0)
  {
    if (h.arrayelements == 0)
      target = h.faces == 0 ? GL_TEXTURE_2D : GL_TEXTURE_CUBE_MAP;
    else
      target = h.faces == 0 ? GL_TEXTURE_2D_ARRAY : GL_TEXTURE_CUBE_MAP_ARRAY;
  }
  else
    target = GL_TEXTURE_3D;
  if (h.pixelwidth == 0 || h.pixelheight == 0 && h.pixeldepth != 0)
    return 0;

  auto temp = tex;

  if (temp == 0)
    glGenTextures(1, &tex);
  glActiveTexture(GL_TEXTURE0 + unitIdx);
  glBindTexture(target, tex);

  auto dataStart = ftell(f) + h.keypairbytes;

  fseek(f, 0, SEEK_END);

  auto dataEnd = ftell(f);
  auto dataSize = dataEnd - dataStart;
  auto data =  new unsigned char[dataSize];

  memset(data, 0, dataSize);
  fseek(f, dataStart, SEEK_SET);
  fread(data, 1, dataSize, f);
  if (h.miplevels == 0)
    h.miplevels = 1;
  switch (target)
  {
    case GL_TEXTURE_1D:
      glTexStorage1D(GL_TEXTURE_1D,
        h.miplevels,
        h.glinternalformat,
        h.pixelwidth);
      glTexSubImage1D(GL_TEXTURE_1D,
        0,
        0,
        h.pixelwidth,
        h.glformat,
        h.glinternalformat,
        data);
      break;
    case GL_TEXTURE_2D:
      if (h.gltype == GL_NONE)
        glCompressedTexImage2D(GL_TEXTURE_2D,
          0,
          h.glinternalformat,
          h.pixelwidth,
          h.pixelheight,
          0,
          420 * 380 / 2,
          data);
      else
      {
        glTexStorage2D(GL_TEXTURE_2D,
          h.miplevels,
          h.glinternalformat,
          h.pixelwidth,
          h.pixelheight);
        {
          auto* ptr = data;
          auto height = h.pixelheight, width = h.pixelwidth;

          glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
          for (uint32_t i = 0; i < h.miplevels; i++)
          {
            glTexSubImage2D(GL_TEXTURE_2D,
              i,
              0,
              0,
              width,
              height,
              h.glformat,
              h.gltype, ptr);
            ptr += height * (size_t)computeStride(h, width, 1);
            if (!(height >>= 1))
              height = 1;
            if (!(width >>= 1))
              width = 1;
          }
        }
      }
      break;
    case GL_TEXTURE_3D:
      glTexStorage3D(GL_TEXTURE_3D,
        h.miplevels,
        h.glinternalformat,
        h.pixelwidth,
        h.pixelheight,
        h.pixeldepth);
      glTexSubImage3D(GL_TEXTURE_3D,
        0,
        0,
        0,
        0,
        h.pixelwidth,
        h.pixelheight,
        h.pixeldepth,
        h.glformat,
        h.gltype,
        data);
      break;
    case GL_TEXTURE_1D_ARRAY:
      glTexStorage2D(GL_TEXTURE_1D_ARRAY,
        h.miplevels,
        h.glinternalformat,
        h.pixelwidth,
        h.arrayelements);
      glTexSubImage2D(GL_TEXTURE_1D_ARRAY,
        0,
        0,
        0,
        h.pixelwidth,
        h.arrayelements,
        h.glformat,
        h.gltype,
        data);
      break;
    case GL_TEXTURE_2D_ARRAY:
      glTexStorage3D(GL_TEXTURE_2D_ARRAY,
        h.miplevels,
        h.glinternalformat,
        h.pixelwidth,
        h.pixelheight,
        h.arrayelements);
      glTexSubImage3D(GL_TEXTURE_2D_ARRAY,
        0,
        0,
        0,
        0,
        h.pixelwidth,
        h.pixelheight,
        h.arrayelements,
        h.glformat,
        h.gltype,
        data);
      break;
    case GL_TEXTURE_CUBE_MAP:
      glTexStorage2D(GL_TEXTURE_CUBE_MAP,
        h.miplevels,
        h.glinternalformat,
        h.pixelwidth,
        h.pixelheight);
      {
        auto fs = computeFaceSize(h);

        for (uint32_t i = 0; i < h.faces; i++)
          glTexSubImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i,
            0,
            0,
            0,
            h.pixelwidth,
            h.pixelheight,
            h.glformat,
            h.gltype,
            data + fs * (size_t)i);
      }
      break;
    case GL_TEXTURE_CUBE_MAP_ARRAY:
      glTexStorage3D(GL_TEXTURE_CUBE_MAP_ARRAY,
        h.miplevels,
        h.glinternalformat,
        h.pixelwidth,
        h.pixelheight,
        h.arrayelements);
      glTexSubImage3D(GL_TEXTURE_CUBE_MAP_ARRAY,
        0,
        0,
        0,
        0,
        h.pixelwidth,
        h.pixelheight,
        h.faces * h.arrayelements,
        h.glformat,
        h.gltype,
        data);
      break;
    default:
      if (temp == 0)
        glDeleteTextures(1, &tex);
      tex = 0;
      goto bad_target;
  }
  if (h.miplevels == 1)
    glGenerateMipmap(target);
bad_target:
  delete []data;
  return tex;
}

} // end namespace file

uint32_t
load(const char* filename, uint32_t tex, uint32_t unitIdx)
{
  auto f = fopen(filename, "rb");

  if (f == nullptr)
    return 0;
  tex = file::read(f, tex, unitIdx);
  fclose(f);
  return tex;
}

} // end namespace ktx

} // end namespace cg
