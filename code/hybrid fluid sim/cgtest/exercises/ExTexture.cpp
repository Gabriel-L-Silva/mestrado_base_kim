#include "graphics/GLImage.h"
#include "ExTexture.h"

#define STRINGIFY(A) "#version 400\n"#A

static const char* vs = STRINGIFY(
out vec2 uv;

void main()
{
  const vec4 positions[] = vec4[](
    vec4(-0.5, -0.5, 0.0, 1.0),
    vec4(+0.5, -0.5, 0.0, 1.0),
    vec4(+0.5, +0.5, 0.0, 1.0),
    vec4(-0.5, +0.5, 0.0, 1.0));
  const vec2 textureCoordinates[] = vec2[](
    vec2(0.0, 0.0),
    vec2(1.0, 0.0),
    vec2(1.0, 1.0),
    vec2(0.0, 1.0));

  gl_Position = positions[gl_VertexID];
  uv = textureCoordinates[gl_VertexID];
}
);

static const char* fs = STRINGIFY(
in vec2 uv;
uniform sampler2D tex;
out vec4 fragmentColor;

void main()
{
  fragmentColor = texture(tex, uv);
}
);

inline GLuint
makeTexture(int width, int height, void* data)
{
  GLuint texture;

  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGB8, width, height);
  glTexSubImage2D(GL_TEXTURE_2D,
    0, // level
    0, // xoffset
    0, // yoffset
    width,
    height,
    GL_RGB, // format
    GL_UNSIGNED_BYTE, // type
    data);
  return texture;
}

inline GLuint
checkerboard(const cg::Color& b, const cg::Color w, int cellSize)
{
  const auto boardSize = cellSize * 16;
  cg::ImageBuffer buffer{boardSize, boardSize};
  const cg::Pixel pattern[]{b, w};

  for (int p = 0, iy = 0, y = 0; y < 16; ++y, iy ^= 1)
    for (int i = 0; i < cellSize; ++i)
      for (int ix = iy, x = 0; x < 16; ++x, ix ^= 1)
        for (int j = 0; j < cellSize; ++j)
          buffer[p++] = pattern[ix];
  return makeTexture(boardSize, boardSize, (void*)buffer.data());
}

ExTexture::ExTexture(int width, int height):
  cg::GLWindow("Exercise: Texture", width, height),
  _program{"ExTexture"}
{
  // do nothing
}

void
ExTexture::initialize()
{
  _program.setShaders(vs, fs).use();
  _texture = checkerboard(cg::Color::red, cg::Color::blue, 32);
  glGenVertexArrays(1, &_vao);
  glBindVertexArray(_vao);
}

void
ExTexture::render()
{
  clear(cg::Color::darkGray);
  glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
}

void
ExTexture::terminate()
{
  glDeleteVertexArrays(1, &_vao);
  glDeleteTextures(1, &_texture);
}
