#include "graphics/GLImage.h"
#include "TextureTestBase.h"

namespace
{

inline uint32_t
checkerboard(const cg::Color& color, const int cellSize = 16)
{
  const auto boardSize = cellSize * 16;
  cg::ImageBuffer buffer{boardSize, boardSize};
  const cg::Pixel pattern[2] = {cg::Color::black, color};

  for (int p = 0, iy = 0, y = 0; y < 16; ++y, iy ^= 1)
    for (int i = 0; i < cellSize; ++i)
      for (int ix = iy, x = 0; x < 16; ++x, ix ^= 1)
        for (int j = 0; j < cellSize; ++j)
          buffer[p++] = pattern[ix];

  GLuint texture;

  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGB8, boardSize, boardSize);
  glTexSubImage2D(GL_TEXTURE_2D,
    0,
    0,
    0,
    boardSize,
    boardSize,
    GL_RGB,
    GL_UNSIGNED_BYTE,
    buffer.data());
  return texture;
}

}

int TextureTestBase::_textureUseCount;
cg::TextureMap TextureTestBase::_defaultTextures;

static int tsCount;

TextureTestBase::TextureSelector::TextureSelector(GLint textureUnitIndex):
  _id{"_texSelector_" + tsCount++},
  _textureUnitIndex{textureUnitIndex}
{
  // do nothing
}

void
TextureTestBase::TextureSelector::initialize()
{
  auto tit = _defaultTextures.begin();
  setTexture(tit->first, tit->second);
}

void
TextureTestBase::TextureSelector::inspect()
{
  constexpr int maxSize = 32;
  char buffer[maxSize];

  ImGui::PushID(_id.c_str());
  snprintf(buffer, maxSize, "%s", _textureName.c_str());
  ImGui::InputText("Texture", buffer, maxSize, ImGuiInputTextFlags_ReadOnly);
  ImGui::SameLine();
  if (ImGui::Button("..."))
    ImGui::OpenPopup("Popup");
  if (ImGui::BeginPopup("Popup"))
  {
    for (auto& t : _defaultTextures)
      if (ImGui::Selectable(t.first.c_str()))
      {
        setTexture(t.first, t.second);
        break;
      }

    auto& textures = cg::Assets::textures();

    if (!textures.empty())
    {
      ImGui::Separator();
      for (auto tit = textures.begin(); tit != textures.end(); ++tit)
        if (ImGui::Selectable(tit->first.c_str()))
        {
          setTexture(tit->first,
            cg::Assets::loadTexture(tit, _textureUnitIndex));
          break;
        }
    }
    ImGui::EndPopup();
  }
  ImGui::PopID();
}

inline void
TextureTestBase::initializeTextures()
{
  if (_textureUseCount++ == 0)
  {
    _defaultTextures["None"] = 0;
    _defaultTextures["Checkerboard"] = checkerboard(cg::Color::yellow, 32);
    cg::Assets::initialize();
  }
}

void
TextureTestBase::initialize()
{
  initializeTextures();
  GLTest::initialize();
}

inline void
TextureTestBase::clearTextures()
{
  if (--_textureUseCount == 0)
  {
    cg::Assets::clear();
    cg::deleteTextures(_defaultTextures);
  }
}

void
TextureTestBase::terminate()
{
  GLTest::terminate();
  clearTextures();
}
