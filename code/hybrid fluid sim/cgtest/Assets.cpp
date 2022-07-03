//[]---------------------------------------------------------------[]
//|                                                                 |
//| Copyright (C) 2018, 2020 Orthrus Group.                         |
//|                                                                 |
//| This software is provided 'as-is', without any express or       |
//| implied warranty. In no event will the authors be held liable   |
//| for any damages arising from the use of this software.          |
//|                                                                 |
//| Permission is granted to anyone to use this software for any    |
//| purpose, including commercial applications, and to alter it and |
//| redistribute it freely, subject to the following restrictions:  |
//|                                                                 |
//| 1. The origin of this software must not be misrepresented; you  |
//| must not claim that you wrote the original software. If you use |
//| this software in a product, an acknowledgment in the product    |
//| documentation would be appreciated but is not required.         |
//|                                                                 |
//| 2. Altered source versions must be plainly marked as such, and  |
//| must not be misrepresented as being the original software.      |
//|                                                                 |
//| 3. This notice may not be removed or altered from any source    |
//| distribution.                                                   |
//|                                                                 |
//[]---------------------------------------------------------------[]
//
// OVERVIEW: Assets.cpp
// ========
// Source file for assets.
//
// Author: Paulo Pagliosa
// Last revision: 31/03/2020

#include "graphics/Application.h"
#include "Assets.h"
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#else
#define NOMINMAX
#include <GL/gl3w.h>
#endif
#include <filesystem>

namespace cg
{ // begin namespace cg

namespace fs = std::filesystem;

void
deleteTextures(TextureMap& textures)
{
  for (auto pair : textures)
    glDeleteTextures(1, &pair.second);
}


/////////////////////////////////////////////////////////////////////
//
// Assets implementation
// ======
static fs::path texPath;
TextureMap Assets::_textures;

void
Assets::initialize()
{
  texPath = Application::assetFilePath("textures/");
  if (fs::is_directory(texPath))
  {
    auto p = fs::directory_iterator(texPath);

    for (auto e = fs::directory_iterator(); p != e; ++p)
      if (fs::is_regular_file(p->status()))
        _textures[p->path().filename().string()] = 0;
  }
}

uint32_t
Assets::loadTexture(const TextureMapIterator& tit, uint32_t unitIdx)
{
  if (tit == _textures.end())
    return 0;

  auto t = tit->second;

  if (t == 0)
  {
    auto filename = (texPath / tit->first).string();

    t = ktx::load(filename.c_str(), 0, unitIdx);
    _textures[tit->first] = t;
  }
  return t;
}

} // end namespace cg
