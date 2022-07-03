#ifndef __TextureTestBase_h
#define __TextureTestBase_h

#include "Assets.h"
#include "GLTest.h"

class TextureTestBase abstract: public GLTest
{
protected:
  class TextureSelector
  {
  public:
    TextureSelector(GLint textureUnitIndex);

    void initialize();

    auto textureUnitIndex() const
    {
      return _textureUnitIndex;
    }

    auto texture() const
    {
      return _texture;
    }

    auto textureName() const
    {
      return _textureName;
    }

    void inspect();

    void bindTexture()
    {
      glActiveTexture(GL_TEXTURE0 + _textureUnitIndex);
      glBindTexture(GL_TEXTURE_2D, _texture);
    }

  private:
    std::string _id;
    GLint _textureUnitIndex;
    uint32_t _texture{};
    std::string _textureName;

    void setTexture(const std::string& textureName, uint32_t texture)
    {
      _textureName = textureName;
      _texture = texture;
      bindTexture();
    }

  }; // TextureSelector

  using GLTest::GLTest;

  void initialize() override;
  void terminate() override;

private:
  static int _textureUseCount;
  static cg::TextureMap _defaultTextures;

  static void initializeTextures();
  static void clearTextures();

}; // TextureTestBase

#endif // __TextureTestBase_h
