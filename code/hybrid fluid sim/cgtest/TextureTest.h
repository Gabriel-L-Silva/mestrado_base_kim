#ifndef __TextureTest_h
#define __TextureTest_h

#include "TextureTestBase.h"

class TextureTest: public TextureTestBase
{
public:
  TextureTest(const char* title = "Texture");

protected:
  void setCurrent() override;
  void initialize() override;
  void run() override;
  void gui() override;

private:
  using Base = TextureTestBase;

  cg::vec2f _scale{1.0f};
  TextureSelector _texSelector{1};

}; // TextureTest

#endif // __TextureTest_h
