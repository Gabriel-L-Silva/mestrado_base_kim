#ifndef __MultiTextureTest_h
#define __MultiTextureTest_h

#include "TextureTestBase.h"

class MultiTextureTest: public TextureTestBase
{
public:
  MultiTextureTest(const char* title = "Multiple Textures");

protected:
  void setCurrent() override;
  void initialize() override;
  void run() override;
  void gui() override;

private:
  using Base = TextureTestBase;

  float _time{};
  cg::vec2f _scale{1.0f};
  TextureSelector _texSelector1{1};
  TextureSelector _texSelector2{2};
  bool _rotationFlag{true};

}; // MultiTextureTest

#endif // __MultiTextureTest_h
