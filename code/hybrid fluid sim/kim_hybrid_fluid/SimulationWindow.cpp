#include "SimulationWindow.h"
#include "graphics/Color.h"

//using namespace jet;

SimulationWindow::SimulationWindow(int width, int height) :
  cg::GLWindow("Dam Breaking Simulation Window", width, height),
  _program{ "GLRenderer" }
{
  // TODO
}

void
SimulationWindow::initialize()
{
  // TODO: seletor de simulação
  //  \_ mostrar owner da simulação + nome da simulação
  // TODO: reset button

  //using namespace jet;

  //_program.setShaders(vertexShader, fragmentShader);
  //_program.setShader(GL_GEOMETRY_SHADER, geometryShader);
  //_program.use();

  //glGenVertexArrays(1, &_vao);
  //glBindVertexArray(_vao);

  //_solver = FlipSolver2::builder()
  //  .withResolution({ 100, 100 })
  //  .withDomainSizeX(1.0)
  //  .makeShared();

  //// Build emitter
  //_surface = jet::Box2::builder()
  //  .withLowerCorner({ 0.0, 0.0 })
  //  .withUpperCorner({ 0.2, 0.8 })
  //  .makeShared();

  //_emitter = VolumeParticleEmitter2::builder()
  //  .withSurface(_surface)
  //  .withSpacing(_particleSpacing)
  //  .withIsOneShot(true)
  //  .makeShared();

  //_solver->setParticleEmitter(_emitter);
  //_solver->update(_frame++);
  //auto data = _solver->particleSystemData();
  //auto n = data->numberOfParticles();
  //_radius = data->radius();
  //_program.setUniform("radius", (float) _radius);
  //// TODO: allow window resize
  //_program.setUniformVec2("viewportSize", cg::vec2f{ 1270, 720 });
  //Logging::setLevel(LoggingLevel::Warn);
  //if (n > 0)
  //{
  //  _positions = new cg::GLBuffer<Vector2D>(n);
  //  _velocities = new cg::GLBuffer<Vector2D>(n);

  //  _positions->bind();
  //  glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 0, 0);
  //  glEnableVertexAttribArray(0);
  //  _velocities->bind();
  //  glVertexAttribPointer(1, 2, GL_DOUBLE, GL_FALSE, 0, 0);
  //  glEnableVertexAttribArray(1);
  //}
}

void
SimulationWindow::gui()
{
  ImGui::Begin("Simulation Controller");
  if (ImGui::Button("Pause"))
    _paused = !_paused;
  if (ImGui::Button("Reset"))
    resetSimulation();
  if (ImGui::Button("Enable/Disable Velocity Color Map"))
    _enableColorMap = !_enableColorMap;
  ImGui::SameLine();
  ImGui::Text("false\0true" + (_enableColorMap * 6));
  if (ImGui::CollapsingHeader("Particle Emmiter"))
  {
    surfaceOptions();
    ImGui::Text("Particle Spacing");
    ImGui::DragFloat("", &_particleSpacing, 0.001f, 0.001, 0.2);
    ImGui::ColorEdit3("Particles Color", (float*)&_particleColor);
    
  }
  ImGui::End();
}

void
SimulationWindow::render()
{
  clear(cg::Color::gray);

  /*if (!_paused)
    _solver->update(_frame++);

  _program.use();
  glBindVertexArray(_vao);

  auto data = _solver->particleSystemData();
  ConstArrayAccessor1<Vector2D> positions = data->positions();
  _positions->bind();
  _positions->setData(positions.data());

  if (_enableColorMap)
  {
    ConstArrayAccessor1<Vector2D> velocities = data->velocities();
    _velocities->bind();
    _velocities->setData(velocities.data());
  }

  auto projectionMatrix = cg::mat4f::perspective(
    60,
    ((float)width()) / height(),
    0.001f,
    100
  );

  auto mvMatrix = cg::mat4f::TRS(
    cg::vec3f{ 0.5f, 0.5f, 2.0f },
    cg::vec3f::null(),
    cg::vec3f{ 1.0f }
  );
  mvMatrix.invert();

  _program.setUniform("use_color_map", _enableColorMap);
  _program.setUniformMat4("projectionMatrix", projectionMatrix);
  _program.setUniformMat4("mvMatrix", mvMatrix);
  _program.setUniformVec4("color", _particleColor);

  glDrawArrays(GL_POINTS, 0, data->numberOfParticles());*/
}

void
SimulationWindow::terminate()
{
}

void
SimulationWindow::resetSimulation()
{
  /*_solver = FlipSolver2::builder()
    .withResolution({ 100, 100 })
    .withDomainSizeX(1.0)
    .makeShared();

  _emitter = VolumeParticleEmitter2::builder()
    .withSurface(_surface)
    .withSpacing(_particleSpacing)
    .withIsOneShot(true)
    .makeShared();

  _solver->setParticleEmitter(_emitter);
  _frame = Frame();
  _solver->update(_frame++);

  auto data = _solver->particleSystemData();
  auto n = data->numberOfParticles();
  _radius = data->radius();
  if (n > 0 && _positions->size() != n)
  {
    _positions->bind();
    _positions->resize(n);
  }*/
}

void
SimulationWindow::surfaceOptions()
{
  /*const char* items[] = {"Box", "Sphere"};
  static const char* current_item = items[0];

  if (ImGui::BeginCombo("##SurfaceType", current_item, ImGuiComboFlags_NoArrowButton))
  {
    for (int n = 0; n < IM_ARRAYSIZE(items); ++n)
    {
      bool selected = (current_item == items[n]);
      if (ImGui::Selectable(items[n], selected))
      {
        if (current_item != items[n])
        {
          if (n == 0)
            _surface = Box2::builder()
              .withLowerCorner({ 0.0, 0.0 })
              .withUpperCorner({ 0.2, 0.8 })
              .makeShared();
          else
            _surface = Sphere2::builder()
              .withCenter({ 0.5, 0.5 })
              .withRadius(0.2)
              .makeShared();
        }
        current_item = items[n];
      }
      if(selected)
        ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }*/
}
