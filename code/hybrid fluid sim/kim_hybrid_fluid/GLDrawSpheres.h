#ifndef __tp_GLDrawSpheres_h
#define __tp_GLDrawSpheres_h

#define STRINGIFY(A) "#version 400\n"#A

  static const char* vertexShader = STRINGIFY(
  layout(location = 0) in vec2 position;
  layout(location = 1) in vec2 velocity;
  layout(location = 2) in float alphas;
  layout(location = 3) in float indices;

  uniform bool use_color_map = false;
  uniform float radius;
  uniform mat4 mvMatrix;
  uniform mat4 projectionMatrix;
  uniform vec2 viewportSize;
  uniform vec4 color;
  out float vRadius;
  out vec4 gColor;
  out vec2 texCoord;

  void main()
  {
    vRadius = radius;

    vec4 p = mvMatrix * vec4(position.xy, 0, 1);
    float d = 2 * vRadius;
    vec4 v = projectionMatrix * vec4(d, d, p.z, p.w);
    vec2 s = viewportSize * v.xy / v.w;

    gl_PointSize = .25f *(s.x + s.y);
    gColor = use_color_map ? color * min(length(velocity), 1.0f) : color;
    gColor.a = alphas;
    gl_Position = projectionMatrix * p;
  }
  );

  static const char* geometryShader = STRINGIFY(
    layout(points) in;
  layout(triangle_strip, max_vertices = 4) out;

  in float vRadius[];
  in vec4 vColor[];
  uniform mat4 projectionMatrix;
  out vec2 texCoord;
  out vec4 gColor;

  void main()
  {
    vec4 p = gl_in[0].gl_Position;
    float pointSize = gl_in[0].gl_PointSize;
    float r = vRadius[0];

    gl_Position = projectionMatrix * (vec4(-r, -r, 0, 0) + p);
    texCoord = vec2(0, 0);
    gColor = vColor[0];
    gl_PointSize = pointSize;
    EmitVertex();
    gl_Position = projectionMatrix * (vec4(+r, -r, 0, 0) + p);
    texCoord = vec2(1, 0);
    gColor = vColor[0];
    gl_PointSize = pointSize;
    EmitVertex();
    gl_Position = projectionMatrix * (vec4(-r, +r, 0, 0) + p);
    texCoord = vec2(0, 1);
    gColor = vColor[0];
    gl_PointSize = pointSize;
    EmitVertex();
    gl_Position = projectionMatrix * (vec4(+r, +r, 0, 0) + p);
    texCoord = vec2(1, 1);
    gColor = vColor[0];
    gl_PointSize = pointSize;
    EmitVertex();
    EndPrimitive();
  }
  );

  static const char* fragmentShader = STRINGIFY(
    in vec2 texCoord;
  in vec4 gColor;
  out vec4 fragmentColor;
  uniform vec4 ambientColor = vec4(0.1, 0.1, 0.1, 1.0);
  uniform vec4 lightColor = vec4(1);
  uniform vec3 lightDirection = vec3(-0.5773, 0.5773, 0.5773);

  vec4 shade(vec4 color, vec3 N)
  {
    return ambientColor + (lightColor * color) * max(dot(lightDirection, N), 0);
  }

  void main()
  {
    /*vec3 N;

    N.xy = texCoord * 2 - vec2(1);

    float s = dot(N.xy, N.xy);

    if (s > 1)
      discard;
    N.z = sqrt(1 - s);*/
    //fragmentColor = shade(gColor, N);
    fragmentColor = gColor;

  }
  );

#endif // __tp_GLDrawSpheres_h