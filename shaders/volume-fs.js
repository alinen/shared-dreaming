var fragmentProgram =`
// relection
#ifdef GL_ES
  precision mediump float;
#endif

varying vec4 v_pos;
uniform float elapsedTime;
uniform sampler2D sphere_info; // width = num_of_spheres * 3 * 4; height = 1
uniform float num_of_spheres;
uniform float size_of_texture;
uniform float threshold;
uniform vec3 left_bs_center;
uniform float left_bs_radius;
uniform vec3 right_bs_center;
uniform float right_bs_radius;

void planeIntersection(in vec3 ray_start, in vec3 ray_dir, in vec3 box_normal, in vec3 box_p1,
  out float plane_intersect)
  {
    float term_1 = dot(ray_start, box_normal);
    float term_2 = dot(box_p1, box_normal);
    float term_3 = dot(ray_dir, box_normal);

    float t = (term_2 - term_1) / term_3;

    plane_intersect = t;
  }

void boxIntersection(in vec3 ray_start, in vec3 ray_dir, out vec3 intersection_point, out float min_t)
{
  vec3 box_normals[5];
  vec3 box_points[5];
  float t;

  box_normals[0] = vec3(0.0, -1.0, 0.0); // top
  box_normals[1] = vec3(0.0, 0.0, -1.0); // back
  box_normals[2] = vec3(1.0, 0.0, 0.0); // left
  box_normals[3] = vec3(-1.0, 0.0, 0.0); // right
  box_normals[4] = vec3(0.0, 1.0, 0.0); // bottom

  box_points[0] = vec3(0.0, 3.0, 0.0); // top
  box_points[1] = vec3(0.0, 0.0, -4.0); // back
  box_points[2] = vec3(3.0, 0.0, 0.0); // left
  box_points[3] = vec3(-3.0, 0.0, 0.0); // right
  box_points[4] = vec3(0.0, -3.0, 0.0); // bottom

  min_t = 100.0;

  for (int i = 0; i < 5; ++i) {
    planeIntersection(ray_start, ray_dir, box_normals[i], box_points[i], t);
    if (t > 0.0)
      min_t = min(min_t, t);
  }

  intersection_point = ray_start + min_t * ray_dir;
}

void computeColor(in vec3 ray_start, in vec3 ray_dir, out vec4 color)
{
  vec3 intersection_point;
  float time;
  boxIntersection(ray_start, ray_dir, intersection_point, time);

  float fraction_of_screen = intersection_point.y + 2.0 / 4.0;

  vec4 top_color = vec4(0.645, 0.616, 0.57, 1.0);
  vec4 bottom_color = vec4(0.15, 0.56, 0.7, 1.0);

  color = top_color * fraction_of_screen + bottom_color * (1.0 - fraction_of_screen);
}

float density1(float a, float b, float rSqr)
{
  return a * exp(-b * rSqr);
}

void getVel(float startIndex, out vec4 vel)
{
  float tex_coord_2 = (startIndex + 1.0)/size_of_texture + 1.0/(2.0 * size_of_texture);
  float tex_coord_y = 0.5;
  vel = texture2D(sphere_info, vec2(tex_coord_2, tex_coord_y));
}

void getColor(float startIndex, out vec4 color)
{
  float tex_coord_3 = (startIndex + 2.0)/size_of_texture + 1.0/(2.0 * size_of_texture);
  float tex_coord_y = 0.5;
  color = texture2D(sphere_info, vec2(tex_coord_3, tex_coord_y));
}

void getPosition(float startIndex, out vec3 center)
{
  float tex_coord_1 = (startIndex + 0.0)/size_of_texture + 1.0/(2.0 * size_of_texture);
  float tex_coord_y = 0.5;

  vec4 pos_rad = texture2D(sphere_info, vec2(tex_coord_1, tex_coord_y));
  float R = 0.05;
  float r = sin(elapsedTime+startIndex*2.0)+1.1;
  float x = r*R*cos(elapsedTime+startIndex); 
  float y = r*R*sin(elapsedTime+startIndex); 
  center = pos_rad.xyz + vec3(x,y,y); 
}

void sphereIntersection(in vec3 ray_start, in vec3 ray_dir, in vec3 center, in float radius, out float t)
{
  vec3 sphere_dir = center - ray_start;            // intersection test
  float sphere_len = length(center - ray_start);            // intersection test
  float projection = dot(sphere_dir, ray_dir); // intersection test
  vec3 dir_perpendicular = sphere_dir - (ray_dir * projection); // intersection test
  float len_dir_perpend = length(dir_perpendicular);       // intersection test
  if (len_dir_perpend > radius) 
  {
    t = -1.0;
    return;
  }

  float intersection_dist = sqrt(radius * radius - len_dir_perpend * len_dir_perpend);
  if (sphere_len > radius) 
  {
    t = projection - intersection_dist;
  } 
  else 
  {
    t = projection + intersection_dist;
  }
}

void blobIntersection(in vec3 ray_start, in vec3 ray_dir, in vec3 center, in float radius, 
    out float t, out float density)
{
  density = 0.0;

  vec3 sphere_dir = center - ray_start;            // intersection test
  float sphere_len = length(center - ray_start);            // intersection test
  float projection = dot(sphere_dir, ray_dir); // intersection test
  if (projection < 0.0)
  {
    t = -1.0;
    return;
  }

  vec3 dir_perpendicular = sphere_dir - (ray_dir * projection); // intersection test
  float len_dir_perpend = length(dir_perpendicular);       // intersection test
  if (len_dir_perpend > radius) 
  {
    t = -1.0;
    return;
  }

  float intersection_dist = sqrt(radius * radius - len_dir_perpend * len_dir_perpend);
  if (sphere_len > radius) 
  {
    t = projection - intersection_dist;
  } 
  else // inside the volume => ignore for now
  {
    t = -1.0;
    return;
  }

  float a = 0.05;
  float b = 1000.0;
  density = 0.0;
  for (float tt = 0.0; tt < 1.0; tt += 0.1) // sum effect of half sphere
  {
    float rSqr = intersection_dist * intersection_dist * tt*tt + len_dir_perpend * len_dir_perpend;
    density += density1(a, b, rSqr); 
  }
  density = density * 2.0; // add the remaining influence
}


void checkSpheres(in vec3 ray_start, in vec3 ray_dir, out float t, out vec4 color)
{
  t = -1.0;

  // Idea: check all intersecting spheres to get density
  // Color pixel based on density
  //
  float closest = 100.0;
  float density = 0.0;
  for (float i = 0.0; i < 500.0; i+=1.0) // need to hardcode loop
  { 
    float startIndex = i * 3.0;
    float radius = 0.05; // ASN TODO: why doesn't this work? pos_rad.w;

    vec3 center;
    getPosition(startIndex, center);

    float hitTime = -1.0;
    float hitDensity = 0.0;
    blobIntersection(ray_start, ray_dir, center, radius, hitTime, hitDensity);
    if (hitTime > 0.0) 
    {
      density += hitDensity;
      if (hitTime < closest) closest = hitTime;
    } 
  }
  color = vec4(1.0, 0.5, 1.0, density);
  if (closest < 100.0) t = closest;
}

void main ()
{
  vec3 camera_pos = vec3(0.0, 0.0, 0.0);
  vec3 sphere_center = vec3(0.0, 0.0, -2.0);
  float radius = 0.10;
  float maxDisplacement = 0.2;

  vec3 view_dir = v_pos.xyz - camera_pos;
  vec3 normalized_view_dir = normalize(view_dir);

  // rotate camera by 45
  float angle1 = 0.0; //elapsedTime * 0.1;
  float angle2 = 0.3;
  // first 3 vars = first column
  mat3 m = mat3(
    cos(angle1), 0, sin(angle1), // first column
    0, 1, 0, // second column
    -sin(angle1), 0, cos(angle1)  // third column
  );
  mat3 m2 = mat3(
    1, 0, 0, // first column
    0, cos(angle2), sin(angle2), // second column
    0, -sin(angle2), cos(angle2)  // third column
  );
  normalized_view_dir = m2 * m * normalized_view_dir;
  camera_pos = m2 * m * (camera_pos - sphere_center) + sphere_center;

  float t = -1.0;
  vec4 hit_color;
  checkSpheres(camera_pos, normalized_view_dir, t, hit_color);

  vec4 color = vec4(1.0,0.0,0.0,1.0);
  computeColor(camera_pos, normalized_view_dir, color);
  if (t > 0.0) 
  {
    gl_FragColor.rgb =  hit_color.rgb * (hit_color.a) + color.rgb * (1.0 - hit_color.a);
    gl_FragColor.a = 1.0;
  }
  else
  {
    gl_FragColor = color;
  }

}
`;
