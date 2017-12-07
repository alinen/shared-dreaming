var PALM_POS = 0;
var CARP = 0;
var MCP = 1;
var PIP = 2;
var DIP = 3;

var squareVertexPositionBuffer;
var mvMatrix = mat4.create();
var pMatrix = mat4.create();
var shaderProgram;
var elapsedTime = 0;
var lastTime = 0;
var numSpheres = 50;
var textureId = 0;
var handTextureId = 1;
var gl;
var system = new Bubbles(50, 100);
var numHandData = 2 * (4*5 + 1); // 2 hands * ( 4 fingerJoints * 5 fingers + palmPos)
var hand = new HandHelper();

function initGL(canvas)
{
   try
   {
      gl = canvas.getContext("experimental-webgl");
      gl.viewportWidth = canvas.width;
      gl.viewportHeight = canvas.height;
      gl.getExtension('OES_texture_float');
   }
   catch (e)
  {
   }
   if (!gl)
   {
      alert("No WebGL for you");
   }
}

function getShader(gl, src, type)
{
   var shader = gl.createShader(type);
   gl.shaderSource(shader, src);
   gl.compileShader(shader);

   if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS))
   {
      alert(gl.getShaderInfoLog(shader));
      return null;
   }

   return shader;
}

function initShader(fragmentShader, vertexShader)
{
   var sp = gl.createProgram();
   gl.attachShader(sp, vertexShader);
   gl.attachShader(sp, fragmentShader);
   gl.linkProgram(sp);

   if (!gl.getProgramParameter(sp, gl.LINK_STATUS))
   {
      alert("Could not initialise shaders");
   }

   sp.vertexPositionAttribute = gl.getAttribLocation(sp, "pos");
   gl.enableVertexAttribArray(sp.vertexPositionAttribute);

   return sp;
}

function initShaders()
{
   var fragmentBackground = getShader(gl, fragmentProgram, gl.FRAGMENT_SHADER);
   var vertexShader = getShader(gl, vertexProgram, gl.VERTEX_SHADER);
   shaderProgram = initShader(fragmentBackground, vertexShader);
}

function initHandRecording()
{
  hand.getFrameData();
}

function initTexture()
{
   textureId = gl.createTexture();

   gl.bindTexture(gl.TEXTURE_2D, textureId);

   gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
   gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
   gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
   gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

   textureSize = system.maxNumSpheres * 3;

   gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, textureSize, 1, 0, gl.RGBA, gl.FLOAT, system.data);
   gl.bindTexture(gl.TEXTURE_2D, null);
}

function initHandTexture()
{
  handTextureId = gl.createTexture();

  gl.bindTexture(gl.TEXTURE_2D, handTextureId);

  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

  handTextureSize = numHandsData * 2;

   gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, handTextureSize, 1, 0, gl.RGBA, gl.FLOAT, hand.data);
   gl.bindTexture(gl.TEXTURE_2D, null);
}

function createGlBuffer(values, itemSize, numItems, type)
{
   var buffer = gl.createBuffer();
   gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
   gl.bufferData(gl.ARRAY_BUFFER, values, type);
   buffer.itemSize = itemSize;
   buffer.numItems = numItems;
   return buffer;
}

function getSquareVertices (vertices)
{
  // GLfloat camera_positon = {0, 0, 2};
  var near_plane = 2;
  var fov = 30;
  var screen = [720, 480];

  var world_height;
  var world_width;

  world_height = 2 * Math.tan(fov);
  world_width = (world_height * screen[0]) / screen[1];

  var p1x = -(world_width/2);
  var p1y = world_height/2;
  var p1z = 0;

  var p2x = world_width/2;
  var p2y = world_height/2;
  var p2z = 0;

  var p3x = world_width/2;
  var p3y = -(world_height/2);
  var p3z = 0;

  var p4x = -(world_width/2);
  var p4y = -(world_height/2);
  var p4z = 0;

  vertices[0] = p1x;
  vertices[1] = p1y;
  vertices[2] = p1z;

  vertices[3] = p2x;
  vertices[4] = p2y;
  vertices[5] = p2z;

  vertices[6] = p3x;
  vertices[7] = p3y;
  vertices[8] = p3z;

  vertices[9] = p4x;
  vertices[10] = p4y;
  vertices[11] = p4z;
}

function initBuffers()
{
   var vertices = [];
   getSquareVertices(vertices);
   mat4.ortho(pMatrix, -1.0, 1.0, -1.0, 1.0, -1.0, 10.0);

   var flatVerticesBuffer = [];
   for (var i = 0; i < vertices.length; i+=3) {
      x = vertices[i+0]
      y = vertices[i+1]
      z = vertices[i+2]
      var tmpVec4 = vec4.create()
      vec4.set(tmpVec4, x, y, z, 1)

      vec4.transformMat4(tmpVec4, tmpVec4, pMatrix)
      flatVerticesBuffer.push(tmpVec4[0], tmpVec4[1], tmpVec4[2]);
   }

   squareVertexPositionBuffer = createGlBuffer(new Float32Array(flatVerticesBuffer), 3, 4, gl.STATIC_DRAW);
}

function drawScene()
{
   gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
   gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

   gl.useProgram(shaderProgram);
   gl.bindBuffer(gl.ARRAY_BUFFER, squareVertexPositionBuffer);
   gl.vertexAttribPointer(shaderProgram.vertexPositionAttribute, squareVertexPositionBuffer.itemSize, gl.FLOAT, false, 0, 0);

   var uelapsedTime = gl.getUniformLocation(shaderProgram, "elapsedTime");
   gl.uniform1f(uelapsedTime, elapsedTime);

   var uTextureSize = gl.getUniformLocation(shaderProgram, "size_of_texture");
   gl.uniform1f(uTextureSize, textureSize);

      var uHandTextureSize = gl.getUniformLocation(shaderProgram, "size_of_hand_texture");
   gl.uniform1f(uHandTextureSize, handTextureSize);

   var uNumSpheres = gl.getUniformLocation(shaderProgram, "num_of_spheres");
   gl.uniform1f(uNumSpheres, numSpheres);

   gl.activeTexture(gl.TEXTURE0);
   gl.bindTexture(gl.TEXTURE_2D, textureId);
   gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, textureSize, 1, 0, gl.RGBA, gl.FLOAT, system.data);

    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, handTextureId);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, handTextureSize, 1, 0, gl.RGBA, gl.FLOAT, hand.data);

   var uSphereInfo = gl.getUniformLocation(shaderProgram, "sphere_info");
   gl.uniform1i(uSphereInfo, 0);

   var uHandInfo = gl.getUniformLocation(shaderProgram, "hand_info");
   gl.uniform1i(uHandInfo, 1);

   gl.drawArrays(gl.TRIANGLE_FAN, 0, squareVertexPositionBuffer.numItems);
}

function tick()
{
   var newTime = new Date().getTime();
   var dt = (newTime - lastTime)*0.001;
   elapsedTime += dt * 4.0;

   system.update(dt);
   hand.update(dt);
   requestAnimFrame(tick);
   drawScene();
   lastTime = newTime
}

function webGLStart()
{
  var canvas = document.getElementById("canvas");
  initHandRecording();
  initGL(canvas);
  initShaders();
  initBuffers();
  initTexture();
  initHandTexture();

  gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);

  lastTime =  new Date().getTime();
  gl.clearColor(0.1, 0.1, 0.1, 1.0);
  tick();
}


