var vertexShaderSource = `
    attribute vec3 aVertexPosition;
    attribute vec4 aColor;

    uniform mat4 uMVMatrix;
    uniform mat4 uPMatrix;
    varying vec4 vColor;
    varying vec4 vPos;

  void main() {
        gl_Position = uPMatrix * uMVMatrix * vec4(aVertexPosition, 1.0);
        gl_PointSize = 10.0;
        vColor = aColor;
        vPos = gl_Position;
  }
`; 

