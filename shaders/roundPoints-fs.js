var fragmentShaderSource_Solid = `
   precision mediump float;
   varying vec4 vColor;
   varying vec4 vPos;

   void main(void) 
   {
      gl_FragColor = vColor;
      vec2 halfSize = 0.5 * vec2(400,500);
      vec2 screenPos = halfSize + vPos.xy/vPos.w * halfSize;
      vec2 pos = gl_FragCoord.xy - screenPos;
      float dist_squared = dot(pos, pos);
      if (dist_squared > 10.0) discard; 
   }
`;
