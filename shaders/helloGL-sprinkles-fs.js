var fragmentProgram = `
   precision mediump float;
   varying float v_depth; 
   void main(void) 
   {
      gl_FragColor = vec4(0.01, 0.01, 0.01, 0.9);
   }
`;

