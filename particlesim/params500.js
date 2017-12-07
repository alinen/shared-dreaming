
class Params
{
   constructor() 
   {
      this.h = 0.15; // particle neighborhood size
      this.r = 0.05; // particle size
      this.dt = 0.05; // time step
      this.rho0 = 1.5; // reference density
      this.k = 0.005; //bulk modulus
      this.mu = 0.001; // viscosity
      this.g = -0.01; // gravity strength
      this.min = vec3.fromValues(-2, -1.25, -2);
      this.max = vec3.fromValues(2, 1.25, 2);
      this.maxacc = 0.5;
      this.kIntersect = 10;
   }  
};

var system = createSystem(498, 500); 

