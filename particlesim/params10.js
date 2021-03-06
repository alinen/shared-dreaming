
class Params
{
   constructor() 
   {
      this.h = 0.6; // particle neighborhood size
      this.r = 0.25; // particle size
      this.dt = 0.1; // time step
      this.rho0 = 0.5; // reference density
      this.k = 0.005; //bulk modulus
      this.mu = 0.05; // viscosity
      this.g = -0.1; // gravity strength
      this.min = vec3.fromValues(-2, -1.25, -2);
      this.max = vec3.fromValues(2, 1.25, 2);
      this.maxacc = 1.5;
      this.kIntersect = 30;
   }  
};

var system = createSystem(10, 25); 

