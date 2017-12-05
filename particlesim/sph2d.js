class Params
{
   constructor() 
   {
      this.h = 16; // particle neighborhood size
      this.r = 4; // particle size
      this.dt = 0.1; // time step
      this.rho0 = 1000; // reference density
      this.k = 800; //bulk modulus
      this.mu = 9000.1; // viscosity
      this.g = 2; // gravity strength
      this.width = 10;
      this.height = 10;
   }  
};

class SPH2D
{
   constructor(numParticles, maxNumParticles, params)
   {
      this.params = params;
      this.numSpheres = numParticles;
      this.maxNumSpheres = maxNumParticles;
      this.mass = 0; // mass
      this.data = null; //positions,vels,rgbs needs to go in data to get rendered;
      this.velocities = [] //vector;
      this.vh = [] //vector;
      this.accelerations = []; //vector
      this.rho = [];  // floats
      this.maxdensity = 0;
      this.paused = true;
      this.intersection = []; // booleans
      this.obstacles = [];
      this.sh = new SphereHelper();

      this.setupSpheres();

      // ASN: the following is holderover from p5 implementation (get rid of it)?
      //positions[0].y = 241;
      //velocities[0].y = -20;
      //damp_reflect(0, new Obstacle(0, this.height, this.width, this.height));
   }

   setupSpheres()
   {
      this.data = new Float32Array(this.maxNumSpheres * 3 * 4)
      for (var i = 0; i < this.data.length; i++)
      {
         this.data[i] = 0.0;
      }
       
      for (var i = 0; i < this.numSpheres; i++)
      {
         this.velocities.push(0,0,0);
         this.vh.push(0,0,0);
         this.accelerations.push(0,0,0);
         this.rho.push(0);
         this.intersection.push(false);
      }
   }

   init()
   {
       this.sh.fromData(0, this.data);
       this.sh.pos[0] = 0.0;
       this.sh.pos[1] = 0.0;
       this.sh.pos[2] = -2.0;
       this.sh.toData(0, this.data);

       /*
      var numcols = 30;
      var margin = this.params.r * 2.2;
      for (var i = 0; i < this.numSpheres; i++)
      {
         var cellj = i % numcols;
         var celli = floor(i/numcols);
         this.positions[i] = createVector(this.params.r+celli*margin, this.params.r+cellj*margin);
         this.positions[i].x += fluidOffset.x;
         this.positions[i].y += fluidOffset.y;
         this.velocities[i] = createVector(0,0);
         this.accelerations[i] = createVector(0,0);
         this.vh[i] = createVector(0,0);
      }
      
      normalize_mass();
      compute_accel();
      leapfrog_start();
      */

   }

    update(dt)  // ASN TODO: Fixed framerate or application framerate?
    {
        if (!this.paused)
        {
            computeAccel();
            leapfrogStep();
        }
    }

   computeDensity()
   {
      maxdensity = 0;
      var h = this.params.h;
      var h2 = h*h;
      var h8 = ( h2*h2 )*( h2*h2 );
      var C = 4 * this.mass / 3.14 / h8;
      for (var i = 0; i < this.numSpheres; i++) this.rho[i] = 0;
      for (var i = 0; i < this.numSpheres; ++i) 
      {
          this.rho[i] += 4 * this.mass / 3.14 / h2;
          for (var j = i+1; j < this.numSpheres; ++j) 
          {
             var dx = this.positions[i].x - this.positions[j].x;
             var dy = this.positions[i].y - this.positions[j].y;
             var r2 = dx*dx + dy*dy;
             var z = h2-r2;
             if (z > 0) 
             {
                var rho_ij = C*z*z*z;
                this.rho[i] += rho_ij;
                this.rho[j] += rho_ij;
             }
          }
      }  
       
      for (var i = 0; i < this.numSpheres; i++) this.maxdensity = max(this.maxdensity, this.rho[i]);
   }

    normalizeMass()
    {
       mass = 1;
       compute_density();
       var rho0 = this.params.rho0;
       var rho2s = 0;
       var rhos = 0;
       for (var i = 0; i < this.numSpheres; ++i) 
       {
          rho2s += (this.rho[i])*(this.rho[i]);
          rhos += this.rho[i];
       }
       this.mass *= ( rho0*rhos / rho2s );
    }

    computeAccel()
    {
      // Unpack basic parameters
      var h = this.params.h;
      var rho0 = this.params.rho0;
      var k = this.params.k;
      var mu = this.params.mu;
      var mass = this.mass;
      var h2 = h*h;

      // Compute density and color
      this.computeDensity();

      // Start with gravity and surface forces
      for (var i = 0; i < this.numSpheres; ++i) 
      {
        this.accelerations[i].x = 0;
        this.accelerations[i].y = this.params.g;
        for (var j = 0; j < obstacles.length; j++)
        {
          var b = capsuleIntersect(this.positions[i], this.velocities[i], obstacles[j]);
          if (b === true)
          {
            // compute repellent force
            var force = capsule_force(this.positions[i], this.velocities[i], obstacles[j]);
            this.accelerations[i].x += force.x;
            this.accelerations[i].y += force.y;        
          }
        }    
      }
      
      // Constants for interaction term
      var C0 = mass / 3.14 / ( (h2)*(h2) );
      var Cp = 15*k;
      var Cv = -40*mu;
      
      // Now compute interaction forces
      for (var i = 0; i < this.numSpheres; ++i) 
      {
        var rhoi = this.rho[i];
        for (var j = i+1; j < this.numSpheres; ++j) 
        {
          var dx = this.positions[i].x-this.positions[j].x;
          var dy = this.positions[i].y-this.positions[j].y;
          var r2 = dx*dx + dy*dy;      
          if (r2 < h2) 
          {
            var rhoj = this.rho[j];
            var q = sqrt(r2)/h;
            var u = 1-q;
            var w0 = C0 * u/rhoi/rhoj;
            var wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
            var wv = w0 * Cv;
            var dvx = this.velocities[i].x - this.velocities[j].x;
            var dvy = this.velocities[i].y - this.velocities[j].y;
            
            this.accelerations[i].x += (wp*dx + wv*dvx);
            this.accelerations[i].y += (wp*dy + wv*dvy);
            
            this.accelerations[j].x -= (wp*dx + wv*dvx);
            this.accelerations[j].y -= (wp*dy + wv*dvy);
          }
        }
      }
    }

    leapfrogStep()
    {
      var maxacc = 100;
      for (var i = 0; i < this.numSpheres; ++i)
      {
        var mag = sqrt(this.accelerations[i].x*this.accelerations[i].x + 
            this.accelerations[i].y*this.accelerations[i].y);
          
        if (mag > maxacc) 
        {
          this.accelerations[i].x = this.accelerations[i].x/mag * maxacc;
          this.accelerations[i].y = this.accelerations[i].y/mag * maxacc;
        }
        this.vh[i].x += this.accelerations[i].x * this.params.dt;
        this.vh[i].y += this.accelerations[i].y * this.params.dt;
      }
      for (var i = 0; i < this.numSpheres; ++i)
      {
        this.velocities[i].x = this.vh[i].x + this.accelerations[i].x * this.params.dt / 2;
        this.velocities[i].y = this.vh[i].y + this.accelerations[i].y * this.params.dt / 2;
      }
      for (var i = 0; i < this.numSpheres; ++i)
      {
        this.positions[i].x += this.vh[i].x * this.params.dt;
        this.positions[i].y += this.vh[i].y * this.params.dt;
      }
      reflectAtBoundaries();
    }

    leapfrogStart()
    {
      for (var i = 0; i < this.numSpheres; ++i)
      {
        this.vh[i].x = this.velocities[i].x + this.accelerations[i].x * this.params.dt / 2;
        this.vh[i].y = this.velocities[i].y + this.accelerations[i].y * this.params.dt / 2;
      }
      for (var i = 0; i < this.numSpheres; ++i)
      {
        this.velocities[i].x += this.accelerations[i].x * this.params.dt;
        this.velocities[i].y += this.accelerations[i].y * this.params.dt;
      }
      for (var i = 0; i < this.numSpheres; ++i)
      {
        this.positions[i].x += this.params.dt * this.vh[i].x;
        this.positions[i].y += this.params.dt * this.vh[i].y;
      }
      reflectAtBoundaries();
    }

    capsuleForce(p, v, obs)
    {
        /* TODO
      var ba = createVector(obs.b.x-obs.a.x, obs.b.y-obs.a.y);
      var pa = createVector(p.x-obs.a.x, p.y-obs.a.y);
      var pb = createVector(p.x-obs.b.x, p.y-obs.b.y);
      
      var len = ba.mag();
      var dir = ba;
      dir.normalize();    
      
      var combinedRadius = obs.r + this.params.r;
        
      var dirT = pa.dot(dir);
      if (dirT > 0 && dirT < len) // might intercept middle part
      {
        var perpx = pa.x - dir.x * dirT; 
        var perpy = pa.y - dir.y * dirT; 
        var perpSq = perpx*perpx + perpy*perpy;
        if (perpSq < combinedRadius*combinedRadius) 
        {
          var scale = (combinedRadius*combinedRadius - perpSq);
          var forcex = perpx * scale;
          var forcey = perpy * scale;
          return createVector(forcex, forcey);
        }
      }
      else
      {
        if (pa.magSq() < combinedRadius*combinedRadius) 
        {
          var scale = combinedRadius*combinedRadius - pa.magSq();
          var forcex = pa.x * scale;
          var forcey = pa.y * scale;
          return createVector(forcex, forcey);
        }
        
        if (pb.magSq() < combinedRadius*combinedRadius)
        {
          var scale = combinedRadius*combinedRadius - pb.magSq();
          var forcex = pb.x * scale;
          var forcey = pb.y * scale;
          return createVector(forcex, forcey);
        }
      }
      return createVector(0,0);
      */
    }

    capsuleIntersect(p, v, obs)
    {
        /*
      var ba = createVector(obs.b.x-obs.a.x, obs.b.y-obs.a.y);
      var pa = createVector(p.x-obs.a.x, p.y-obs.a.y);
      var pb = createVector(p.x-obs.b.x, p.y-obs.b.y);
      
      var len = ba.mag();
      var dir = ba;
      dir.normalize();    
      
      var combinedRadius = obs.r + this.params.r;
      var dirT = pa.dot(dir);
      if (dirT > 0 && dirT < len) // might intercept middle part
      {
        var perpx = pa.x - dir.x * dirT; 
        var perpy = pa.y - dir.y * dirT; 
        var perpSq = perpx*perpx + perpy*perpy;
        if (perpSq <combinedRadius*combinedRadius) return true;
      }
      else
      {
        if (pa.magSq() < combinedRadius*combinedRadius) return true;
        if (pb.magSq() < combinedRadius*combinedRadius) return true;
      }*/
      return false;
    }

    // return time wihen v intercepts our obstacle
    intersection(p, v, obs)
    {
        /*
      var v_perp = createVector(-v.y, v.x); //<>//
      var ba = createVector(obs.b.x-obs.a.x, obs.b.y-obs.a.y);
      var pa = createVector(p.x-obs.a.x, p.y-obs.a.y);
      
      var result = new HitObject();
      var test2_numerator = p5.Vector.dot(pa, v_perp);
      var test2_denominator = p5.Vector.dot(ba, v_perp);
      if (abs(test2_denominator) < 0.0001)
      {
        result.success = false;
        return result;
      }
      
      var t_segment = test2_numerator / test2_denominator;
      if (t_segment < 0 || t_segment > 1.0) 
      {
        result.success = false;
        return result;    
      }
      
      var test1_numerator = ba.x*pa.y + ba.y*pa.x;
      var test1_denominator = p5.Vector.dot(ba, v_perp);
      if (abs(test1_denominator) < 0.0001)
      {
        result.success = false;
        return result;
      }
      var t_ray = test1_numerator / test1_denominator;
      result.success = t_ray > 0;
      result.t = t_ray;
      
      //p5.Vector test = new p5.Vector(p.x + t_ray*v.x, p.y + t_ray*v.y, 0);
      //p5.Vector test2 = new p5.Vector(obs.a.x + t_segment*ba.x, obs.a.y + t_segment*ba.y, 0);
      return result;
      */
    }

    dampReflect(which, obs)
    {
        /*
      // Coefficient of resitiution
      var DAMP = 0.75;
      
      // bounce back along normal direction
      var obsDir = p5.Vector.sub(obs.b, obs.a);   //<>//
      var obsN = createVector(-obsDir.y, obsDir.x);
      obsDir.normalize();
      obsN.normalize();
      
      var velN = p5.Vector.dot(obsN, this.velocities[which]);
      var velT = p5.Vector.dot(obsDir, this.velocities[which]);
      var relfectVelx = -velN * DAMP * (obsN.x) + velT * (obsDir.x);
      var relfectVely = -velN * DAMP * (obsN.y) + velT * (obsDir.y); //<>//
      
      var vhN = p5.Vector.dot(obsN, this.vh[which]);
      var vhT = p5.Vector.dot(obsDir, this.vh[which]);
      var relfectVhx = -vhN * DAMP * obsN.x + vhT * obsDir.x;
      var relfectVhy = -vhN * DAMP * obsN.y + vhT * obsDir.y;
      
      this.velocities[which].x = relfectVelx;
      this.velocities[which].y = relfectVely;
      this.vh[which].x = relfectVhx;
      this.vh[which].y = relfectVhy;
      
      // Damp the velocities
      this.velocities[which].x *= DAMP; 
      this.vh[which].x *= DAMP;
      
      // Reflect the position and velocity
      var diff = p5.Vector.sub(this.positions[which], obs.a);
      var flipDist = diff.dot(obsN); // (p-a).dot(obsN) 
      this.positions[which].x -= 2*flipDist*obsN.x;
      this.positions[which].y -= 2*flipDist*obsN.y;    
      */
    }

    reflectAtBoundaries()
    {
        /*
      // Boundaries of the computational domain

      var XMIN = 0.0;
      var XMAX = this.width;
      var YMIN = 0.0;
      var YMAX = this.height;
      for (var i = 0; i < this.numSpheres; ++i) 
      {
        if (this.positions[i].x < XMIN) damp_reflect(i, new Obstacle(0,0, 0, this.height));
        if (this.positions[i].x > XMAX) damp_reflect(i, new Obstacle(this.width, 0, this.width, this.height));
        if (this.positions[i].y > YMAX) damp_reflect(i, new Obstacle(0, this.height, this.width, this.height));
      }
      */
    }

    clearObstacles()
    {
        obstacles = [];
    }

    pushObstacle(start, end)
    {
        // ASN TODO: Fix for webGL/new world size
        /*
        newObstacle.b = createVector(mouseX, mouseY);
        var dx = newObstacle.a.x - newObstacle.b.x;
        var dy = newObstacle.a.y - newObstacle.b.y;
        var len = dx*dx + dy*dy;
        if (len > 100)
        {
            obstacles.push(newObstacle);
        }*/
    }

    popObstacle()
    {
        if (this.obstacles.length > 0) 
        {
            obstacles.splice(obstacles.length - 1, 1);
        }
    }

};

/* how to draw
  for (var i = 0; i < this.numSpheres; i++)
  {
    var p = this.positions[i];
    var c = this.rho[i]/this.maxdensity;
    if (this.intersection[i])
    {
      fill(c*255, 255, 0);
    }
    else
    {
      fill(c*255, 0, 0);
    }
    ellipse(p.x, p.y, this.params.r*2, this.params.r*2);
  }
  
  stroke(255);
  for (var i = 0; i < obstacles.length; i++)
  {
    line(obstacles[i].a.x, obstacles[i].a.y, obstacles[i].b.x, obstacles[i].b.y);
  }
 */ 

function createSystem(numSpheres, maxSpheres)
{
    return new SPH2D(numSpheres, maxSpheres);
}

