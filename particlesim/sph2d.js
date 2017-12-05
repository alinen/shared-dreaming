class Params
{
   constructor() 
   {
      this.h = 0.5; // particle neighborhood size
      this.r = 0.15; // particle size
      this.dt = 0.05; // time step
      this.rho0 = 1000; // reference density
      this.k = 800; //bulk modulus
      this.mu = 9000.1; // viscosity
      this.g = -0.1; // gravity strength
      this.minx = -2;
      this.maxx = 2;
      this.miny = -2;
      this.maxy = 2;
      this.maxacc = 100;
   }  
};

class Obstacle
{
    constructor(xstart, ystart, xend, yend)
    {
        this.a = vec3.fromValues(xstart, ystart, 0.0);
        this.b = vec3.fromValues(xend, yend, 0.0);
        console.log("OBS "+this.a);
        console.log("OBS "+this.b);
    }
}

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
      this.paused = false;
      this.intersection = []; // booleans
      this.obstacles = [];
      this.sh = new SphereHelper();

      this.setupSpheres();
      this.init();

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
         this.velocities.push(vec3.create());
         this.vh.push(vec3.create());
         this.accelerations.push(vec3.create());
         this.rho.push(0);
         this.intersection.push(false);
      }
   }

   init()
   {

      var numcols = 5;
      var margin = this.params.r * 2.2;
      for (var i = 0; i < this.numSpheres; i++)
      {
         var cellj = i % numcols;
         var celli = Math.floor(i/numcols);

         this.sh.fromData(i, this.data);
         this.sh.radius = this.params.r; 
         this.sh.pos[0] = -2.0 + this.params.r+cellj*margin;
         this.sh.pos[1] = this.params.r+celli*margin;
         this.sh.pos[2] = -2.0;
         this.sh.toData(i, this.data);
         
         vec3.set(this.velocities[i], 0,0,0);
         vec3.set(this.accelerations[i], 0,0,0);
         vec3.set(this.vh[i], 0,0,0);
      }
      
      this.normalizeMass();
      this.computeAccel();
      this.leapfrogStart();
   }

    update(dt)  // ASN TODO: Fixed framerate or application framerate?
    {
        if (!this.paused)
        {
            this.computeAccel();
            this.leapfrogStep();
        }
    }

   computeDensity()
   {
      this.maxdensity = 0;
      var h = this.params.h;
      var h2 = h*h;
      var h8 = ( h2*h2 )*( h2*h2 );
      var C = 4 * this.mass / 3.14 / h8;
      let posi = vec3.create();
      let posj = vec3.create();
      let dir = vec3.create();
      for (var i = 0; i < this.numSpheres; i++) this.rho[i] = 0;
      for (var i = 0; i < this.numSpheres; ++i) 
      {
          this.rho[i] += 4 * this.mass / 3.14 / h2;
          this.sh.fromData(i, this.data);
          vec3.set(posi, this.sh.pos[0], this.sh.pos[1], this.sh.pos[2]);

          for (var j = i+1; j < this.numSpheres; ++j) 
          {
             this.sh.fromData(j, this.data);
             vec3.set(posj, this.sh.pos[0], this.sh.pos[1], this.sh.pos[2]);

             // dir = posi - posj
             vec3.sub(dir, posi, posj);
             var r2 = vec3.sqrLen(dir); 
             var z = h2-r2;
             if (z > 0) 
             {
                var rho_ij = C*z*z*z;
                this.rho[i] += rho_ij;
                this.rho[j] += rho_ij;
             }
          }
      }  
       
      for (var i = 0; i < this.numSpheres; i++) this.maxdensity = Math.max(this.maxdensity, this.rho[i]);
   }

    normalizeMass()
    {
       this.mass = 1;
       this.computeDensity();
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
        vec3.set(this.accelerations[i], 0, this.params.g, 0);
        /* ASN TODO
        for (var j = 0; j < obstacles.length; j++)
        {
          var b = this.capsuleIntersect(this.positions[i], this.velocities[i], obstacles[j]);
          if (b === true)
          {
            // compute repellent force
            var force = this.capsuleForce(this.positions[i], this.velocities[i], obstacles[j]);
            this.accelerations[i].x += force.x;
            this.accelerations[i].y += force.y;        
          }
        }    
        */
      }
      
      // Constants for interaction term
      /*
      var C0 = mass / 3.14 / ( (h2)*(h2) );
      var Cp = 15*k;
      var Cv = -40*mu;
      
      // Now compute interaction forces
      let posi = vec3.create();
      let posj = vec3.create();
      let dir = vec3.create();
      let dv = vec3.create();
      let scrap1 = vec3.create();
      let scrap2 = vec3.create();
      for (var i = 0; i < this.numSpheres; ++i) 
      {
        var rhoi = this.rho[i];
        this.sh.fromData(i, this.data); // todo: make this support vec3 directly?
        vec3.set(posi, this.sh.pos[0], this.sh.pos[1], this.sh.pos[2]);

        for (var j = i+1; j < this.numSpheres; ++j) 
        {
          this.sh.fromData(j, this.data); // todo: make this support vec3 directly?
          vec3.set(posj, this.sh.pos[0], this.sh.pos[1], this.sh.pos[2]);

          vec3.sub(dir, posi, posj); // dir = posi - posj
          var r2 = vec3.sqrLen(dir);
          if (r2 < h2) 
          {
            var rhoj = this.rho[j];
            var q = Math.sqrt(r2)/h;
            var u = 1-q;
            var w0 = C0 * u/rhoi/rhoj;
            var wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
            var wv = w0 * Cv;

            // dv = vel_i - vel_j
            vec3.sub(dv, this.velocities[i], this.velocities[j]);
            
            // accel_i += wp * dir + wv * dv
            // accel_j -= wp * dir + wv * dv
            vec3.scale(scrap1, dir, wp);
            vec3.scale(scrap2, dv, wv);
            vec3.add(scrap1, scrap1, scrap2);
            vec2.add(this.accelerations[i], this.accelerations[i], scrap1);
            vec2.sub(this.accelerations[j], this.accelerations[i], scrap1);
          }
        }
      }*/
    }

    leapfrogStep()
    {
        let scrap = vec3.create();
        for (var i = 0; i < this.numSpheres; ++i)
        {
            var mag = vec3.len(this.accelerations[i]);
            if (mag > this.params.maxacc) 
            {
                vec3.normalize(this.accelerations[i], this.accelerations[i]);
                vec3.scale(this.accelerations[i], this.accelerations[i], this.params.maxacc);
            }
            vec3.scale(scrap, this.accelerations[i], this.params.dt);
            vec3.add(this.vh[i], this.vh[i], this.accelerations[i]);
       }
       for (var i = 0; i < this.numSpheres; ++i)
       {
            vec3.scale(scrap, this.accelerations[i], this.params.dt * 0.5);
            vec3.add(this.velocities[i], this.vh[i], scrap);
       }
       for (var i = 0; i < this.numSpheres; ++i)
       {
            vec3.scale(scrap, this.vh[i], this.params.dt);

            this.sh.fromData(i, this.data);
            this.sh.pos[0] += scrap[0];
            this.sh.pos[1] += scrap[1];
            this.sh.pos[2] += scrap[2];
            this.sh.toData(i, this.data);
       }
       this.reflectAtBoundaries();
    }

    leapfrogStart()
    {
        var scrap = vec3.create();
        for (var i = 0; i < this.numSpheres; i++)
        {
            // vh = vel + 0.5 * dt * accel
            vec3.scale(scrap, this.accelerations[i], this.params.dt * 0.5);
            vec3.add(this.vh[i], this.velocities[i], scrap);
        }
        for (var i = 0; i < this.numSpheres; ++i)
        {
            // vel = vel + dt * accel
            vec3.scale(scrap, this.accelerations[i], this.params.dt);
            vec3.add(this.velocities[i], this.velocities[i], scrap);
        }
        for (var i = 0; i < this.numSpheres; ++i)
        {
            // pos = pos + dt * vh
            vec3.scale(scrap, this.vh[i], this.params.dt);

            this.sh.fromData(i, this.data);
            this.sh.pos[0] += scrap[0];
            this.sh.pos[1] += scrap[1];
            this.sh.pos[2] += scrap[2];
            this.sh.toData(i, this.data);
        }
        this.reflectAtBoundaries();
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
      // Coefficient of resitiution
      var DAMP = 0.75;
      
      // bounce back along normal direction
      var obsDir = vec3.create();
      vec3.sub(obsDir, obs.b, obs.a); 

      var obsN = vec3.create();
      vec3.set(obsN, -obsDir[1], obsDir[0], 0.0);

      vec3.normalize(obsDir, obsDir);
      vec3.normalize(obsN, obsN);
      
      var velN = vec3.dot(obsN, this.velocities[which]);
      var velT = vec3.dot(obsDir, this.velocities[which]);
      var relfectVelx = -velN * DAMP * (obsN[0]) + velT * (obsDir[0]);
      var relfectVely = -velN * DAMP * (obsN[1]) + velT * (obsDir[1]); 
      
      var vhN = vec3.dot(obsN, this.vh[which]);
      var vhT = vec3.dot(obsDir, this.vh[which]);
      var relfectVhx = -vhN * DAMP * obsN[0] + vhT * obsDir[0];
      var relfectVhy = -vhN * DAMP * obsN[1] + vhT * obsDir[1];
      
      this.velocities[which][0] = relfectVelx;
      this.velocities[which][1] = relfectVely;
      this.vh[which][0] = relfectVhx;
      this.vh[which][1] = relfectVhy;
      console.log(this.vh[which]);
      
      // Damp the velocities
      this.velocities[which][0] *= DAMP; 
      this.vh[which][0] *= DAMP;
      
      // Reflect the position and velocity
      this.sh.fromData(which, this.data);
      var diffx = this.sh.pos[0] - obs.a[0];
      var diffy = this.sh.pos[1] - obs.a[1];
      var diffz = this.sh.pos[2] - obs.a[2];
      var flipDist = diffx * obsN[0] + diffy * obsN[1] + diffz * obsN[2]; // (p-a).dot(obsN) 
      this.sh.pos[0] -= 2*flipDist*obsN[0];
      this.sh.pos[1] -= 2*flipDist*obsN[1];    
      this.sh.toData(which, this.data);
      console.log("REFLECT "+this.sh.pos);

    }

    reflectAtBoundaries()
    {
      // Boundaries of the computational domain
      var XMIN = this.params.minx;
      var XMAX = this.params.maxx;
      var YMIN = this.params.miny;
      var YMAX = this.params.maxy;
      for (var i = 0; i < this.numSpheres; ++i) 
      {
          this.sh.fromData(i, this.data);
          // ASN TODO: Create and save walls
          if (this.sh.pos[0] < XMIN) this.dampReflect(i, new Obstacle(XMIN, YMIN, XMIN, YMAX));
          if (this.sh.pos[0] > XMAX) this.dampReflect(i, new Obstacle(XMAX, YMIN, XMAX, YMAX));
          if (this.sh.pos[1] < YMIN) this.dampReflect(i, new Obstacle(XMIN, YMIN, XMAX, YMIN));
      }
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
    return new SPH2D(numSpheres, maxSpheres, new Params());
}

