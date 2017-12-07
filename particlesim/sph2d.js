class Obstacle
{
    constructor(start, end, radius)
    {
        this.a = start;
        this.b = end
        this.r = radius;
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
      this.obstacles = [];
      this.sh = new SphereHelper();

      // Boundaries of the computational domain
      var leftWallStart = vec3.fromValues(this.params.min[0], this.params.min[1], 0); // todo: make a plane
      var leftWallEnd = vec3.fromValues(this.params.min[0], this.params.max[1], 0); // todo 
      this.leftWall = new Obstacle(leftWallStart, leftWallEnd, 0);
            
      var rightWallStart = vec3.fromValues(this.params.max[0], this.params.min[1], 0); // todo: make a plane
      var rightWallEnd = vec3.fromValues(this.params.max[0], this.params.max[1], 0); // todo 
      this.rightWall = new Obstacle(rightWallStart, rightWallEnd, 0);

      var bottomWallStart = vec3.fromValues(this.params.min[0], this.params.min[1], 0); // todo: make a plane
      var bottomWallEnd = vec3.fromValues(this.params.max[0], this.params.min[1], 0); // todo 
      this.bottomWall = new Obstacle(bottomWallStart, bottomWallEnd, 0);

      var topWallStart = vec3.fromValues(this.params.min[0], this.params.max[1], 0); // todo: make a plane
      var topWallEnd = vec3.fromValues(this.params.max[0], this.params.max[1], 0); // todo 
      this.topWall = new Obstacle(topWallStart, topWallEnd, 0);
      
      // create some obstacles to test
      this.pushObstacle(vec3.fromValues(-2,1,-2), vec3.fromValues(2,-2,-2));
      
      this.setupSpheres();
      this.init();
      this.reflectAtBoundaries();
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
      }
   }

   init()
   {

      var numcols = 18;
      var margin = this.params.r * 2.2;
      for (var i = 0; i < this.numSpheres; i++)
      {
         var cellj = i % numcols;
         var celli = Math.floor(i/numcols);

         this.sh.fromData(i, this.data);
         this.sh.radius = this.params.r; 
         this.sh.pos[0] = -2.0 + this.params.r+cellj*margin;
         this.sh.pos[1] = this.params.r-celli*margin+1;
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

            //this.sh.fromData(0, this.data);
            //console.log(this.sh.pos[1]);
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
      var spherePos = vec3.create();
      for (var i = 0; i < this.numSpheres; ++i) 
      {
        this.sh.fromData(i, this.data);
        vec3.set(spherePos, this.sh.pos[0], this.sh.pos[1], this.sh.pos[2]);
        
        // put gravity at the center
        var force = vec3.create();
        var r = vec3.len(spherePos);
        vec3.scale(force, spherePos, this.params.g); 
        vec3.set(this.accelerations[i], force[0], force[1], 0);

        //vec3.set(this.accelerations[i], 0, this.params.g, 0);
        //vec3.set(this.accelerations[i], 0, 0, 0);

        for (var j = 0; j < this.obstacles.length; j++)
        {
            // compute repellent force
            var force = this.capsuleForce(spherePos, this.velocities[i], this.obstacles[j]);
            vec3.add(this.accelerations[i], this.accelerations[i], force);
        }    
      }
      
      // Constants for interaction term
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
      }
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
      var ba = vec3.create(); 
      var pa = vec3.create(); 
      var pb = vec3.create(); 
      
      vec3.sub(ba, obs.b, obs.a);
      vec3.sub(pa, p, obs.a);
      vec3.sub(pb, p, obs.b);
      
      var len = vec3.len(ba);
      var dir = vec3.create();
      vec3.normalize(dir, ba);
      
      var combinedRadius = obs.r + this.params.r;
        
      var force = vec3.create();
      var dirT = vec3.dot(pa, dir);
      if (dirT > 0 && dirT < len) // might intercept middle part
      {
          // perp = pa - dir * dirT
          var perp = vec3.create();
          vec3.scale(perp, dir, -dirT); 
          vec3.add(perp, pa, perp); 
          var perpSq = vec3.sqrLen(perp);
          if (perpSq < combinedRadius*combinedRadius) 
          {
              var scale = this.params.kIntersect * (combinedRadius*combinedRadius - perpSq);
              vec3.scale(force, perp, scale);
          }
      }
      else
      {
          var magsq = vec3.sqrLen(pa);
          if (magsq < combinedRadius*combinedRadius) 
          {
            var scale = this.params.kIntersect * (combinedRadius*combinedRadius - magsq);
            vec3.scale(force, pa, scale);
          }
        
          magsq = vec3.sqrLen(pb);
          if (magsq < combinedRadius*combinedRadius)
          {
            var scale = this.params.kIntersect * (combinedRadius*combinedRadius - magsq);
            vec3.scale(force, pb, scale);
          }
      }

      vec3.scale(force, force, -1.0);
      return force;
    }

    capsuleIntersect(p, v, obs)
    {
      var ba = vec3.create(); 
      var pa = vec3.create(); 
      var pb = vec3.create(); 
      
      vec3.sub(ba, obs.b, obs.a);
      vec3.sub(pa, p, obs.a);
      vec3.sub(pb, p, obs.b);
      
      var len = vec3.len(ba);
      var dir = vec3.create();
      vec3.normalize(dir, ba);
      
      var combinedRadius = obs.r + this.params.r;
        
      var dirT = vec3.dot(pa, dir);
      if (dirT > 0 && dirT < len) // might intercept middle part
      {
          // perp = pa - dir * dirT
          var perp = vec3.create();
          vec3.scale(perp, dir, -dirT); 
          vec3.add(perp, pa, perp); 
          var perpSq = vec3.sqrLen(perp);
          if (perpSq < combinedRadius*combinedRadius) 
          {
              return true;
          }
      }
      else
      {
          var magsq = vec3.sqrLen(pa);
          if (magsq < combinedRadius*combinedRadius) return true;
        
          magsq = vec3.sqrLen(pb);
          if (magsq < combinedRadius*combinedRadius) return true;
      }
      return false;
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

    }

    reflectAtBoundaries()
    {
        for (var i = 0; i < this.numSpheres; ++i) 
        {
            this.sh.fromData(i, this.data);
            if (this.sh.pos[0] < this.params.min[0]) 
            {
                this.dampReflect(i, this.leftWall); 
            }
            if (this.sh.pos[0] > this.params.max[0]) 
            {
                this.dampReflect(i, this.rightWall);
            }
            if (this.sh.pos[1] < this.params.min[1]) 
            {
                this.dampReflect(i, this.bottomWall); 
            }
            if (this.sh.pos[1] > this.params.max[1]) 
            {
                this.dampReflect(i, this.topWall); 
            }
        }
    }

    clearObstacles()
    {
        this.obstacles = [];
    }

    pushObstacle(start, end, r = 0.15)
    {
        this.obstacles.push(new Obstacle(start,end, r));
    }

    popObstacle()
    {
        if (this.obstacles.length > 0) 
        {
            this.obstacles.splice(this.obstacles.length - 1, 1);
        }
    }

};

function createSystem(numSpheres, maxSpheres)
{
    return new SPH2D(numSpheres, maxSpheres, new Params());
}

