class Obstacle
{
    constructor(start, end, radius)
    {
        this.a = start;
        this.b = end
        this.r = radius;
    }
}

class HandParticles
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
      var leftWallStart = vec3.fromValues(this.params.min[0], this.params.min[1], 0); 
      var leftWallEnd = vec3.fromValues(this.params.min[0], this.params.max[1], 0); 
      this.leftWall = new Obstacle(leftWallStart, leftWallEnd, 0);
            
      var rightWallStart = vec3.fromValues(this.params.max[0], this.params.min[1], 0); 
      var rightWallEnd = vec3.fromValues(this.params.max[0], this.params.max[1], 0); 
      this.rightWall = new Obstacle(rightWallStart, rightWallEnd, 0);

      var bottomWallStart = vec3.fromValues(this.params.min[0], this.params.min[1], 0); 
      var bottomWallEnd = vec3.fromValues(this.params.max[0], this.params.min[1], 0); 
      this.bottomWall = new Obstacle(bottomWallStart, bottomWallEnd, 0);

      var topWallStart = vec3.fromValues(this.params.min[0], this.params.max[1], 0); 
      var topWallEnd = vec3.fromValues(this.params.max[0], this.params.max[1], 0); 
      this.topWall = new Obstacle(topWallStart, topWallEnd, 0);
      
      this.setupSpheres();
      this.initialized = false;
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

   init(hand)
   {
      if (!hand.currentFrame) return false;

      /*
      for (var i = i; i < this.numSpheres; i++)
      {
         this.sh.fromData(i, this.data);
         this.sh.radius = this.params.r; 
         this.sh.pos[0] = -2.0 + this.params.r+cellj*margin;
         this.sh.pos[1] = this.params.r-celli*margin+1;
         this.sh.pos[2] = -2.0;
         this.sh.toData(i, this.data);
         
         vec3.set(this.velocities[i], 0,0,0);
         vec3.set(this.accelerations[i], 0,0,0);
         vec3.set(this.vh[i], 0,0,0);
      }*/
      
      //this.reflectAtBoundaries();
      //this.computeAccel(hand);
      //this.leapfrogStart();
      return true;
   }

   writePosition(idx, pos)
   {
       if (pos[2] > -0.001) return; // hack for missing data
       this.sh.fromData(idx, this.data);
       this.sh.pos[0] = pos[0];
       this.sh.pos[1] = pos[1];
       this.sh.pos[2] = pos[2];
       //console.log(pos[0]+" "+pos[1]+" "+pos[2]);
       this.sh.toData(idx, this.data);
   }

    update(dt, hand)  // ASN TODO: Fixed framerate or application framerate?
    {
      if (!this.initialized)
      {
        this.initialized = this.init(hand);
      }

      if (!this.paused && this.initialized && hand.getCurrentFrame())
      {
            this.writePosition(0, hand.wristPosition('left'));
            this.writePosition(1, hand.fingerJoint('left', 'index', 3));
            this.writePosition(2, hand.fingerJoint('left', 'middle', 3));
            this.writePosition(3, hand.fingerJoint('left', 'ring', 3));
            this.writePosition(4, hand.fingerJoint('left', 'pinky', 3));
            this.writePosition(5, hand.fingerJoint('left', 'thumb', 3));

            this.writePosition(6, hand.wristPosition('right'));
            this.writePosition(7, hand.fingerJoint('right', 'index', 3));
            this.writePosition(8, hand.fingerJoint('right', 'middle', 3));
            this.writePosition(9, hand.fingerJoint('right', 'ring', 3));
            this.writePosition(10, hand.fingerJoint('right', 'pinky', 3));
            this.writePosition(11, hand.fingerJoint('right', 'thumb', 3));
       }
       //this.computeAccel(hand);
       //this.leapfrogStep();
    }

    computeAccel(hand)
    {
      var spherePos = vec3.create();
      for (var i = 0; i < numHandData && i < this.numSpheres; ++i) 
      {
        if (hand.getCurrentFrame())
        {
          hand.fromData(i);
          var handPos = vec3.fromValues(hand.pos[0], hand.pos[1], hand.pos[2]);

          this.sh.fromData(i, this.data);
          vec3.set(spherePos, this.sh.pos[0], this.sh.pos[1], this.sh.pos[2]);

          var distance = vec3.distance(handPos, spherePos);
          //if (distance > 0.1)
          {
            this.sh.pos[0] = hand.pos[0];
            this.sh.pos[1] = hand.pos[1];
            this.sh.pos[2] = hand.pos[2];
            this.sh.toData(i, this.data);

            //console.log(this.sh.pos[0]);

            vec3.set(this.accelerations[i], 0, 0, 0);
          }
          /*
          else
          {
            var force = vec3.create();
            vec3.sub(force, handPos, spherePos); 
            vec3.scale(force, force, this.params.kIntersect); 
            vec3.set(this.accelerations[i], force[0], force[1], force[2]);
          }*/

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
};

function createSystem(numSpheres, maxSpheres)
{
    return new HandParticles(numSpheres, maxSpheres, new Params());
}

