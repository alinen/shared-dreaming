// Simple particle system which tracks leap hand
// In this example, we create additional particles at random positions around the hands
class HandParticlesTrack1
{
   constructor(numParticles, maxNumParticles, params)
   {
      this.params = params;
      this.numSpheres = numParticles;
      this.maxNumSpheres = maxNumParticles;
      this.data = null; //positions,vels,rgbs needs to go in data to get rendered;
      this.sh = new SphereHelper();
      this.min = new Float32Array(3);
      this.max = new Float32Array(3);

      this.setupSpheres();
      this.initialized = false;
   }

   setupSpheres()
   {
      this.data = new Float32Array(this.maxNumSpheres * 3 * 4)
      for (var i = 0; i < this.data.length; i++)
      {
         this.data[i] = 20.0;
      }
   }

   init(hand)
   {
      if (!hand.currentFrame) return false;

      for (var i = i; i < this.numSpheres; i++)
      {
         this.sh.fromData(i, this.data);
         this.sh.radius = 0.5;
         this.sh.pos[0] = -2.0 + this.params.r+cellj*margin;
         this.sh.pos[1] = this.params.r-celli*margin+1;
         this.sh.pos[2] = -2.0;
         //this.sh.rgb[0] = 1.0; // not used
         //this.sh.rgb[1] = 0.0;
         //this.sh.rgb[2] = 0.0;
         this.sh.toData(i, this.data);
      }
      
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

   placeAlongSegment(idx, start, end, num)
   {
      for (var i = 0; i < num; i++)
      {
          var alpha = i/(1.0*num);
          var pos = vec3.fromValues(start[0], start[1], start[2]);
          var dir = vec3.create();
          vec3.sub(dir, end, start);
          vec3.scale(dir, dir, alpha);
          vec3.add(pos, pos, dir); 
          this.writePosition(idx++, pos); 

          this.min[0] = Math.min(this.min[0], pos[0]);
          this.min[1] = Math.min(this.min[1], pos[1]);
          this.min[2] = Math.min(this.min[2], pos[2]);
          
          this.max[0] = Math.max(this.max[0], pos[0]);
          this.max[1] = Math.max(this.max[1], pos[1]);
          this.max[2] = Math.max(this.max[2], pos[2]);
      }
      return idx;
   }

   updateFinger(idx, hand, which, finger, num)
   {
      idx = this.placeAlongSegment(idx, hand.fingerJoint(which, finger, 1), hand.fingerJoint(which, finger, 2), num);
      idx = this.placeAlongSegment(idx, hand.fingerJoint(which, finger, 2), hand.fingerJoint(which, finger, 3), num);
      this.writePosition(idx++, hand.fingerJoint(which, finger, 3));
      return idx;
   }

   updateHand(idx, hand, which)
   {
      idx = this.updateFinger(idx, hand, which, 'index', 10);
      idx = this.updateFinger(idx, hand, which, 'middle', 10);
      idx = this.updateFinger(idx, hand, which, 'ring', 10);
      idx = this.updateFinger(idx, hand, which, 'pinky', 10);
      idx = this.updateFinger(idx, hand, which, 'thumb', 10);

      var src = vec3.fromValues(0,-2,2);
      idx = this.placeAlongSegment(idx, src, hand.fingerJoint(which, 'thumb', 1), 100);
      //idx = this.placeAlongSegment(idx, src, hand.fingerJoint(which, 'middle', 1), 20);
      //idx = this.placeAlongSegment(idx, src, hand.fingerJoint(which, 'ring', 1), 20);
      //idx = this.placeAlongSegment(idx, src, hand.fingerJoint(which, 'pinky', 1), 20);
      //idx = this.placeAlongSegment(idx, src, hand.fingerJoint(which, 'index', 1), 20);
      idx = this.placeAlongSegment(idx, hand.fingerJoint(which, 'thumb', 1), hand.fingerJoint(which, 'index', 1), 10);
      idx = this.placeAlongSegment(idx, hand.fingerJoint(which, 'thumb', 1), hand.fingerJoint(which, 'middle', 1), 10);
      idx = this.placeAlongSegment(idx, hand.fingerJoint(which, 'thumb', 1), hand.fingerJoint(which, 'ring', 1), 10);
      idx = this.placeAlongSegment(idx, hand.fingerJoint(which, 'thumb', 1), hand.fingerJoint(which, 'pinky', 1), 10);
      return idx;
   }

   update(dt, hand)  // ASN TODO: Fixed framerate or application framerate?
   {
      if (!this.initialized)
      {
          this.initialized = this.init(hand);
      }

      if (!this.paused && this.initialized && hand.getCurrentFrame())
      {
          this.min[0] = this.min[1] = this.min[2] = 99.0;
          this.max[0] = this.max[1] = this.max[2] = -99.0;

          var idx = 0;
          idx = this.updateHand(idx, hand, 'left');
          idx = this.updateHand(idx, hand, 'right');
          if (idx > this.maxNumSpheres) console.log("WARNING: "+idx+" > "+this.maxNumSpheres);

          console.log(this.min);
          console.log(this.max);
       }
   }
};

function createSystem(numSpheres, maxSpheres)
{
    return new HandParticlesTrack1(numSpheres, maxSpheres, new Params());
}

