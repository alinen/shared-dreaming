// Simple particle system which tracks leap hand
class HandParticlesTrack
{
   constructor(numParticles, maxNumParticles, params)
   {
      this.params = params;
      this.numSpheres = numParticles;
      this.maxNumSpheres = maxNumParticles;
      this.data = null; //positions,vels,rgbs needs to go in data to get rendered;
      this.sh = new SphereHelper();

      this.setupSpheres();
      this.initialized = false;
   }

   setupSpheres()
   {
      this.data = new Float32Array(this.maxNumSpheres * 3 * 4)
      for (var i = 0; i < this.data.length; i++)
      {
         this.data[i] = -20.0;
      }
   }

   init(hand)
   {
      if (!hand.currentFrame) return false;

      for (var i = i; i < this.numSpheres; i++)
      {
         this.sh.fromData(i, this.data);
         this.sh.radius = 0.1;
         this.sh.pos[0] = -2.0 + this.params.r+cellj*margin;
         this.sh.pos[1] = this.params.r-celli*margin+1;
         this.sh.pos[2] = -2.0;
         this.sh.rgb[0] = 1.0;
         this.sh.rgb[1] = 0.0;
         this.sh.rgb[2] = 0.0;
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

    update(dt, hand)  // ASN TODO: Fixed framerate or application framerate?
    {
      if (!this.initialized)
      {
        this.initialized = this.init(hand);
      }

      if (!this.paused && this.initialized && hand.getCurrentFrame())
      {
            var idx = 0;
            this.writePosition(idx++, hand.fingerJoint('left', 'index', 1));
            this.writePosition(idx++, hand.fingerJoint('left', 'index', 2));
            this.writePosition(idx++, hand.fingerJoint('left', 'index', 3));

            this.writePosition(idx++, hand.fingerJoint('left', 'middle', 1));
            this.writePosition(idx++, hand.fingerJoint('left', 'middle', 2));
            this.writePosition(idx++, hand.fingerJoint('left', 'middle', 3));

            this.writePosition(idx++, hand.fingerJoint('left', 'ring', 1));
            this.writePosition(idx++, hand.fingerJoint('left', 'ring', 2));
            this.writePosition(idx++, hand.fingerJoint('left', 'ring', 3));

            this.writePosition(idx++, hand.fingerJoint('left', 'pinky', 1));
            this.writePosition(idx++, hand.fingerJoint('left', 'pinky', 2));
            this.writePosition(idx++, hand.fingerJoint('left', 'pinky', 3));

            this.writePosition(idx++, hand.fingerJoint('left', 'thumb', 1));
            this.writePosition(idx++, hand.fingerJoint('left', 'thumb', 2));
            this.writePosition(idx++, hand.fingerJoint('left', 'thumb', 3));

            this.writePosition(idx++, hand.fingerJoint('right', 'index', 1));
            this.writePosition(idx++, hand.fingerJoint('right', 'index', 2));
            this.writePosition(idx++, hand.fingerJoint('right', 'index', 3));

            this.writePosition(idx++, hand.fingerJoint('right', 'middle', 1));
            this.writePosition(idx++, hand.fingerJoint('right', 'middle', 2));
            this.writePosition(idx++, hand.fingerJoint('right', 'middle', 3));

            this.writePosition(idx++, hand.fingerJoint('right', 'ring', 1));
            this.writePosition(idx++, hand.fingerJoint('right', 'ring', 2));
            this.writePosition(idx++, hand.fingerJoint('right', 'ring', 3));

            this.writePosition(idx++, hand.fingerJoint('right', 'pinky', 1));
            this.writePosition(idx++, hand.fingerJoint('right', 'pinky', 2));
            this.writePosition(idx++, hand.fingerJoint('right', 'pinky', 3));

            this.writePosition(idx++, hand.fingerJoint('right', 'thumb', 1));
            this.writePosition(idx++, hand.fingerJoint('right', 'thumb', 2));
            this.writePosition(idx++, hand.fingerJoint('right', 'thumb', 3));
       }
    }
};

function createSystem(numSpheres, maxSpheres)
{
    return new HandParticlesTrack(numSpheres, maxSpheres, new Params());
}

