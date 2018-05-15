class SphereHelper
{
   constructor()
   {
      this.pos = [0.0,0.0,0.0];
      this.vel = [0.0,0.0,0.0];
      this.rgb = [0.0,0.0,0.0];
      this.radius = 0.0;
   }

   // ASN: maybe avoid this copying
   fromData(idx, data)
   {
     var i = idx * 4;
     this.pos[0] = data[i + 0];
     this.pos[1] = data[i + 1];
     this.pos[2] = data[i + 2];
     this.radius = data[i + 3];
   }

   toData(sphereIndex, data)
   {
     var i = sphereIndex * 4;

     // pos
     data[i + 0] = this.pos[0]; // x
     data[i + 1] = this.pos[1]; // y
     data[i + 2] = this.pos[2]; // z
     data[i + 3] = this.radius; // radius
   }

   update(idx, dt, data, hand)
   {
     var wind = 0;
     if (hand !== undefined && hand.getCurrentFrame())
     {
       var wrist = hand.wristPosition('left');
       if (wrist[2] < 0.0001)
       {
         wind = wrist[0] * 0.5;
       }
     }

     this.fromData(idx, data);

     var finished = false;
     if (this.pos[1] > 2.0) 
     {
        var y = -1.5 - this.radius - getRandom(0.05, 0.1); // y + radius
        this.pos[1] = y; // put this sphere behind the sphere it's following
        //console.log("reset: "+y +" "+previ
     }
     this.pos[0] = this.vel[0] + Math.sin(elapsedTime+idx)*0.15 + wind; // hack, vel[0] has starting x
     this.pos[1] += this.vel[1] * dt;
     this.pos[2] += this.vel[2] * dt;
     this.toData(idx, data);

     //console.log(this.pos);

     return finished;
   }
}

