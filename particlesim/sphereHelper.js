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
     var i = idx * 3 * 4;
     this.pos[0] = data[i + 0];
     this.pos[1] = data[i + 1];
     this.pos[2] = data[i + 2];
     this.radius = data[i + 3];

     this.vel[0] = data[i + 4];
     this.vel[1] = data[i + 5];
     this.vel[2] = data[i + 6];

     this.rgb[0] = data[i + 8];
     this.rgb[1] = data[i + 9];
     this.rgb[2] = data[i + 10];
   }

   toData(sphereIndex, data)
   {
     var i = sphereIndex * 3 * 4;

     // pos
     data[i + 0] = this.pos[0]; // x
     data[i + 1] = this.pos[1]; // y
     data[i + 2] = this.pos[2]; // z
     data[i + 3] = this.radius; // radius
     // velocity
     data[i + 4] = this.vel[0]; // v x
     data[i + 5] = this.vel[1]; // v y
     data[i + 6] = this.vel[2]; // v z
     data[i + 7] = 0.0; // free for now
     // color
     data[i + 8] = this.rgb[0]; // r
     data[i + 9] = this.rgb[1]; // g
     data[i + 10] = this.rgb[2]; // b
     data[i + 11] = 1.0; // free for now      
   }

   update(idx, dt, data)
   {
     this.fromData(idx, data);

     var finished = false;
     if (this.pos[1] > 2.0) 
     {
        var y = -1.5 - this.radius - getRandom(0.05, 0.1); // y + radius
        this.pos[1] = y; // put this sphere behind the sphere it's following
        //console.log("reset: "+y +" "+previ
     }
     this.pos[0] = this.vel[0] + Math.sin(elapsedTime+idx)*0.15; // hack, vel[0] has starting x
     this.pos[1] += this.vel[1] * dt;
     this.pos[2] += this.vel[2] * dt;
     this.toData(idx, data);

     return finished;
   }
}

