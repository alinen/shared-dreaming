class HandHelper
{
  constructor()
  {
    this.pos = [0.0, 0.0, 0.0];
    this.radius = 0.0;
    this.rgba = [0.0, 0.0, 0.0, 0.0];
  }

  fromData(idx, data)
  {
    var i = idx * 2 * 4;
    this.pos[0] = data[i + 0];
    this.pos[1] = data[i + 1];
    this.pos[2] = data[i + 2];
    this.radius = data[i + 3];

    this.rgba[0] = data[i + 4];
    this.rgba[1] = data[i + 5];
    this.rgba[2] = data[i + 6];
    this.rgba[3] = data[i + 7];
  }

  toData(handIndex, data)
  {
    var i = handIndex * 2 * 4;

    data[i + 0] = this.pos[0];
    data[i + 1] = this.pos[1];
    data[i + 2] = this.pos[2];
    data[i + 3] = this.radius;

    data[i + 4] = this.rgb[0];
    data[i + 5] = this.rgb[1];
    data[i + 6] = this.rgb[2];
    data[i + 7] = this.rgb[3];
  }

  update(idx, dt, data)
  {
    this.fromData(idx, data);

    // parse JSON object /////////////

    var finished = false;
    // if (this.pos[1] > 2.0)
    // {
    //   var y = -1.5 - this.radius - getRandom(0.05, 0.1); // y + radius
    //   this.pos[1] = y; // put this sphere behind the sphere it's following
    //   //console.log("reset: "+y +" "+previ
    // }
    // this.pos[0] = this.vel[0] + Math.sin(elapsedTime+idx)*0.15; // hack, vel[0] has starting x
    // this.pos[1] += this.vel[1] * dt;
    // this.pos[2] += this.vel[2] * dt;
    // this.toData(idx, data);

    return finished;
  }
}