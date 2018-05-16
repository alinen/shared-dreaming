class Bubbles
{
  constructor(numSpheres, maxNumSpheres)
  {
    this.maxNumSpheres = maxNumSpheres;
    this.numSpheres = numSpheres;
    this.textureSize = maxNumSpheres * 3;
    this.data = [];
    this.sh = new SphereHelper();

    this.setupSpheres();
  }

  update(dt)
  {
    for (var i = 0; i < this.numSpheres; ++i)
    {
      this.sh.update(i, dt, this.data);
    }
  }

  setupSpheres()
  {
    this.data = new Float32Array(this.maxNumSpheres * 3 * 4)
    for (var i = 0; i < this.data.length; i++)
    {
      this.data[i] = 0.0;
    }

    var prevR = 0.0;
    var averageX = -2.0;
    var posY = -3.0;

    this.sh.vel = [0.0, 1.0, 0.0];
    this.sh.pos = [0.0, 0.0, -2.0];
    this.sh.rgb = [0.5, 0.0, 0.0];
    this.sh.radius = 0.25;

    for (var i = 0; i < this.numSpheres; ++i)
    {
      this.sh.radius = getRandom(0.1, 0.3);
      posY -= (prevR + this.sh.radius + getRandom(0.05, 0.1));
      this.sh.vel[0] = getRandom(0, 0.1) + averageX;
      this.sh.vel[1] = getRandom(0.8, 1.5);
      this.sh.pos[1] = posY;
      this.sh.pos[2] = -2.0;
      this.sh.rgb[0] = getRandom(0, 1);
      this.sh.rgb[1] = getRandom(0, 1);
      this.sh.rgb[2] = getRandom(0, 1);
      prevR = this.sh.radius;
      this.sh.toData(i, this.data);

      if (i > 0 && i % 10 == 0)
      {
        averageX += 1.0 + getRandom(-0.1, 0.1);
        prevR = 0.0;
        posY = -3.0;
        console.log(averageX);
      }
    }
  }
}

function createSystem(numSpheres, maxSpheres)
{
  return new Bubbles(numSpheres, maxSpheres);
}

