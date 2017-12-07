class HandHelper
{
  constructor(jsonName = undefined) {
    this.palmPos = [0.0, 0.0, 0.0];
    this.radius = 0.0;
    this.rgba = [0.0, 0.0, 0.0, 0.0];

    this.data = [];
    this.frames = null;
    this.frameIndex = 0;
    this.radius = 0.1;
    this.rgba = [1,1,1,1];

    this.init();

    if (jsonName !== undefined)
    {
        this.getFrameData(jsonName);
    }
  }

  init() {
    this.data = new Float32Array(numHandData * 2 * 4);
    for (var i = 0; i < this.data.length; i++) this.data[i] = 0.0;
  }

  toData(idx, data, pos, rgba, radius) {
    var i = idx * 2 * 4;

    data[i + 0] = pos[0]
    data[i + 1] = pos[1]
    data[i + 2] = pos[2]

    data[i + 3] = radius;

    data[i + 4] = rgba[0];
    data[i + 5] = rgba[1];
    data[i + 6] = rgba[2];
    data[i + 7] = rgba[3];
  }

  normalizeValues(positionValues) {
// this.center (3) [0, 200, 0]
// this.boundingBox (3) [235.247, 235.247, 147.751]
    var worldCenter = [0.0, 0.0, -2.0];
    var normX = ((positionValues[0] - this.center[0]) / this.boundingBox[0]);
    var normY = ((positionValues[1] - this.center[1]) / this.boundingBox[1]);
    var normZ = ((positionValues[2] - this.center[2]) / this.boundingBox[2]);

    // rotate hand = - 90
    //var angle = -90 * 180.0/Math.PI;
    //var x = Math.cos(angle) * normX - Math.sin(angle) * normY;
    //var y = Math.sin(angle) * normX + Math.sin(angle) * normY;

    normX = normX * 2.0 + worldCenter[0];
    normY = normY * 2.0 + worldCenter[1];
    normZ = normZ * 2.0 + worldCenter[1];

    return [normX, normY, normZ];
  }

  update(dt) {

    var nameMap = ["thumb", "index", "middle", "ring", "pinky"];
    if (this.frames && this.frameIndex < this.frames.length) {

        var frame = this.frames[this.frameIndex];

        var idx = 0;
        if (Object.keys(frame['right']).length !== 0)
        {
            var palmPos = this.normalizeValues(frame['right']['palmPosition']);
            this.toData(idx++, this.data, palmPos, this.rgba, this.radius);

            for (var i = 0; i < 5; i++)
            {
                for (var j = 0; j < 4; j++) 
                {
                    var pos = this.normalizeValues(frame['right'][nameMap[i]][j]);
                    this.toData(idx++, this.data, pos, this.rgba, this.radius);
                }
            }
        }

        if (Object.keys(frame['left']).length !== 0)
        {
            idx = 4 * 5 + 1;
            var palmPos = this.normalizeValues(frame['left']['palmPosition']);
            this.toData(idx++, this.data, palmPos, this.rgba, this.radius);
            console.log("Left PALM "+palmPos);

            for (var i = 0; i < 5; i++)
            {
                for (var j = 0; j < 4; j++) 
                {
                    var pos = this.normalizeValues(frame['left'][nameMap[i]][j]);
                    this.toData(idx++, this.data, pos, this.rgba, this.radius);
                }
            }
        }

        this.frameIndex = (this.frameIndex + 1) % this.frames.length;
    }
  }

  getFrameData(jsonName) {
    var url = jsonName; 
    var xhr = new XMLHttpRequest();
    var _this = this;

    xhr.onreadystatechange = function() {
      if (xhr.readyState === xhr.DONE) {
        if (xhr.status === 200 || xhr.status === 0) {
          if (xhr.responseText) {
            _this.frameData = JSON.parse(xhr.responseText);
            _this.frames = _this.frameData.frames;
            _this.center = _this.frameData.center;
            _this.boundingBox = _this.frameData.boundingBox;
            _this.update(0, []);
          } else {
            console.error('Leap Playback: "' + url + '" seems to be unreachable or the file is empty.');
          }
        } else {
          console.error('Leap Playback: Couldn\'t load "' + url + '" (' + xhr.status + ')');
        }
      }
    }
    xhr.open("GET", url, true);
    xhr.send(null);
  }
}
