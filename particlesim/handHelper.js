class HandHelper
{
  constructor() {
    this.palmPos = [0.0, 0.0, 0.0];
    this.radius = 0.0;
    this.rgba = [0.0, 0.0, 0.0, 0.0];

    this.data = [];
    this.frames = null;
    this.frameIndex = 0;
    this.radius = 0.1;
    this.rgba = [1,1,1,1];

    this.init();
  }

  init() {
    this.data = new Float32Array(numHandsData * 2 * 4);
    for (int i = 0; i < this.data.length; i++) this.data[i] = 0.0;
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
    var normX = (positionValues.x - this.center.x) / this.boundingBox.x
    var normY = (positionValues.y - this.center.y) / this.boundingBox.y
    var normZ = (positionValues.z - this.center.z) / this.boundingBox.z
    return [normX, normY, normZ];
  }

  update(dt) {

    var nameMap = ["thumb", "index", "middle", "ring", "pinky"];
    if (this.frames && this.frameIndex < this.frames.length) {

        var frame = frames[this.frameIndex];

        var idx = 0;
        var palmPos = this.normalizeValues(frame['right']['palmPosition']);
        this.toData(idx++, this.data, palmPos, this.rgba, this.radius);

        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 4; j++) 
            {
                var pos = this.normalizeValues(frame['right'][nameMap[i]][j];
                this.toData(idx++, this.data, pos, this.rgba, this.radius);
            }
        }

        palmPos = this.normalizeValues(frame['left']['palmPosition']);
        this.toData(idx++, this.data, palmPos, this.rgba, this.radius);

        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 4; j++) 
            {
                var pos = this.normalizeValues(frame['left'][nameMap[i]][j];
                this.toData(idx++, this.data, pos, this.rgba, this.radius);
            }
        }

        this.frameIndex = (this.frameIndex + 1) % this.frames.length;
    }
  }

  getFrameData() {
    var url = "https://raw.githubusercontent.com/alinen/shared-dreaming/master/recordings_pinch-57fps-57fps.json";
    var xhr = new XMLHttpRequest();
    var _this = this;

    xhr.onreadystatechange = function() {
      if (xhr.readyState === xhr.DONE) {
        if (xhr.status === 200 || xhr.status === 0) {
          if (xhr.responseText) {
            _this.frameData = JSON.parse(xhr.responseText);
            _this.frames = _this.frameData.frames;
            _this.center = _this.frameData.center;
            _this.boundingBox = _this.frameData.bbox;
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
