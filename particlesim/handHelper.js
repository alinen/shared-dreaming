class HandHelper
{
  constructor() {
    this.palmPos = [0.0, 0.0, 0.0];
    this.radius = 0.0;
    this.rgba = [0.0, 0.0, 0.0, 0.0];

    this.data = [];
    this.metadata;
    this.timeDelta;
    this.frames;
    this.frameIndex = 1;

    this.setDefaultValues();
  }

  setDefaultValues() {
    this.data = new Float32Array(numHands * 2 * 4);
    var i = 0;
    this.data[i + 0] = 1.0;
    this.data[i + 1] = 0.0;
    this.data[i + 2] = 0.0;

    this.data[i + 3] = 1.0;

    this.data[i + 4] = 1.0;
    this.data[i + 5] = 0.0;
    this.data[i + 6] = 0.0;
    this.data[i + 7] = 1.0;
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

  update(dt, data) {
    if (this.frames && this.frames[this.frameIndex++]) {
      if (this.frameIndex > this.metadata.frames) {
        this.frameIndex = 1;
      }
      if (this.frames[this.frameIndex]) {
        // this.rgba = this.frames[this.frameIndex][2][0][4]
        var rgba = this.frames[this.frameIndex][2][0][4].map(function(n) { return n/255.0 % 1.0});
        console.log('this.rgba', this.rgba)

        var palm = this.normalizeValues(this.frames[this.frameIndex][2][0][4]);
        this.toData(RT_PALM, this.data, palm, rgba, radius);



        this.radius = 1.0;
      }
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
            _this.metadata = _this.frameData.metadata;
            _this.timeDelta = _this.metadata.frames / _this.metadata.frameRate / 30;
            _this.frames = _this.frameData.frames;
            _this.center = _this.frames[1][4][0];
            _this.boundingBox = _this.frames[1][4][1];

            _this.update(1, []);
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