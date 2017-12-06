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

    this.jsIsDumb();
  }

  jsIsDumb() {
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

  fromData(idx, data) {
    var i = idx * 2 * 4;

    this.palmPos[0] = this.data[i + 0];
    this.palmPos[1] = this.data[i + 1];
    this.palmPos[2] = this.data[i + 2];
    this.radius = this.data[i + 3];

    this.rgba[0] = this.data[i + 4];
    this.rgba[1] = this.data[i + 5];
    this.rgba[2] = this.data[i + 6];
    this.rgba[3] = this.data[i + 7];
  }

  toData(idx, data) {
    var i = idx * 2 * 4;

    // this.data[i + 0] = 0.1;
    // this.data[i + 1] = 0.3;
    // this.data[i + 2] = 0.3;

    this.data[i + 0] = 1.0;
    this.data[i + 1] = 0.0;
    this.data[i + 2] = 0.0;

    this.data[i + 3] = 1.0;

    this.data[i + 4] = 1.0;
    this.data[i + 5] = 0.0;
    this.data[i + 6] = 0.0;
    this.data[i + 7] = 1.0;
    // this.data[i + 7] = this.rgba[3];
  }

  update(dt, data) {
    if (this.frames && this.frames[this.frameIndex++]) {
      if (this.frames[this.frameIndex]) {
        this.palmPos = this.frames[this.frameIndex][2][0][4];
        this.radius = 1.0;
        this.rgba = [0.0, 0.0, 0.0, 1.0];
        this.toData(0);
        // console.log('this.data', this.data)
      }
    }
  }

  // getPalmData() {

  // }

  normalizeHandData() {

  }

  getNextHandFrame() {

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
            // initHandTexture();
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