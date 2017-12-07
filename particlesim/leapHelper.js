// to make working with angles easy
var Deg2Rad = Math.PI / 180;
var Rad2Deg = 1 / Deg2Rad;
var debugPrint = false;

function vectorToString(v) {
  return '(' + v[0] + ', ' + v[1] + ', ' + v[2] + ')';
}

var startdate = new Date();
var startTime = startdate.getTime() * 0.001;


Leap.loop({

// frame callback is run before individual frame components
frame: function(frame) {
  var myFrame  = {
    "timestamp": 1234,
    "right": {},
    "left": {}
  };
  var handString = '';

  var date = new Date();
  var newTime = date.getTime() * 0.001;
  var timestamp = newTime - startTime;

  myFrame.timestamp = timestamp;
  framesJSONobj.center = frame.interactionBox.center
  framesJSONobj.boundingBox = frame.interactionBox.size

  for (var i = 0; i < frame.hands.length; i++) {
    var hand = frame.hands[i];
    var currentHand = "";

    if ("right" == hand.type) {
      myFrame.right["palmPosition"] = hand.palmPosition
      currentHand = "right"
    }
    if ("left" == hand.type) {
      myFrame.left["palmPosition"] = hand.palmPosition
      currentHand = "left"
    };

    if (debugPrint) {
      handString += "Hand ID: " + hand.id + "<br />";
      handString += "Direction: " + vectorToString(hand.direction, 2) + "<br />";
      handString += "Palm normal: " + vectorToString(hand.palmNormal, 2) + "<br />";
      handString += "Palm position: " + vectorToString(hand.palmPosition) + " mm<br />";
      handString += "Palm velocity: " + vectorToString(hand.palmVelocity) + " mm/s<br />";
      handString += "Sphere center: " + vectorToString(hand.sphereCenter) + " mm<br />";
      handString += "Sphere radius: " + hand.sphereRadius.toFixed(1) + " mm<br />";
    }
    for (var j = 0; j < hand.pointables.length; j++) {
      var pointable = hand.pointables[j];

      var nameMap = ["thumb", "index", "middle", "ring", "pinky"];
      var fingerName = nameMap[pointable.type];

      myFrame[currentHand][fingerName] = []

      myFrame[currentHand][fingerName].push(pointable.carpPosition)
      myFrame[currentHand][fingerName].push(pointable.mcpPosition)
      myFrame[currentHand][fingerName].push(pointable.pipPosition)
      myFrame[currentHand][fingerName].push(pointable.dipPosition)
      if (debugPrint) {
        handString += "Pointable Name: " + fingerName + "<br />";
        handString += "Pointable ID: " + pointable.id + "<br />";
        handString += "Belongs to hand with ID: " + pointable.handId + "<br />";
        handString += "Length: " + pointable.length.toFixed(1) + " mm<br />";
        handString += "Width: "  + pointable.width.toFixed(1) + " mm<br />";
        handString += "Direction: " + vectorToString(pointable.direction, 2) + "<br />";
        handString += "Tip position: " + vectorToString(pointable.tipPosition) + " mm<br />";
        handString += "Tip velocity: " + vectorToString(pointable.tipVelocity) + " mm/s<br />";
        handString += "Carp position: " + vectorToString(pointable.carpPosition) + " mm<br />";
        handString += "MCP position: " + vectorToString(pointable.mcpPosition) + " mm<br />";
        handString += "DIP position: " + vectorToString(pointable.dipPosition) + " mm<br />";
        handString += "PIP position: " + vectorToString(pointable.pipPosition) + " mm<br />";
      }
    }
  }
  out.innerHTML = handString;
  framesJSONobj.frames.push(myFrame);
},

// hand callbacks are run once for each hand in the frame
hand: function(hand){
  out.innerHTML += "Hand: " + hand.id + ' &nbsp;roll: ' + Math.round(hand.roll() * Rad2Deg) + '°<br/>'
}

});

function onClick() {
   var txt = JSON.stringify(framesJSONobj);
   var blob = new Blob([txt], { type: "text/json;charset=utf-8"})
   var a = document.createElement("a");
   var url = URL.createObjectURL(blob);
   a.href = url;
   var d = new Date()
   var date = d.getTime()
   a.download = "leapMotion_" + date + ".json";
   document.body.appendChild(a);
   a.click();
   setTimeout(function() {
     document.body.removeChild(a);
     window.URL.revokeObjectURL(url);
   }, 0);
}

