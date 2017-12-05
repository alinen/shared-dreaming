import {Socket, Presence} from "phoenix"

class Presences {
  constructor() {
    this.presences = {};
  }

  syncState(state) {
    this.presences = Presence.syncState(this.presences, state);
  }

  syncDiff(diff) {
    this.presences = Presence.syncDiff(this.presences, diff);
  }
}

class App {
  constructor({ token, channel, presences }) {
    this.token = token;
    this.channel = channel; // communication channel where to publish changes to other clients
    this.presences = presences; // Collection of current users connected

    // All internal variables should be declared here
  }

  setupChannel() {
    // Get the current state of the presences list after join
    this.channel.on("presence_state", state => this.presences.syncState(state))

    // Get the current state of the presences list after any change
    this.channel.on("presence_diff", diff => this.presences.syncDiff(diff))

    this.channel.on("dream:update", payload => this.update(payload))

    this.channel.join()
      .receive("ok", resp => { console.log("Joined successfully", resp) })
      .receive("error", resp => { console.log("Unable to join", resp) });
  }

  start() {
    console.log("Starting the app for user " + this.token);
    this.setupChannel();

    // Initialize the whole webGL setup

    // Let's simulate the tick. THIS IS JUST A DEMO
    window.setInterval(this.ping, 2000, this.channel);
  }

  update(payload) {
    if (payload.user_id != this.token) {
      console.log(`An update from ${payload.user_id}: ${payload.data}`);
    }
  }

  ping(channel) {
    channel.push("dream:update", { data: "PING" });
  }

};

let run = () => {
  let token = window.token;
  let socket = new Socket("/socket", {params: { token }});
  socket.connect();

  // Now that you are connected, you can join channels with a topic:
  let channel = socket.channel("dream", {});
  let presences = new Presences();

  let app = new App({ token, channel, presences });
  app.start();
};

window.onload = run;
