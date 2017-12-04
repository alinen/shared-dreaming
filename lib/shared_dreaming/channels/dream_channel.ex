defmodule SharedDreamingWeb.DreamChannel do
  use Phoenix.Channel

  # TODO: Handle reconnections
  # TODO: User tokens

  def join("dream", message, socket) do
    Process.flag(:trap_exit, true)
    send(self(), {:after_join, message})

    {:ok, socket}
  end

  def handle_info({:after_join, _msg}, socket) do
    broadcast! socket, "info", %{msg: "New user joining"}
    push socket, "join", %{status: "connected"}

    {:noreply, socket}
  end
end
