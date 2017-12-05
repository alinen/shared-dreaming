defmodule SharedDreamingWeb.DreamChannel do
  use Phoenix.Channel
  alias SharedDreamingWeb.Presence

  @limit 5

  # TODO: Handle reconnections

  def join("dream", _msg, socket) do
    Process.flag(:trap_exit, true)

    list = Presence.list(socket)
    size = list |> Map.size

    if size >= @limit do
      {:error, %{reason: "too many users connected"}}
    else
      send(self(), :after_join)

      {:ok, socket}
    end
  end

  def handle_info(:after_join, socket) do
    push socket, "presence_state", Presence.list(socket)
    {:ok, _} = Presence.track(socket, socket.assigns.token, %{
      online_at: inspect(System.system_time(:seconds))
    })

    {:noreply, socket}
  end

  def handle_in("dream:update", msg, socket) do
    broadcast! socket, "dream:update", %{user_id: socket.assigns.token, data: msg["data"]}

    {:noreply, socket}
  end
end
