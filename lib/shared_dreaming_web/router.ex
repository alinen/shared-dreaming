defmodule SharedDreamingWeb.Router do
  use SharedDreamingWeb, :router

  pipeline :browser do
    plug :accepts, ["html"]
    plug :fetch_session
    plug :fetch_flash
    plug :protect_from_forgery
    plug :put_secure_browser_headers
    plug :auth_user
    plug :put_user_token
  end

  pipeline :api do
    plug :accepts, ["json"]
  end

  scope "/", SharedDreamingWeb do
    pipe_through :browser # Use the default browser stack

    get "/", PageController, :index
  end

  # Other scopes may use custom stacks.
  # scope "/api", SharedDreamingWeb do
  #   pipe_through :api
  # end

  defp auth_user(conn, _) do
    current_user = Plug.Conn.get_session(conn, :current_user)
    if current_user do
      conn
    else
      length = 16
      id = :crypto.strong_rand_bytes(length) |> Base.encode64 |> binary_part(0, length)
      Plug.Conn.put_session(conn, :current_user, id)
    end
  end

  defp put_user_token(conn, _) do
    if current_user = Plug.Conn.get_session(conn, :current_user) do
      assign(conn, :token, current_user)
    else
      conn
    end
  end
end
