# Source me. Selects JDK 21 and a TornadoVM 7.0.1 SDK; downloads the SDK into .tornado/ if none is given.
#   JAVA_HOME          a JDK 21 (default: $JAVA_HOME if it is 21, else sdkman's 21.0.2-open)
#   TORNADO_SDK        an existing TornadoVM 7.0.1 JDK 21 SDK (optional)
#   TORNADO_BACKEND    cuda (default) or opencl, for the downloaded SDK
TORNADO_VERSION=7.0.1
_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -z "${JAVA_HOME:-}" ]] || ! "$JAVA_HOME/bin/java" -version 2>&1 | grep -q '"21\.'; then
  JAVA_HOME="$HOME/.sdkman/candidates/java/21.0.2-open"
fi
if ! "$JAVA_HOME/bin/java" -version 2>&1 | grep -q '"21\.'; then
  echo "error: TornadoVM $TORNADO_VERSION (jdk21 build) needs JDK 21; set JAVA_HOME to one" >&2
  return 1 2>/dev/null || exit 1
fi
export JAVA_HOME

if [[ -z "${TORNADO_SDK:-}" ]]; then
  _backend="${TORNADO_BACKEND:-cuda}"
  _name="tornadovm-$TORNADO_VERSION-jdk21-$_backend-linux-amd64"
  TORNADO_SDK="$_root/.tornado/tornadovm-$TORNADO_VERSION-$_backend"
  if [[ ! -x "$TORNADO_SDK/bin/tornado" ]]; then
    echo "Downloading $_name ..." >&2
    mkdir -p "$_root/.tornado"
    curl -fL --progress-bar -o "$_root/.tornado/$_name.tar.gz" \
      "https://github.com/beehive-lab/TornadoVM/releases/download/v$TORNADO_VERSION/$_name.tar.gz"
    tar -xzf "$_root/.tornado/$_name.tar.gz" -C "$_root/.tornado"
    rm "$_root/.tornado/$_name.tar.gz"
    if [[ ! -x "$TORNADO_SDK/bin/tornado" ]]; then
      TORNADO_SDK="$(dirname "$(find "$_root/.tornado" -path '*/bin/tornado' | head -1)")/.."
      TORNADO_SDK="$(cd "$TORNADO_SDK" && pwd)"
    fi
  fi
fi
# the launcher prefers TORNADOVM_HOME over TORNADO_SDK: set both, or a stale TORNADOVM_HOME wins
export TORNADO_SDK TORNADOVM_HOME="$TORNADO_SDK"
export PATH="$TORNADO_SDK/bin:$JAVA_HOME/bin:$PATH"
