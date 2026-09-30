#!/usr/bin/env bash
# Runs the JUnit tests under the TornadoVM runtime, so that the GPU paths are exercised
# (a plain `mvn test` has no TornadoVM runtime and only checks the CPU fallbacks).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/scripts/env.sh"
LAUNCHER="$ROOT/target/junit-platform-console-standalone.jar"
if [[ ! -f "$LAUNCHER" ]]; then
  (cd "$ROOT" && mvn -q -B dependency:copy -Dartifact=org.junit.platform:junit-platform-console-standalone:1.10.2 \
     -DoutputDirectory=target -Dmdep.stripVersion=true >/dev/null)
fi
(cd "$ROOT" && mvn -q -B test-compile)
exec "$ROOT/scripts/run.sh" -Dbiojava.tornado.expectGpu=true org.junit.platform.console.ConsoleLauncher \
  --disable-banner --details=tree --scan-classpath "$ROOT/target/test-classes" "$@"
