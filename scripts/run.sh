#!/usr/bin/env bash
# Runs a main class of this project under the TornadoVM runtime.
#   scripts/run.sh [-Dprop=value ...] <main-class> [args...]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/scripts/env.sh"
if [[ ! -f "$ROOT/target/cp.txt" || "$ROOT/pom.xml" -nt "$ROOT/target/cp.txt" ]]; then
  (cd "$ROOT" && mvn -q -B dependency:build-classpath -Dmdep.outputFile=target/cp.txt -Dmdep.includeScope=test \
     -Dmdep.excludeArtifactIds=tornado-api >/dev/null)
fi
JVM=()
while [[ $# -gt 0 && "$1" == -* ]]; do JVM+=("$1"); shift; done
CP="$ROOT/target/classes:$ROOT/target/test-classes:$(cat "$ROOT/target/cp.txt")"
exec tornado --jvm="-Xmx24g -Dtornado.device.memory=20GB ${JVM[*]:-}" -cp "$CP" "$@"
