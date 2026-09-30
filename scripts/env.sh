# Source me. Pins the JDK and TornadoVM SDK used by the scripts.
# Note: the `tornado` launcher prefers TORNADOVM_HOME over TORNADO_SDK, so both are set.
export JAVA_HOME="${JAVA_HOME_21:-$HOME/.sdkman/candidates/java/21.0.2-open}"
export TORNADO_SDK="${TORNADO_SDK_OVERRIDE:-$HOME/.sdkman/candidates/tornadovm/7.0.1-jdk21-cuda}"
export TORNADOVM_HOME="$TORNADO_SDK"
export PATH="$TORNADO_SDK/bin:$JAVA_HOME/bin:$PATH"
