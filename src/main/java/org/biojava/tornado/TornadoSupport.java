package org.biojava.tornado;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Decides whether the TornadoVM path may be used at all.
 * <p>
 * The GPU path is used only when the program runs under the TornadoVM runtime (the {@code tornado} launcher) and the
 * system property {@code biojava.tornado} is not {@code off}. Otherwise every accelerated class uses its exact CPU
 * path, so the same code runs everywhere.
 */
public final class TornadoSupport {

	private static final Logger logger = LoggerFactory.getLogger(TornadoSupport.class);

	/** System property: {@code off} disables the GPU path, {@code force} ignores the size thresholds. */
	public static final String PROPERTY = "biojava.tornado";

	private static volatile Boolean runtimeAvailable;

	private TornadoSupport() {
	}

	/** @return true if the GPU path is enabled and a TornadoVM runtime is present */
	public static boolean isEnabled() {
		if ("off".equalsIgnoreCase(System.getProperty(PROPERTY))) {
			return false;
		}
		return isRuntimeAvailable();
	}

	/** @return true if the size thresholds should be ignored (tests, benchmarks) */
	public static boolean isForced() {
		return "force".equalsIgnoreCase(System.getProperty(PROPERTY));
	}

	/**
	 * @param kernel name of the accelerated operation (e.g. "asa")
	 * @param work a size measure of the problem
	 * @param threshold the minimum size for which the GPU is worth it
	 * @return true if the GPU path should be taken for a problem of the given size
	 */
	public static boolean useGpu(String kernel, long work, long threshold) {
		return isEnabled() && !FAILED.contains(kernel) && (isForced() || work >= threshold);
	}

	private static final java.util.Set<String> FAILED = java.util.concurrent.ConcurrentHashMap.newKeySet();

	private static boolean isRuntimeAvailable() {
		Boolean available = runtimeAvailable;
		if (available == null) {
			try {
				// the runtime classes are only resolvable when launched through the TornadoVM SDK
				Class.forName("uk.ac.manchester.tornado.runtime.TornadoCoreRuntime");
				available = true;
			} catch (Throwable t) {
				logger.info("TornadoVM runtime not found, using the CPU paths");
				available = false;
			}
			runtimeAvailable = available;
		}
		return available;
	}

	/**
	 * Called when a GPU execution of the given operation failed (e.g. the device lacks FP64): logs it and disables
	 * the GPU path of that operation, and only that one, for the rest of the run.
	 */
	public static void disable(String kernel, Throwable cause) {
		if (FAILED.add(kernel)) {
			logger.warn("TornadoVM execution of {} failed, using its CPU path from now on: {}", kernel,
					cause.toString());
			logger.debug("TornadoVM failure", cause);
		}
	}
}
