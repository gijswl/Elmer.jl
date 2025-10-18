import Conda

# Ensure pip interoperability is enabled, and install elmer-circuitbuilder from the pypi channel
Conda.pip_interop(true)
Conda.pip("install", "elmer-circuitbuilder")