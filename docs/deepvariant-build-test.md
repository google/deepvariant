## Building DeepVariant from sources

DeepVariant comes with scripts to build it on Ubuntu 14 and 16, with
Ubuntu 16 recommended.  It can likely be built and run on other unix-based
systems with some minimal modifications to these scripts.  One way to get
access to a machine running Ubuntu is through a cloud computing platform
like Google Cloud Engine; see for example the instructions at the beginning
of the [DeepVariant Quick Start](deepvariant-quick-start.md).

To build and test DeepVariant, run the following commands:

```shell
./build-prereq.sh

./build_and_test.sh
```

At the end of the output of that last command, you should see a summary
message like "Executed 55 out of 55 tests: 55 tests pass." along with the
message "Target //deepvariant:binaries up-to-date:" followed by a list of
the just-built deepvariant binaries.

## Preparing a machine to run DeepVariant

The following command should be run on any machine on which
you wish run DeepVariant:

```shell
./run-prereq.sh
```

## settings.sh

Advanced users may want to edit the settings.sh file before building.
It contains options for configuring TensorFlow, CUDA, GPU usage, etc.
