## Building DeepVariant from sources

DeepVariant comes with scripts to build it on Ubuntu 20.04. It can likely be
built and run on other unix-based systems with some minimal modifications to
these scripts. One way to get access to a machine running Ubuntu is through a
cloud computing platform like Google Cloud Engine.

First install the [Google Cloud SDK](https://cloud.google.com/sdk/downloads),
because we will need to use its `gsutil` command to fetch some dependencies.

The `build-prereq.sh` command below will install a number of system packages to
fulfill DeepVariant's prerequisites (using apt-get and pip, invoked via sudo).
This commands also downloads and builds TensorFlow and CLIF from source.

Run the following commands to install prerequisites, build the DeepVariant
programs, and then run tests.

```shell
./build-prereq.sh

./build_and_test.sh
```

At the end of the output of that last command, you should see a summary message
like "Executed 55 out of 55 tests: 55 tests pass." along with the message
"Target //deepvariant:binaries up-to-date:" followed by a list of the just-built
deepvariant binaries.

## Preparing a machine to run DeepVariant

The following command should be run on any machine on which you wish run
DeepVariant, since there are runtime dependencies, such as Python packages like
numpy and Tensorflow to be installed:

```shell
./run-prereq.sh
```

## Configuring the build

Advanced users may want to edit the settings.sh file before building. It
contains options for configuring TensorFlow, CUDA, GPU usage, etc.
