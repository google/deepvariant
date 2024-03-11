# Getting Started with GCP

DeepVariant doesn't require GCP, but if you want to use it, these are some
instructions that we found to be useful when getting started.

## Set up a Google Cloud account

To get started using DeepVariant on Google Cloud Platform (GCP), you first need
to set up an account and a project to contain your cloud resources.

*   If you do not have an account yet, you should create one at
    [cloud.google.com](https://cloud.google.com). You should then [enable
    billing for your
    account](https://support.google.com/cloud/answer/6288653?hl=en) but note
    that if your account is new, [you receive $300 of free
    credit](https://cloud.google.com/free/). Once your cloud account is set up,
    you should be able to log in to the [Cloud
    Console](https://console.cloud.google.com) to view or administer your cloud
    resources.

*   From the Cloud Console, [set up a
    project](https://cloud.google.com/resource-manager/docs/creating-managing-projects)
    to house all of the cloud resources (storage, compute, services) that you
    will associate with your use of DeepVariant. For example, if your
    organization is AcmeCorp, you might call your project
    `acmecorp-deepvariant`.

*   Finally, please visit the ["Compute Engine" page on Cloud
    Console](https://console.cloud.google.com/compute). You don't need to create
    Compute Engine instances at this time, but simply visiting this page will
    initialize your compute engine "service account" so that we can authorize
    it.

(As you progress in your use of Google Cloud Platform, you will likely find it
useful to create a [Cloud
Organization](https://cloud.google.com/resource-manager/docs/creating-managing-organization)
to house your projects. Here are some [best
practices](https://cloud.google.com/docs/enterprise/best-practices-for-enterprise-organizations)
for organizating cloud projects for an enterprise.)

## Install the Google Cloud SDK

The Google Cloud SDK comes with two very useful command line utilities that you
can use on your local workstation---`gcloud`, which lets you administer your
cloud resources, and `gsutil`, which lets you manage and transfer data to Google
Cloud Storage buckets. We will make use of these tools in the following
instructions. To install the Cloud SDK, [follow the installation instructions
here](https://cloud.google.com/sdk/downloads).

The final step in the installation process (`gcloud init`) will have you
authenticate via your web browser and select a default [zone and
region](https://cloud.google.com/compute/docs/regions-zones/regions-zones) for
your cloud resources, which you can choose based on your location and regional
hardware availability.

NOTE: Not all zones are equipped with GPUs, so if you want to use GPUs for your
project, please take note of the availability listing
[here](https://cloud.google.com/compute/docs/gpus/).

To verify that the installation and authentication succeeded, run

```shell
gcloud auth list
```

and verify that your account email address is printed.

## Starting a Compute Engine instance

A simple way to access compute on GCP is Google Compute Engine. Compute Engine
instances can be sized to meet computational and storage needs for your project.

Before we get started, [ensure you have adequate quota
provisioned](https://cloud.google.com/compute/quotas) so that you can get all
the CPUs/GPUs that you need. To start with, you might want to request quota for
64 CPUs and 2 GPUs in your zone.

DeepVariant can make use of multiple CPU cores and (currently, a single) GPU
device. For this "quick start" guide, let's allocate an 8-core non-preemptible
instance in your default zone with a single GPU, running Ubuntu 20.04, with a
disk of reasonable size for modest work with genomic data. From our local
command line, we do:

```shell
gcloud beta compute instances create "${USER}-deepvariant-quickstart" \
  --scopes "compute-rw,storage-full,cloud-platform"  \
  --image-family ubuntu-2204-lts --image-project ubuntu-os-cloud \
  --machine-type n1-standard-8  \
  --boot-disk-size=200GB \
  --zone us-west1-b \
  --accelerator type=nvidia-tesla-k80,count=1 --maintenance-policy TERMINATE --restart-on-failure
```

NOTE: To create an instance *without GPU*, simply omit the last line from the
command.

Check that the instance has been created and started:

```shell
gcloud compute instances list
```

which should produce output like:

```
NAME                    ZONE        MACHINE_TYPE    PREEMPTIBLE   INTERNAL_IP  EXTERNAL_IP     STATUS
[USER]-deepvariant-quickstart  us-west1-b  n1-standard-8                 10.138.0.4   35.185.203.59   RUNNING
```

Then connect to your instance via SSH:

```shell
gcloud compute ssh --zone us-west1-b "${USER}-deepvariant-quickstart"
```

You should land at a shell prompt in your new instance!

NOTE: All of these steps can also be completed from the Cloud Console, if you
prefer. Consult [this
guide](https://cloud.google.com/compute/docs/quickstart-linux), but be sure to
choose Ubuntu 20.04 as your image, as DeepVariant has not been tested on other
Linux distributions.

For more information about getting started with Compute Engine, see:

*   [Compute Engine instance creation in `gcloud`
    manual](https://cloud.google.com/sdk/gcloud/reference/compute/instances/create)
*   [Reference to machine
    sizes/types](https://cloud.google.com/compute/docs/machine-types)
