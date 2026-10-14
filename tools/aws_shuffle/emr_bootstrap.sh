#!/bin/bash
# EMR on EC2 bootstrap action: install the (small) Python deps of
# shuffle_tfrecords_spark.py on every node. No Beam, no TensorFlow needed.
set -euo pipefail
sudo python3 -m pip install --no-cache-dir google-crc32c psutil boto3
