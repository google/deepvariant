# Build local hap.py v0.3.9

## Context

`pkrusche/hap.py:v0.3.9` no longer works.

For reproducibility, here is how to build one that you can use, which should
still be v0.3.9.

## To build

1.  **Create a file named `Dockerfile`** with the content from
    [tools/Dockerfile.happy_v0.3.9](../tools/Dockerfile.happy_v0.3.9).

2. **Make sure the file is in an empty directory:**

```
mkdir docker_build_happy
cp /path/to/Dockerfile.happy_v0.3.9 docker_build_happy/Dockerfile
cd docker_build_happy
```

3.  **Build the Docker image** using the following command:

```bash
docker build -t local_happy:v0.3.9 .
```

After the build process completes, the image `local_happy:v0.3.9` will be
available in the local Docker image repository.

This takes quite a while to build. You might want to save it somewhere to reuse
later.

4. **Run the binary**

```bash
docker run local_happy:v0.3.9 /usr/local/bin/som.py
```
