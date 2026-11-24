# containers

Recipes for container images of this project's software, built from `pixi.lock`, so the container and a local pixi
environment run exactly the same packages:

- `Dockerfile` — a pixi environment (`ENV_NAME`, default `default`) copied into a small Ubuntu image; the project
  folder is mounted at `/work`, not copied in. Build and run commands are at the top of the file.
- `apptainer.def` — the same image for clusters that use Apptainer.

Pixi alone is enough to reproduce the results; the images are an extra for systems without pixi. Built images
(`*.tar.gz`, `*.sif`) are not committed.
