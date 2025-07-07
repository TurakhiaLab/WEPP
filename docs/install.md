WEPP offers multiple installation methods. Using a [Docker](https://docs.docker.com/engine/install/) is recommended to prevent any conflict with existing packages.

1. Docker image from DockerHub
2. Dockerfile 
3. Shell Commands 

!!!Note
     ⚠️The Docker image is currently built for the `linux/amd64` platform. While it can run on `arm64` systems (e.g., Apple Silicon or Linux aarch64) via emulation, this may lead to reduced performance.


## **Option-1: Install via DockerHub** <a name=dockerhub></a> 
The Docker image includes all dependencies required to run WEPP.

**Step 1:** Get the image from DockerHub 
```bash
docker pull pranavgangwar/wepp:latest
```
**Step 2:** Start and run Docker container. The command below will take you inside the docker container with WEPP already installed.
```bash
# -p <host_port>:<container_port> → Maps container port to a port on your host (Accessing Dashboard, NOT needed otherwise)
# Replace <host_port> with your desired local port (e.g., 80 or 8080)
# Use this command if your datasets can be downloaded from the Web
docker run -it -p 80:80 pranavgangwar/wepp:latest

# Use this command if your datasets are present in your current directory
docker run -it -p 80:80 -v "$PWD":/WEPP -w /WEPP pranavgangwar/wepp:latest
```
**Step 3:** Confirm proper working by running the following command. This should print WEPP's help menu.
```bash
snakemake test --cores 1 --use-conda
```

All set to try the [examples](quickstart.md#example).


## **Option-2: Install via Dockerfile** <a name=dockerfile></a> 
The Dockerfile contains all dependencies required to run WEPP.

**Step 1:** Clone the repository
```bash
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git 
cd WEPP
```
**Step 2:** Build a Docker Image
```bash
cd docker
docker build -t wepp . 
cd ..
```
**Step 3:** Start and run Docker container. The command below will take you inside the docker container with the view of the current directory.
```bash
# -p <host_port>:<container_port> → Maps container port to a port on your host (Accessing Dashboard, NOT needed otherwise)
# Replace <host_port> with your desired local port (e.g., 80 or 8080)
docker run -it -p 80:80 -v "$PWD":/workspace -w /workspace wepp
```

All set to try the [examples](quickstart.md#example).


## **Option-3: Install via Shell Commands (requires sudo access)** <a name=script></a>
Users without sudo access are advised to install WEPP via [Docker Image](#dockerhub).

**Step 1:** Clone the repository
```bash
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git
cd WEPP
```
**Step 2:** Install dependencies (might require sudo access)
WEPP depends on the following common system libraries, which are typically pre-installed on most development environments:
```text
- wget
- curl
- pip
- build-essential 
- python3-pandas
- pkg-config
- zip
- cmake 
- libtbb-dev
- libprotobuf-dev
- protobuf-compiler
- snakemake
- conda
- nodejs(v18+)
- nginx
```

For Ubuntu users with sudo access, if any of the required libraries are missing, you can install them with:
```bash
sudo apt-get update
sudo apt-get install -y wget pip curl python3-pip build-essential python3-pandas pkg-config zip cmake libtbb-dev libprotobuf-dev protobuf-compiler snakemake nginx
```

Note: WEPP expects the `python` command to be available. If your system only provides python3, you can optionally set up a symlink:
```bash
update-alternatives --install /usr/bin/python python /usr/bin/python3 1
```

If you do not have Node.js v18 or higher installed, follow these steps to install Node.js v22:
```bash
# Update and install prerequisites
apt-get install -y curl gnupg ca-certificates

# Add NodeSource Node.js 22 repo
curl -fsSL https://deb.nodesource.com/setup_22.x | bash -

# Install Node.js 22
apt-get install -y nodejs
```

```bash
# Install Yarn package manager globally
npm install -g yarn
# Install TaxoniumTools Python package
pip install taxoniumtools
```

If your system doesn't have Conda, you can install it with:
```bash
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/download/24.11.3-2/Miniforge3-24.11.3-2-Linux-x86_64.sh"
bash Miniforge3.sh -b -p "${HOME}/conda"

source "${HOME}/conda/etc/profile.d/conda.sh"
source "${HOME}/conda/etc/profile.d/mamba.sh"
```

All set to try the [examples](quickstart.md#example).