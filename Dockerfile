FROM ubuntu:22.04

WORKDIR /app/COMPAS

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    g++ \
    libhdf5-serial-dev \
    libboost-all-dev \
    libgsl-dev \
    python3 \
    python3-pip \
    zip \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install numpy pyyaml

# Copy only the source required to compile COMPAS
COPY src/ src/

RUN mkdir obj bin logs

ENV COMPAS_ROOT_DIR /app/COMPAS

# Compile COMPAS
RUN cd src && make DOCKER_BUILD=1 -j $(nproc)

# Run COMPAS
# CMD [ "python", "src/pythonSubmitDefault.py" ]