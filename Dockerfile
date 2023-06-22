FROM ubuntu:18.04

WORKDIR /app/COMPAS

RUN apt-get update && apt-get install -y \
    g++ \
    libhdf5-serial-dev \
    libboost-all-dev \
    libgsl-dev \
    python3 \
    python3-pip \
    zip

# RUN pip install numpy awscli
RUN pip3 install numpy
RUN pip3 install pyyaml

# Copy only the source required to compile COMPAS
COPY src/ src/

RUN mkdir obj bin logs

ENV COMPAS_ROOT_DIR /app/COMPAS

# Compile COMPAS
RUN cd src && make -f Makefile.docker -j $(nproc)

# Run COMPAS
# CMD [ "python", "src/pythonSubmitDefault.py" ] 
