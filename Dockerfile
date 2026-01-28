FROM ubuntu:focal

RUN apt-get update && apt-get install -y sudo

# Ensure tzdata installation (a dependency) is non-interactive.
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
RUN apt-get install -y software-properties-common wget build-essential cmake \
    git python3-dev xxd libgmp-dev libssl-dev libtbb-dev gdb
RUN apt-get install -y vim
RUN apt-get install -y --no-install-recommends libntl-dev
RUN apt-get install -y iproute2
RUN apt-get install -y valgrind

WORKDIR /home
RUN git clone https://github.com/nlohmann/json.git \
  && cd json \
  && mkdir build && cd build \
  && cmake -DCMAKE_BUILD_TYPE=Release .. \
  && make \
  && make install

WORKDIR /home
RUN wget https://archives.boost.io/release/1.82.0/source/boost_1_82_0.tar.bz2
RUN tar --bzip2 -xf boost_1_82_0.tar.bz2 && rm boost_1_82_0.tar.bz2
RUN cd /home/boost_1_82_0 && ./bootstrap.sh --with-python=/usr/bin/python3
RUN cd /home/boost_1_82_0 && ./b2 install

WORKDIR /home
RUN git clone https://github.com/emp-toolkit/emp-tool.git \
  && cd emp-tool \
  && cmake -DCMAKE_BUILD_TYPE=Release . \
  && make \
  && make install

CMD ["bash"]
