FROM dorowu/ubuntu-desktop-lxde-vnc:xenial

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    g++ freeglut3-dev libglm-dev \
    && rm -rf /var/lib/apt/lists/*

COPY . /root/Desktop
RUN chmod 755 /root/Desktop/docker-entrypoint.sh

ENTRYPOINT ["/root/Desktop/docker-entrypoint.sh"]