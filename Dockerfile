FROM dorowu/ubuntu-desktop-lxde-vnc:xenial

RUN apt-get update \
    && apt-get install -y g++ freeglut3-dev libglm-dev
