#!/bin/bash

# Compiles projects
cd /root/Desktop/ball_in_a_cube

g++ -std=c++11 main.cpp -lGL -lGLU -lglut -o main

cd /root/Desktop/bezier_curve

g++ -std=c++11 main.cpp -lGL -lGLU -lglut -o main

cd /root/Desktop/image_processing

g++ -std=c++11 main.cpp -lGL -lGLU -lglut -o main

cd /root/Desktop/bounce_unlimited

g++ -std=c++11 main.cpp -lGL -lGLU -lglut -o main

/startup.sh