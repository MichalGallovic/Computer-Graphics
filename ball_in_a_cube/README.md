### Ball in Cube & Fractals

Basic OpenGL drawing with camera movement and fractals drawn on the inner walls of cube.

<p align="center"><img src="https://raw.githubusercontent.com/MichalGallovic/Computer-Graphics/master/ball_in_a_cube/assets/preview.gif"></p>

#### Controls
W S A D - Movement

c - Toggle camera mode between wall / follow mode

+/- Add/Remove square tiles on the floor

Mouse move - Look around when in "follow mode"

ESC - exit

#### Build and run

##### Linux

Install dependencies
```
apt-get update && apt-get install -y g++ freeglut3-dev libglm-dev
```

Build
```
g++ -std=c++11 main.cpp -lGL -lGLU -lglut -o main
```

Run binary
```
./main
```

##### Mac
Mac comes with OpenGL and GLUT preinstalled. To install `g++` you can do `brew install gcc`

Build
```
g++ -std=c++11 main.cpp -L/System/Library/Frameworks -framework GLUT -framework OpenGL -o main
```

Run binary
```
./main
```