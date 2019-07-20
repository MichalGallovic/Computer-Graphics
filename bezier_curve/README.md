### Bezier curve

Draws bezier curve and redraws shape when new points are added or existing points moved. Uses [Casteljau's algorithm](https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm).

<p align="center"><img src="https://raw.githubusercontent.com/MichalGallovic/Computer-Graphics/master/bezier_curve/assets/preview.gif"></p>


#### Build and run

##### Linux

Install dependencies
```
apt-get update && apt-get install -y g++ freeglut3-dev libglm-dev
```

Build
```
g++ -std=c++11 main.cpp -lGL -lGLU -lGLUT -o main
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