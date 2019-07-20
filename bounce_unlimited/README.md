### Bounce unlimited

3D Platform game Bounce Unlimited 💯

<p align="center"><img src="https://raw.githubusercontent.com/MichalGallovic/Computer-Graphics/master/bounce_unlimited/assets/preview.gif"></p>

Collect points running over `red squares`. 
Every 5 points enemies level up. 
You only have 3 lives.
If you collide with enemy you loose 1 life.

Ball controls:
W S A D - movement
Space - jump
P - Super POWER!

Game controls:
R - start/reset game
I - instructions
ESC - quit game

Physics engine, Collision detection and vector operation is based on book [Game Physics Engine Development](https://www.goodreads.com/book/show/1501484.Game_Physics_Engine_Development_With_CDROM_)

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