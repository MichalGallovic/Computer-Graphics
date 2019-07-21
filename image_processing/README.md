### Image processing

Based on the Image of [Chichen Itza pyramid](https://en.wikipedia.org/wiki/Chichen_Itza) program tries to reconstruct its 3D image. 

It applies [Sobel filter](https://en.wikipedia.org/wiki/Sobel_operator) to detect edges of the object so it can remove background from the image and calculate slope of the pyramid walls.

<p align="center"><img src="https://raw.githubusercontent.com/MichalGallovic/Computer-Graphics/master/image_processing/assets/preview.png"></p>


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

Run binary with argument pointing to the image
```
./main "$(pwd)/assets/Wall.bmp"
```

##### Mac
Mac comes with OpenGL and GLUT preinstalled. To install `g++` you can do `brew install gcc`

Build
```
g++ -std=c++11 main.cpp -L/System/Library/Frameworks -framework GLUT -framework OpenGL -o main
```

Run binary with argument pointing to the image
```
./main "$(pwd)/assets/Wall.bmp"
```

#### Licence
[MIT license](https://opensource.org/licenses/MIT)