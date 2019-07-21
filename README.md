### Computer graphics

Code from computer graphics course I've taken during my university studies. It includes basic opengl commands, image processing using convolutional filters, drawing bezier curves and my final project - game with physics, rigid bodies and collision detection.

Code was written in C++ using OpenGL/GLM and GLUT libraries.

For more information checkout README files of each subproject.

#### Build and run

Instructions on how to build and run each repository can be found in each folder.

#### Preview project using Docker

You can also preview the project using Docker:

```
# Change directory to the root of the project

cd Computer-Graphics

# Build docker image
docker build -t cg/preview:0.1.0 .

# Run docker image
docker run --rm -d -p 6080:80 -p 5900:5900 -e RESOLUTION=1920x1080 cg/preview:0.1.0

# Stop docker iamges
docker stop <container_id>
```

When docker is running, you can interact with Ubuntu and run projects using browser on `localhost:6080` or through VNC software on `localhost:5900`. All projects are pre-built on `Desktop`.

#### Bounce unlimited

3D Platform game Bounce Unlimited 💯

<p align="center"><img src="https://raw.githubusercontent.com/MichalGallovic/Computer-Graphics/master/bounce_unlimited/assets/preview.gif"></p>

Currently there is segmentation fault issue in this project and I did not have time to fix it. Though when built on mac without docker it works.

#### Image processing

Based on the Image of [Chichen Itza pyramid](https://en.wikipedia.org/wiki/Chichen_Itza) program tries to reconstruct its 3D image. 

<p align="center"><img src="https://raw.githubusercontent.com/MichalGallovic/Computer-Graphics/master/image_processing/assets/preview.png"></p>

#### Bezier curve

Draws bezier curve and redraws shape when new points are added or existing points moved. Uses [Casteljau's algorithm](https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm).

<p align="center"><img src="https://raw.githubusercontent.com/MichalGallovic/Computer-Graphics/master/bezier_curve/assets/preview.gif"></p>

#### Ball in cube & Fractals

Basic OpenGL drawing with camera movement and fractals drawn on the inner walls of cube.

<p align="center"><img src="https://raw.githubusercontent.com/MichalGallovic/Computer-Graphics/master/ball_in_a_cube/assets/preview.gif"></p>