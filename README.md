# TRANSFORMERS

![build status](https://gitlab.mrt.uni-karlsruhe.de/MRT/transformers/badges/master/build.svg)
![coverage report](https://gitlab.mrt.uni-karlsruhe.de/MRT/transformers/badges/master/coverage.svg)

Concatenating spatial transforms is tedious and error prone.
This project aims on making life easier while transforming.
Features are:
* Cheese cake (aka baseline): A coordinate transform that wraps transform types with coordinate system and tests at compile or runtime if those coordinate systems are concatenated correctly.
* Cherry 1: Those transforms can be registered in a topology graph from which you can extract shortest paths (aka topology path) of transformations at certain time instances.
* Cherry 2: Those topology paths can be given to the transformation graph which will return the corresponding transformation.
* Cherry 3: Support for temporaly static or dynamic transformations.

### [Doxygen documentation](http://MRT.pages.mrt.uni-karlsruhe.de/transformers/doxygen/index.html)
### [Coverage report](http://MRT.pages.mrt.uni-karlsruhe.de/transformers/coverage/index.html)

## Usage

See unittests (more detail coming up.)

## Todo
* Support for ceres-solver
* Optimization in transformation graph at compile time.
