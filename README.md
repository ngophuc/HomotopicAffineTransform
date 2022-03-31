# HomotopicAffineTransform
Source code of JMIV Paper: Homotopic Affine Transformations in the 2D Cartesian Grid

## Compilation

cd Sources

mkdir build

cd build

cmake .. -DDGtal_DIR=DGTAL_DIR [-DCMAKE_BUILD_TYPE=Release / Debug]

make -j4

## Execution
./HomotopicAffineTransform -i ../Samples/ball.pgm --a11 1.5 --a12 0.2 --a21 0.5 --a22 1.2

./HomotopicAffineTransform -i ../Samples/ball.pgm -o ../Results/ -s --a11 1.5 --a12 0.2 --a21 0.5 --a22 1.2

./HomotopicAffineTransform -i ../Samples/ball.pgm -o ../Results/ -s --a11 1.5 --a12 0.2 --a21 0.5 --a22 1.2 -m 0

./HomotopicAffineTransform -i ../Samples/ball.pgm -o ../Results/ -s --a11 1.5 --a12 0.2 --a21 0.5 --a22 1.2 -m 1

## Help
./HomotopicAffineTransform -h
Apply homotopic affine transforma on a given image.
 Example:
 	 HomotopicAffineTransform --input <PgmFileName> --outputDir <OutputDir> -s --a11 1.5 --a12 0.2 --a21 0.5 --a22 1.2 --tx 0 --ty 0 -m <0|1>

Usage: ./HomotopicAffineTransform [OPTIONS] 1

Positionals:
  1 TEXT:FILE REQUIRED                  Input file.

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Input file.
  -o,--output TEXT                      Output directory (default ./).
  -m,--model INT                        Transformation model: Majority vote (0, default), Gauss (1).
  -s,--save                             Save all steps (defaut no).
  --a11 FLOAT=1                         affine transform parameter (default 1.0)
  --a12 FLOAT=0                         affine transform parameter (default 0.0)
  --a21 FLOAT=0                         affine transform parameter (default 0.0)
  --a22 FLOAT=1                         affine transform parameter (default 1.0)
  --tx FLOAT=0                          X component of translation vector (default 0.0)
  --ty FLOAT=0                          Y component of translation vector (default 0.0)
